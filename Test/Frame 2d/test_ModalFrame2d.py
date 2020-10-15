# Three dimensional Frame: Eigenvalue Analysis & Effective Modal Mass Participation Ratios
# REFERENCES:
# used in verification by SAP2000 and SeismoStruct:
# SAP2000 Integrated Finite Element Analysis and Design of Structures, Verification Manual,
# Computers and Structures, 2009. Example 1-022.
# SeismoStruct, Verification Report, 2020. Example 11.

def test_ModalAnalysis_2DFrame():

    print("=====================================================================")
    print("ModalFrame2d: ModalAnalysis of 2D Frame Structure, using nodal masses")        

    import openseespy.opensees as op
    from openseespy.postprocessing.ModalAnalysis import ModalAnalysis
    from math import asin, sqrt

    # number of eigenvalues to calculate
    numEigen = 8

    op.wipe()
    op.model('Basic', '-ndm', 2, '-ndf', 3)

    # properties

    #    units kip, ft

    numBay = 2
    numFloor = 7

    bayWidth = 360.0
    storyHeights = [162.0, 162.0, 156.0, 156.0, 156.0, 156.0, 156.0]

    E = 29500.0
    massX = 0.49
    M = 0.
    coordTransf = "Linear"  # Linear, PDelta, Corotational
    massType = "-cMass"  # -lMass, -cMass

    beams = ['W24X160', 'W24X160', 'W24X130', 'W24X130', 'W24X110', 'W24X110', 'W24X110']
    eColumn = ['W14X246', 'W14X246', 'W14X246', 'W14X211', 'W14X211', 'W14X176', 'W14X176']
    iColumn = ['W14X287', 'W14X287', 'W14X287', 'W14X246', 'W14X246', 'W14X211', 'W14X211']
    columns = [eColumn, iColumn, eColumn]

    WSection = {
        'W14X176': [51.7, 2150.],
        'W14X211': [62.1, 2670.],
        'W14X246': [72.3, 3230.],
        'W14X287': [84.4, 3910.],
        'W24X110': [32.5, 3330.],
        'W24X130': [38.3, 4020.],
        'W24X160': [47.1, 5120.]
    }

    nodeTag = 1


    # procedure to read
    def ElasticBeamColumn(eleTag, iNode, jNode, sectType, E, transfTag, M, massType):
        found = 0

        prop = WSection[sectType]

        A = prop[0]
        I = prop[1]
        op.element('elasticBeamColumn', eleTag, iNode, jNode, A, E, I, transfTag, '-mass', M, massType)


    # add the nodes
    #  - floor at a time
    yLoc = 0.
    for j in range(0, numFloor + 1):

        xLoc = 0.
        for i in range(0, numBay + 1):
            op.node(nodeTag, xLoc, yLoc)
            xLoc += bayWidth
            nodeTag += 1

        if j < numFloor:
            storyHeight = storyHeights[j]

        yLoc += storyHeight

    # fix first floor
    op.fix(1, 1, 1, 1)
    op.fix(2, 1, 1, 1)
    op.fix(3, 1, 1, 1)

    # rigid floor constraint & masses
    nodeTagR = 5
    nodeTag = 4
    for j in range(1, numFloor + 1):
        for i in range(0, numBay + 1):

            if nodeTag != nodeTagR:
                op.equalDOF(nodeTagR, nodeTag, 1)
            else:
                op.mass(nodeTagR, massX, 0.0, 0.0)

            nodeTag += 1

        nodeTagR += numBay + 1

    # add the columns
    # add column element
    op.geomTransf(coordTransf, 1)
    eleTag = 1
    for j in range(0, numBay + 1):

        end1 = j + 1
        end2 = end1 + numBay + 1
        thisColumn = columns[j]

        for i in range(0, numFloor):
            secType = thisColumn[i]
            ElasticBeamColumn(eleTag, end1, end2, secType, E, 1, M, massType)
            end1 = end2
            end2 += numBay + 1
            eleTag += 1

    # add beam elements
    for j in range(1, numFloor + 1):
        end1 = (numBay + 1) * j + 1
        end2 = end1 + 1
        secType = beams[j - 1]
        for i in range(0, numBay):
            ElasticBeamColumn(eleTag, end1, end2, secType, E, 1, M, massType)
            end1 = end2
            end2 = end1 + 1
            eleTag += 1

    # calculate eigenvalues & print results
    numEigen = 7
    T, Mratios, Mfactors, Mtots = ModalAnalysis(numEigen, outname='OpenSeespy', pflag=0)

    # determine PASS/FAILURE of test
    ok = 0

    #               SAP2000   SeismoStruct
    comparisonResults = [[1.2732, 0.4313, 0.2420, 0.1602, 0.1190, 0.0951, 0.0795],
                         [1.2732, 0.4313, 0.2420, 0.1602, 0.1190, 0.0951, 0.0795]]
    print("\n\nComparisons of Periods [sec]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))

    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.5f}{:>15.4f}{:>15.4f}'.format(i + 1, T[i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(T[i] - resultOther) > 9.99e-5:
            ok - 1

    #               SAP2000   SeismoStruct
    comparisonResults = [[79.96,11.34,4.18,2.12,1.41,0.68,0.31],
                         [79.962562, 11.336245, 4.181011, 2.115034, 1.414566, 0.679965, 0.310617]]

    print("\n\nComparisons for Effective Modal Mass Participating in U\u2081 [%]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))

    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.5f}{:>15.4f}{:>15.4f}'.format(i + 1, Mratios[1][i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(Mratios[1][i] - resultOther) > 9.99e-5:
            ok - 1

    assert ok == 0