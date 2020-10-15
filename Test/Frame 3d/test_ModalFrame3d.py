# Three dimensional Frame: Eigenvalue Analysis & Effective Modal Mass Participation Ratios
# REFERENCES:
# used in verification by SAP2000 and SeismoStruct:
# SAP2000 Integrated Finite Element Analysis and Design of Structures, Verification Manual,
# Computers and Structures, 2009. Example 1-024.
# SeismoStruct, Verification Report, 2020. Example 12.

def test_ModalAnalysis_3DFrame():
    
    print("=====================================================================")
    print("ModalFrame3d: ModalAnalysis of 3D Frame Structure, using nodal masses")    
    
    import openseespy.opensees as op
    from openseespy.postprocessing.ModalAnalysis import ModalAnalysis
    import numpy as np
    
    # Basic Units
    m = 1.0
    kN = 1.0
    sec = 1.0
    LunitTXT = 'meter'
    FunitTXT = 'kN'
    TunitTXT = 'sec'
    # Constants
    pi = np.pi
    g = 9.81*m/sec**2
    # Length
    mm = m/1000.0
    cm = m/100.0
    inch = 25.4*mm
    ft = 12.0*inch
    # Area
    m2 = m**2
    cm2 = cm**2
    mm2 = mm**2
    inch2 = inch**2
    ft2 = ft**2
    # First Moment of Area
    m3 = m**3
    cm3 = cm**3
    mm3 = mm**3
    inch3 = inch**3
    ft3 = ft**3
    # Second Moment of Area
    m4 = m**4
    cm4 = cm**4
    mm4 = mm**4
    inch4 = inch**4
    ft4 = ft**4
    # Force
    N = kN/1000.0
    kip = kN*4.448221615
    # Moment
    kNm = kN*m
    # Stress (kN/m2 or kPa)
    Pa = N/(m2)
    kPa = Pa*1.0e3
    MPa = Pa*1.0e6
    GPa = Pa*1.0e9
    ksi = 6.8947573*MPa
    psi = 1e-3*ksi
    # Angles
    degrees = pi/180.0
    
    op.wipe()
    op.model('basic', '-ndm', 3, '-ndf', 6)
    
    # number of eigenvalues to calculate
    numEigen = 4
    
    # Frame grid
    Xs = [0,35*ft,70*ft]
    Ys = [0,25*ft,50*ft]
    Zs = [0,13*ft,26*ft]
    
    # Center of mass at each floor
    Xcm = 38*ft
    Ycm = 27*ft
    
    # Lumped floor masses
    massX = 6.2112*kip*sec**2/ft
    massY = 6.2112*kip*sec**2/ft
    
    # Distributed Loading
    wload = -10*kip/ft
    
    # Beam section properties
    v = 0.2
    Eb = 500000*kip/ft2
    Gb = Eb/(2*(1+v))
    Ab = 5*ft2
    Iyb = 2.61*ft4
    Izb = 1.67*ft4
    Jb = 0
    
    # Column section properties
    v = 0.2
    Ec = 350000*kip/ft2
    Gc = Ec/(2*(1+v))
    Ac = 4*ft2
    Izc = 1.25*ft4
    Iyc = 1.25*ft4
    Jc = 0
    
    # Define Transformation Tags
    ColTransf = 1
    BeamXTransf = 2
    BeamYTransf = 3
    op.geomTransf('Linear', ColTransf, 0, 1, 0)
    op.geomTransf('Linear', BeamXTransf, 0, 0, 1)
    op.geomTransf('Linear', BeamYTransf, 0, 0, 1)
    
    # define NODAL COORDINATES
    storey = 0
    
    for k in range(len(Zs)):
        
        no = 1
        constrained = []
        
        for i in range(len(Xs)):
            for j in range(len(Ys)):
                nodeID = int(str(no)+'00'+str(storey))
                op.node(nodeID,Xs[i],Ys[j],Zs[k])
                if k == 0:
                    op.fix(nodeID,1,1,1,1,1,1)     
                else:
                    constrained.append(nodeID)
                no += 1
                
        if k != 0: # Center of mass
            nodeID = int(str(no)+'00'+str(storey))
            op.node(nodeID,Xcm,Ycm,Zs[k])
            op.mass(nodeID,massX,massY,0,0,0,0)
            # constrain other dofs that don't belong to rigid diaphragm control
            op.fix(nodeID,0,0,1,1,1,0)		
            op.rigidDiaphragm(3,nodeID,*constrained)
        storey += 1
    
    # Define Columns
    colTag = '00'
    for no in range(1,(len(Xs))*(len(Ys))+1):
        for storey in range(1,len(Zs)):
            nodeI = int(str(no)+'00'+str(storey-1))
            nodeJ = int(str(no)+'00'+str(storey))     
            eleTag = int(str(no)+colTag+str(storey))
            op.element('elasticBeamColumn',eleTag,nodeI,nodeJ,Ac, Ec, Gc, Jc, Iyc, Izc, ColTransf)
    
    beamEles = []
    # Define Beams in X axis
    beamXtag = '01'    
    for storey in range(1,len(Zs)):
        no = 1
        for i in range(len(Xs)-1):
            for j in range(1,len(Ys)+1):
                nodeI = int(str(i*len(Ys)+j)+'00'+str(storey))
                nodeJ = int(str((i+1)*len(Ys)+j)+'00'+str(storey))
                eleTag = int(str(no)+beamXtag+str(storey))
                beamEles.append(eleTag)
                op.element('elasticBeamColumn',eleTag,nodeI,nodeJ,Ab, Eb, Gb, Jb, Iyb, Izb, BeamXTransf)
                no+=1
    
    # Define Beams in Y axis
    beamYtag = '02'    
    for storey in range(1,len(Zs)):
        no = 1
        for i in range(len(Xs)):
            for j in range(1,len(Ys)):
                nodeI = int(str(i*len(Ys)+j)+'00'+str(storey))
                nodeJ = int(str(i*len(Ys)+(j+1))+'00'+str(storey))
                eleTag = int(str(no)+beamYtag+str(storey))
                beamEles.append(eleTag)
                op.element('elasticBeamColumn',eleTag,nodeI,nodeJ,Ab, Eb, Gb, Jb, Iyb, Izb, BeamYTransf)
                no+=1
    
    T, Mratios, Mfactors, Mtots = ModalAnalysis(numEigen, outname='OpenSeespy', pflag=1)
    
    # determine PASS/FAILURE of test
    ok = 0
    
    #               SAP2000   SeismoStruct
    comparisonResults = [[0.227062,0.215633,0.073345,0.072005],
                         [0.22706191,0.21563345,0.07334548,0.07200536]]
    print("\n\nComparisons of Periods [sec]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.5f}{:>15.4f}{:>15.4f}'.format(i + 1, T[i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(T[i] - resultOther) > 1e-5:
            ok - 1
    
    comparisonResults = [[90.258,0.19,9.388,0.164],
                         [90.257941,0.189952,9.38773,0.164377]]
    print("\n\nComparisons for Effective Modal Mass Participating in U\u2081 [%]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.4f}{:>15.4f}{:>15.4f}'.format(i + 1, Mratios[1][i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(Mratios[1][i] - resultOther) > 1e-3:
            ok - 1
            
    comparisonResults = [[0.192,91.046,0.151,8.612],
                         [0.191706,91.046194,0.150526,8.611574]]
    print("\n\nComparisons for Effective Modal Mass Participating in U\u2082 [%]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.4f}{:>15.4f}{:>15.4f}'.format(i + 1, Mratios[2][i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(Mratios[2][i] - resultOther) > 1e-3:
            ok - 1

    assert ok == 0