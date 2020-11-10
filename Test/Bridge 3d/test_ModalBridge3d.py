# Three Dimensional Bridge structure: Eigenvalue Analysis & Effective Modal 
# Mass Participation Ratios
# This example is a made up example just to illustrate that ModalAnalysis()
# can be used even if distributed element masses are used

def test_ModalAnalysis_3DBridge():

    print("=====================================================================")
    print("ModalBridge3d: ModalAnalysis of 3D Bridge Structure, using element masses")        

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
    # First Moment of Area
    m3 = m**3
    cm3 = cm**3
    mm3 = mm**3
    inch3 = inch**3
    # Second Moment of Area
    m4 = m**4
    cm4 = cm**4
    mm4 = mm**4
    inch4 = inch**4
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
    
    #  Elastic Material Properties
    v = 0.2
    E = 30*GPa
    G = E/(2*(1+v))
    
    # Deck Properties
    DeckA = 7.85**m2
    DeckJ = 16.58*m4
    DeckIy = 6.76*m4
    DeckIz = 33.6*m4
    DeckMass = 20
    
    # start nodes
    X0 = [0.0, 40*m, 90*m, 140*m]
    Y0 = [0.0, 0.0, 0.0, 0.0]
    Z0 = [0.0, 0.0, 0.0, 0.0] 
    # end nodes
    X1 = [40.0, 90*m, 140*m, 180*m]
    Y1 = [0.0, 0.0, 0.0, 0.0]
    Z1 = [0.0, 0.0, 0.0, 0.0]
    
    # Pier Properties
    PierA = 3.14*m2
    PierJ = 1.57*m4
    PierIz = 0.785*m4
    PierIy = 0.785*m4
    PierMass = 8
    
    PierH  = [12*m,4*m,8*m]
    
    # number of eigenvalues to calculate
    numEigen = 8
    
    op.wipe()
    op.model('basic', '-ndm', 3, '-ndf', 6)
    masstype = 'Lumped' # 'Consistent', 'Lumped'
    
    DeckNodes = []
    DeckTransfTag = 1
    PierTransfTag = 2
    
    op.geomTransf('Linear', DeckTransfTag, 0, 0, 1)
    op.geomTransf('Linear', PierTransfTag, 0, 1, 0)
    
    for i in range(len(X0)):
        # Deck nodes and elements
        StartNode_d = int(str(i+1)+'00')
        EndNode_d = int(str(i+2)+'00')
        if i == 0:
            op.node(StartNode_d, X0[i], Y0[i], Z0[i])
            op.fix(StartNode_d,0,1,1,1,0,1)
            
        op.node(EndNode_d, X1[i], Y1[i], Z1[i])
        if i == len(X0)-1:
            op.fix(EndNode_d,0,1,1,1,0,1)
        
        eleNodes = [StartNode_d, EndNode_d]
        eleTag = int(str(i+1)+'00')
        
        if masstype == 'Lumped':
            op.element('elasticBeamColumn', eleTag, *eleNodes, DeckA, E, G, DeckJ, DeckIy, DeckIz, DeckTransfTag,'-mass',DeckMass)
        elif masstype == 'Consistent':
            op.element('elasticBeamColumn', eleTag, *eleNodes, DeckA, E, G, DeckJ, DeckIy, DeckIz, DeckTransfTag,'-mass',DeckMass,'-cMass')
            
        # Pier nodes and elements        
        if i!=len(X0)-1:
            StartNode_p = int(str(i+1)+'01')
            EndNode_p = EndNode_d
            op.node(StartNode_p, X1[i], Y1[i], Z1[i]-PierH[i])
            op.fix(StartNode_p,1,1,1,1,1,1)
            eleNodes = [StartNode_p, EndNode_p]
            eleTag = int(str(i+1)+'01')
            if masstype == 'Lumped':
                op.element('elasticBeamColumn', eleTag, *eleNodes, PierA, E, G, PierJ, PierIz, PierIy, PierTransfTag,'-mass',PierMass)
            elif masstype == 'Consistent':
                op.element('elasticBeamColumn', eleTag, *eleNodes, PierA, E, G, PierJ, PierIz, PierIy, PierTransfTag,'-mass',PierMass,'-cMass')
    
    T, Mratios, Mfactors, Mtots = ModalAnalysis(numEigen, outname='OpenSeespy', pflag=1)
    
    # determine PASS/FAILURE of test
    ok = 0
    
    #               SAP2000   SeismoStruct
    comparisonResults = [[0.348192,0.271303,0.211800,0.141027,0.107349,0.068943,0.055869,0.055094],
                         [0.34819207,0.27130325,0.21179972,0.14102723,0.10734891,0.06894083,0.05585407,0.05509375]]
    print("\n\nComparisons of Periods [sec]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.5f}{:>15.4f}{:>15.4f}'.format(i + 1, T[i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(T[i] - resultOther) > 1e-5:
            ok - 1
    
    comparisonResults = [[0.000000,0.000000,98.865000,0.000000,0.160000,0.000957,0.023000,0.923000],
                         [0.000000,0.000000,98.865005,0.000000,0.160053,0.000957,0.024196,0.922046]]
    print("\n\nComparisons for Effective Modal Mass Participating in U\u2081 [%]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.4f}{:>15.4f}{:>15.4f}'.format(i + 1, Mratios[1][i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(Mratios[1][i] - resultOther) > 1e-3:
            ok - 1
            
    comparisonResults = [[33.129000,49.987000,0.000000,16.884000,0.000000,0.000000,0.000000,0.000000],
                         [33.128723,49.986901,0.000000,16.884376,0.000000,0.000000,0.000000,0.000000]]
    print("\n\nComparisons for Effective Modal Mass Participating in U\u2082 [%]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.4f}{:>15.4f}{:>15.4f}'.format(i + 1, Mratios[2][i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(Mratios[2][i] - resultOther) > 1e-3:
            ok - 1
            
    comparisonResults = [[0.000000,0.000000,0.000111,0.000000,0.000008,32.735000,31.414000,0.972000],
                         [0.000000,0.000000,0.000111,0.000000,0.000008,32.753245,31.603705,1.007302]]
    print("\n\nComparisons for Effective Modal Mass Participating in U\u2083 [%]:")
    print('{:>10}{:>15}{:>15}{:>15}'.format('Mode', 'OpenSees', 'SAP2000', 'SeismoStruct'))
    
    # formatString {%10s%15.5f%15.4f%15.4f}
    for i in range(0, numEigen):
        print('{:>10}{:>15.4f}{:>15.4f}{:>15.4f}'.format(i + 1, Mratios[3][i], comparisonResults[0][i], comparisonResults[1][i]))
        resultOther = comparisonResults[0][i]
        if abs(Mratios[3][i] - resultOther) > 1e-3:
            ok - 1

    assert ok == 0