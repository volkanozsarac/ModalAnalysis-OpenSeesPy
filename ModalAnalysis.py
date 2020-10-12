def ModalAnalysis(numEigen, outname=None, pflag=1):
    """
    Details
    -------
        This script will perform a modal analysis on OpenSeespy model.
    
    Information
    -----------
        Author: Volkan Ozsarac, Earthquake Engineering PhD Candidate
        Affiliation: University School for Advanced Studies IUSS Pavia
        e-mail: volkanozsarac@iusspavia.it
    
    References
    ----------
        Chopra, A.K. 2012. Dynamics of Structures: Theory and 
        Applications to Earthquake Engineering, Prentice Hall.

    Notes
    -----
        Total mass is obtained by summing the masses assigned to the
        unrestrained dofs.
        
    Parameters
    ----------
    numEigen : int
        Number of eigenvalues to calculate.
    outname  : str, optional (The default is None)
        if not None, the modal properties for the first numEigen modes
        will be printed into outname.csv.
    pflag    : int (1 or 0)
        flag to print output information on screen

    Returns
    -------
    T        : numpy.ndarray
        Period array for the first numEigen modes.
    Mratios  : dictionary
        Effective modal mass participation ratios for the first numEigen modes.
    Mfactors : dictionary
        Modal particpation factors for the first numEigen modes.
    Mtots    : dictionary
        Total mass of the structure

    """
    import numpy as np
    import openseespy.opensees as op
    import sys

    op.wipeAnalysis()
    op.system('FullGeneral')
    op.analysis('Transient')

    # Extract the Mass Matrix
    op.integrator('GimmeMCK',1.0,0.0,0.0)
    op.analyze(1,0.0) 
    # Number of equations in the model
    N = op.systemSize()         # Has to be done after analyze
    Mmatrix = op.printA('-ret') # Or use op.printA('-file','M.out')
    Mmatrix = np.array(Mmatrix) # Convert the list to an array
    Mmatrix.shape = (N,N)       # Make the array an NxN matrix
    print( '\n************************************************************', \
          '\nExtracting the mass matrix, ignore the warnings...')
        
    # Determine maximum number of DOFs/node used in the system
    NDF = 0
    for node in op.getNodeTags():
        temp = len(op.nodeDOFs(node))
        if temp > NDF: NDF = temp

    # Rerrange the mass matrix in accordance with nodelist obtained from getNodeTags() to obtain Muu
    DOFs = []       # List containing indices of unrestrained DOFs from the global mass matrix
    used = {}       # Dictionary with nodes and associated dof indices used in unrestrained mass matrix
    ldict = {}      # Dictionary containing influence vectors
    Mratios = {}     # Dictionary containing modal masses
    Mfactors = {}    # Dictionary containing modal participation factors
    for i in range(1,NDF+1):
        ldict[i] = np.zeros([N,1])
        Mratios[i] = np.zeros(numEigen)
        Mfactors[i] = np.zeros(numEigen)
        
    idx = 0 # new dof index to use in unrestrained mass matrix
    for node in op.getNodeTags():               # start iterating over each node
        used[node] = []                         # create the list of unrestrained dofs used for the current dof
        ndf = len(op.nodeDOFs(node))            # number of DOFs used for the the current node
        for j in range(ndf):                    # iterate over used DOFs for the current node
            temp = op.nodeDOFs(node)[j]         # get idx of this dof, if -1 (restrained)
            if temp not in DOFs and temp >= 0:  # if it is unrestrained and not in DOFs list
                DOFs.append(temp)               # add to the unrestrained DOFs list
                used[node].append(j+1)          # for current node add to the used dof list (new idx in Muu)
                ldict[j+1][idx,0] = 1           # Assign 1 in influence vector if dof is in dir j
                idx += 1

    # Get the unrestrained part of mass matrix, Muu
    Mmatrix = Mmatrix[DOFs,:][:,DOFs]           

    # Calculate the total masses in unrestrained dofs
    Mtots = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
    for i in range(1,NDF+1):
        Mtots[i] = (ldict[i].T@Mmatrix@ldict[i])[0,0]

    op.wipeAnalysis()
    listSolvers = ['-genBandArpack','-fullGenLapack','-symmBandLapack']
    ok = 1  
    for s in listSolvers:
        print("Using %s as solver..." % s[1:])
        try:
            eigenValues = op.eigen(s,numEigen)
            catchOK = 0
            ok = 0
        except: 
            catchOK = 1
        
        if catchOK==0:
            for i in range(numEigen):
                if eigenValues[i] < 0: 
                    ok = 1
            if ok==0: 
                print('Eigenvalue analysis is completed.')
                break
    if ok!=0:
        print("Error on eigenvalue something is wrong...")
        sys.exit()
    else:  
        Lambda = np.asarray(eigenValues)
        Omega = Lambda**0.5
        T = 2*np.pi/Omega
        frq = 1/T

    for mode in range(1,numEigen+1):
        idx = 0
        phi = np.zeros([N,1])                         # Eigen vectors
        for node in used:
            for dof in used[node]:
                phi[idx,0]=op.nodeEigenvector(node,mode,dof)
                idx += 1
                
        phi = phi/(phi.T@Mmatrix@phi)**0.5  # Normalize the eigen vector by modal mass
        Mn = phi.T@Mmatrix@phi              # Modal mass (should always be equal to 1)

        for j in range(1,NDF+1):
            if Mtots[j] != 0:                              # Check if any mass is assigned
                Ln = phi.T@Mmatrix@ldict[j]                # Modal excitation factor
                Mnstar = (Ln**2/Mn)[0,0]                   # Effective modal mass
                Mfactors[j][mode-1] = Ln/Mn                # Modal participation factor
                Mratios[j][mode-1] = (Mnstar/Mtots[j]*100) # Effective modal mass participation ratio [%]
    
    for j in range(1,7):
        try: 
            Mratios[j]
        except: 
            Mratios[j]  = np.zeros(numEigen)
            Mfactors[j] = np.zeros(numEigen)

    # Calculate cumulative modal mass participation ratio
    sM1 = np.cumsum(Mratios[1]); sM2 = np.cumsum(Mratios[2]); sM3 = np.cumsum(Mratios[3])
    sM4 = np.cumsum(Mratios[4]); sM5 = np.cumsum(Mratios[5]); sM6 = np.cumsum(Mratios[6])      

    # Print results to the .csv file
    if outname != None:
        with open(outname+'.csv','w', encoding='utf-32') as f:
            f.write('Modal Periods and Frequencies\n')
            f.write('Mode,T [sec],f [Hz],\u03C9 [rad/sec],\u03BB [rad\u00b2/sec\u00b2]\n')
            for mode in range(numEigen):      
                f.write('%s,%s,%s,%s,%s\n' \
                      % ("{:.0f}".format(mode+1), "{:.3f}".format(T[mode]), "{:.3f}".format(frq[mode]), \
                         "{:.2f}".format(Omega[mode]), "{:.2f}".format(Lambda[mode])))

            f.write('\nTotal Mass of the Structure\n')
            f.write('M\u2081 [%],M\u2082 [%],M\u2083 [%],M\u2084 [%],M\u2085 [%],M\u2086 [%]\n')
            f.write('%s,%s,%s,%s,%s,%s\n' \
                    % ( "{:.2f}".format(Mtots[1]), "{:.2f}".format(Mtots[2]), "{:.2f}".format(Mtots[3]), \
                        "{:.2f}".format(Mtots[4]), "{:.2f}".format(Mtots[5]), "{:.2f}".format(Mtots[6])))

            f.write('\nModal Mass Participation Factors\n') 
            f.write('Mode,\u0393\u2081,\u0393\u2082,\u0393\u2083,\u0393\u2084,\u0393\u2085,\u0393\u2086\n')               
            for mode in range(numEigen):
                f.write('%s,%s,%s,%s,%s,%s,%s\n' % ("{:.0f}".format(mode+1), \
                    "{:.3f}".format(Mfactors[1][mode]), "{:.3f}".format(Mfactors[2][mode]), "{:.3f}".format(Mfactors[3][mode]), \
                    "{:.3f}".format(Mfactors[4][mode]), "{:.3f}".format(Mfactors[5][mode]), "{:.3f}".format(Mfactors[6][mode])))  

            f.write('\nEffective Modal Mass Participation Ratios\n') 
            f.write('Mode,U\u2081 [%],U\u2082 [%],U\u2083 [%],U\u2084 [%],U\u2085 [%],U\u2086 [%]\n')               
            for mode in range(numEigen):
                f.write('%s,%s,%s,%s,%s,%s,%s\n' % ("{:.0f}".format(mode+1), \
                    "{:.2f}".format(Mratios[1][mode]), "{:.2f}".format(Mratios[2][mode]), "{:.2f}".format(Mratios[3][mode]), \
                    "{:.2f}".format(Mratios[4][mode]), "{:.2f}".format(Mratios[5][mode]), "{:.2f}".format(Mratios[6][mode])))  

            f.write('\nCumulative Effective Modal Mass Participation Ratios\n') 
            f.write('Mode,\u2211U\u2081 [%],\u2211U\u2082 [%],\u2211U\u2083 [%],\u2211U\u2084 [%],\u2211U\u2085 [%],\u2211U\u2086 [%]\n')               
            for mode in range(numEigen):
                f.write('%s,%s,%s,%s,%s,%s,%s\n' % ("{:.0f}".format(mode+1), \
                    "{:.2f}".format(sM1[mode]), "{:.2f}".format(sM2[mode]), "{:.2f}".format(sM3[mode]), \
                    "{:.2f}".format(sM4[mode]), "{:.2f}".format(sM5[mode]), "{:.2f}".format(sM6[mode])))  

    # Print modal analysis results to the screen
    if pflag == 1:
        print('\nModal Periods and Frequencies')
        print('%4s|%8s|%10s|%12s|%12s' \
              % ('Mode', 'T [sec]','f [Hz]','\u03C9 [rad/sec]', '\u03BB [rad\u00b2/sec\u00b2]'))
        for mode in range(numEigen):      
            print('%4s|%8s|%10s|%12s|%12s' \
                  % ("{:.0f}".format(mode+1), "{:.3f}".format(T[mode]), "{:.3f}".format(frq[mode]), \
                     "{:.2f}".format(Omega[mode]), "{:.2f}".format(Lambda[mode])))

        print('\nTotal Mass of the Structure')
        print('%8s|%8s|%8s|%8s|%8s|%8s' \
              % ('M\u2081','M\u2082','M\u2083','M\u2084','M\u2085','M\u2086'))
        print('%8s|%8s|%8s|%8s|%8s|%8s' \
                % ( "{:.1f}".format(Mtots[1]), "{:.1f}".format(Mtots[2]), "{:.1f}".format(Mtots[3]), \
                    "{:.1f}".format(Mtots[4]), "{:.1f}".format(Mtots[5]), "{:.1f}".format(Mtots[6])))

        print('\nModal Mass Participation Factors') 
        print('%4s|%7s|%7s|%7s|%7s|%7s|%7s' \
            % ('Mode','\u0393\u2081','\u0393\u2082','\u0393\u2083','\u0393\u2084','\u0393\u2085','\u0393\u2086') )             
        for mode in range(numEigen):
            print('%4s|%7s|%7s|%7s|%7s|%7s|%7s' % ("{:.0f}".format(mode+1), \
                "{:.3f}".format(Mfactors[1][mode]), "{:.3f}".format(Mfactors[2][mode]), "{:.3f}".format(Mfactors[3][mode]), \
                "{:.3f}".format(Mfactors[4][mode]), "{:.3f}".format(Mfactors[5][mode]), "{:.3f}".format(Mfactors[6][mode])))  

        print('\nEffective Modal Mass Participation Ratios [%]') 
        print('%4s|%6s|%6s|%6s|%6s|%6s|%6s' \
            % ('Mode','U\u2081','U\u2082','M\u2083','U\u2084','U\u2085','U\u2086') )              
        for mode in range(numEigen):
            print('%4s|%6s|%6s|%6s|%6s|%6s|%6s' % ("{:.0f}".format(mode+1), \
                "{:.2f}".format(Mratios[1][mode]), "{:.2f}".format(Mratios[2][mode]), "{:.2f}".format(Mratios[3][mode]), \
                "{:.2f}".format(Mratios[4][mode]), "{:.2f}".format(Mratios[5][mode]), "{:.2f}".format(Mratios[6][mode])))  

        print('\nCumulative Effective Modal Mass Participation Ratios [%]') 
        print('%4s|%6s|%6s|%6s|%6s|%6s|%6s' \
            % ('Mode','\u2211U\u2081','\u2211U\u2082','\u2211U\u2083','\u2211U\u2084','\u2211U\u2085','\u2211U\u2086') )              
        for mode in range(numEigen):
            print('%4s|%6s|%6s|%6s|%6s|%6s|%6s' % ("{:.0f}".format(mode+1), \
                "{:.2f}".format(sM1[mode]), "{:.2f}".format(sM2[mode]), "{:.2f}".format(sM3[mode]), \
                "{:.2f}".format(sM4[mode]), "{:.2f}".format(sM5[mode]), "{:.2f}".format(sM6[mode])))  

    return T, Mratios, Mfactors, Mtots
