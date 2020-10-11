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
        unrestrained dofs joints.
        
    Parameters
    ----------
    numEigen : int
        Number of eigenvalues to calculate.
    outname  : str, optional (The default is None)
        if not None, the modal properties for first numEigen modes
        will be printed into outname.csv.
    pflag    : int (1 or 0)
        flag to print output information on screen

    Returns
    -------
    T      : numpy.ndarray
        Period array for the first numEigen modes.
    Mratio : dictionary
        Mass participation ratios for first numEigen modes.
    Mtots  : dictionary
        Total mass of the structure

    """
    import numpy as np
    import openseespy.opensees as op
    import sys
    import time

    # Determine maximum number of DOFs/node used in the system
    NDF = 0
    for node in op.getNodeTags():
        temp = len(op.nodeDOFs(node))
        if temp > NDF: NDF = temp

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
    time.sleep(0.5)
    print( '\n************************************************************', \
          '\nExtracting the mass matrix, ignore the previous warnings...')   
    # Rerrange the mass matrix in accordance with nodelist obtained from getNodeTags() to obtain Muu
    DOFs = []       # List containing indices of unrestrained DOFs from the global mass matrix
    used = {}       # Dictionary with nodes and associated dof indices used in unrestrained mass matrix
    ldict = {}      # Dictionary containing influence vectors
    Mratio = {}     # Dictionary containing modal masses
    for i in range(1,NDF+1):
        ldict[i] = np.zeros([N,1])      # Dictionary containing influence vectors
        Mratio[i] = np.zeros(numEigen)  # Dictionary containing modal masses
        
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
                if eigenValues[i] < 0: ok = 1
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
        
        for j in range(1,NDF+1):
            if Mtots[j] != 0:                          # Check if any mass is assigned
                Mn = phi.T@Mmatrix@phi                 # Modal mass
                Ln = phi.T@Mmatrix@ldict[j]            # Effective modal mass
                Mnstar = (Ln**2/Mn/Mtots[j]*100)[0,0]  # Normalised effective modal mass participating [%]
            else: Mnstar = 0                           # No mass is assigned, insert 0, otherwise results is nan
            Mratio[j][mode-1] = Mnstar                 # Save the modal mass participation ratio
    
    for j in range(1,7):
        try: 
            Mratio[j]
        except: 
            Mratio[j] = np.zeros(numEigen)

    # Calculate cumulative modal mass participation ratio
    sM1 = np.cumsum(Mratio[1]); sM2 = np.cumsum(Mratio[2]); sM3 = np.cumsum(Mratio[3])
    sM4 = np.cumsum(Mratio[4]); sM5 = np.cumsum(Mratio[5]); sM6 = np.cumsum(Mratio[6])      

    # Print results to the .csv file
    if outname != None:
        with open(outname+'.csv','w') as f:
            f.write('Eigenvalue Analysis\n')
            f.write('Mode,T [sec],f [Hz],\u03C9 [rad/sec],\u03BB [rad\u00b2/sec\u00b2]\n')
            for mode in range(numEigen):      
                f.write('%s,%s,%s,%s,%s\n' \
                      % ("{:.0f}".format(mode+1), "{:.3f}".format(T[mode]), "{:.3f}".format(frq[mode]), \
                         "{:.2f}".format(Omega[mode]), "{:.2f}".format(Lambda[mode])))

            f.write('\nTotal Mass of the Structure\n')
            f.write('M1 [%],M2 [%],M3 [%],M4 [%],M5 [%],M6 [%]\n')
            f.write('%s,%s,%s,%s,%s,%s\n' \
                    % ( "{:.2f}".format(Mtots[1]), "{:.2f}".format(Mtots[2]), "{:.2f}".format(Mtots[3]), \
                        "{:.2f}".format(Mtots[4]), "{:.2f}".format(Mtots[5]), "{:.2f}".format(Mtots[6])))

            f.write('\nModal Mass Participation Ratio\n') 
            f.write('Mode,M1 [%],M2 [%],M3 [%],M4 [%],M5 [%],M6 [%]\n')               
            for mode in range(numEigen):
                f.write('%s,%s,%s,%s,%s,%s,%s\n' % ("{:.0f}".format(mode+1), \
                    "{:.2f}".format(Mratio[1][mode]), "{:.2f}".format(Mratio[2][mode]), "{:.2f}".format(Mratio[3][mode]), \
                    "{:.2f}".format(Mratio[4][mode]), "{:.2f}".format(Mratio[5][mode]), "{:.2f}".format(Mratio[6][mode])))  

            f.write('\nCumulative Modal Mass Participation Ratio\n') 
            f.write('Mode,\u2211M1 [%],\u2211M2 [%],\u2211M3 [%],\u2211M4 [%],\u2211M5 [%],\u2211M6 [%]\n')               
            for mode in range(numEigen):
                f.write('%s,%s,%s,%s,%s,%s,%s\n' % ("{:.0f}".format(mode+1), \
                    "{:.2f}".format(sM1[mode]), "{:.2f}".format(sM2[mode]), "{:.2f}".format(sM3[mode]), \
                    "{:.2f}".format(sM4[mode]), "{:.2f}".format(sM5[mode]), "{:.2f}".format(sM6[mode])))  

    # Print modal analysis results to the screen
    if pflag == 1:
        print('\nEigenvalue Analysis Results')
        print('Mode| T [sec] | f [Hz] | \u03C9 [rad/sec] | \u03BB [rad\u00b2/sec\u00b2]')
        for mode in range(numEigen):      
            print('%3s |%7s  |%6s  |%9s    |%9s' \
                  % ("{:.0f}".format(mode+1), "{:.3f}".format(T[mode]), "{:.3f}".format(frq[mode]), \
                     "{:.2f}".format(Omega[mode]), "{:.2f}".format(Lambda[mode])))

        print('\nTotal Mass of the Structure')
        print('%8s|%8s|%8s|%8s|%8s|%8s' \
              % ('M1','M2','M3','M4','M5','M6'))
        print('%8s|%8s|%8s|%8s|%8s|%8s' \
                % ( "{:.1f}".format(Mtots[1]), "{:.1f}".format(Mtots[2]), "{:.1f}".format(Mtots[3]), \
                    "{:.1f}".format(Mtots[4]), "{:.1f}".format(Mtots[5]), "{:.1f}".format(Mtots[6])))

        print('\nModal Mass Participation Ratio') 
        print('%4s|%8s|%8s|%8s|%8s|%8s|%8s' \
            % ('Mode','M1 [%]','M2 [%]','M3 [%]','M4 [%]','M5 [%]','M6 [%]') )              
        for mode in range(numEigen):
            print('%4s|%8s|%8s|%8s|%8s|%8s|%8s' % ("{:.0f}".format(mode+1), \
                "{:.2f}".format(Mratio[1][mode]), "{:.2f}".format(Mratio[2][mode]), "{:.2f}".format(Mratio[3][mode]), \
                "{:.2f}".format(Mratio[4][mode]), "{:.2f}".format(Mratio[5][mode]), "{:.2f}".format(Mratio[6][mode])))  

        print('\nModal Mass Participation Ratio') 
        print('%4s|%8s|%8s|%8s|%8s|%8s|%8s' \
            % ('Mode','\u2211M1 [%]','\u2211M2 [%]','\u2211M3 [%]','\u2211M4 [%]','\u2211M5 [%]','\u2211M6 [%]') )              
        for mode in range(numEigen):
            print('%4s|%8s|%8s|%8s|%8s|%8s|%8s' % ("{:.0f}".format(mode+1), \
                "{:.2f}".format(sM1[mode]), "{:.2f}".format(sM2[mode]), "{:.2f}".format(sM3[mode]), \
                "{:.2f}".format(sM4[mode]), "{:.2f}".format(sM5[mode]), "{:.2f}".format(sM6[mode])))  

    return T, Mratio, Mtots
