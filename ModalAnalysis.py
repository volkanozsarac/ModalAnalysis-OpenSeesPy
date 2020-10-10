def ModalAnalysis(numEigen, outname=None):
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
        If values of mass participation ratio in a direction is -1
        the total unrestrained mass in this dof is equal to 0.
        
    Parameters
    ----------
    numEigen : int
        Number of eigenvalues to calculate.
    outname  : str, optional (The default is None)
        if not None, the modal properties for first numEigen modes
        will be printed into outname.csv.

    Returns
    -------
    T     : numpy.ndarray
        Period array for the first numEigen modes.
    Mdict : dictionary
        dictionary containing mass participation ratios for first numEigen
        modes in all 6 dofs.

    """
    import numpy as np
    import openseespy.opensees as op
    import sys

    # Determine maximum number of DOFs/node used in the system
    NDF = 0
    for node in op.getNodeTags():
        temp = len(op.nodeDOFs(node))
        if temp > NDF: NDF = temp

    print('Extracting the mass matrix, ignore the following warnings...', \
          '\n************************************************************')    
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

    # Rerrange the mass matrix in accordance with nodelist obtained from getNodeTags() to obtain Muu
    DOFs = []       # List containing indices of unrestrained DOFs from the global mass matrix
    used = {}       # Dictionary with nodes and associated dof indices used in unrestrained mass matrix
    ldict = {}      # Dictionary containing influence vectors
    Mdict = {}      # Dictionary containing modal masses
    for i in range(1,NDF+1):
        ldict['l'+str(i)] = np.zeros([N,1])     # influence vector (df i)
        Mdict['M'+str(i)] = np.zeros(numEigen)  # Modal masses (df i)
        
    idx = 0 # new dof index to use in unrestrained mass matrix
    for node in op.getNodeTags():               # start iterating over each node
        used[node] = []                         # create the list of unrestrained dofs used for the current dof
        ndf = len(op.nodeDOFs(node))            # number of DOFs used for the the current node
        for j in range(ndf):                    # iterate over used DOFs for the current node
            temp = op.nodeDOFs(node)[j]         # get idx of this dof, if -1 (restrained)
            if temp not in DOFs and temp >= 0:  # if it is unrestrained and not in DOFs list
                DOFs.append(temp)               # add to the unrestrained DOFs list
                used[node].append(j+1)          # for current node add to the used dof list (new idx in Muu)
                ldict['l'+str(j+1)][idx,0] = 1  # Assign 1 in influence vector if dof is in dir j
                idx += 1

    Mmatrix = Mmatrix[DOFs,:][:,DOFs]           # Unrestrained mass matrix or Muu

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
        Lamda = np.asarray(eigenValues)
        Omega = Lamda**0.5
        T = 2*np.pi/Omega
        frq = 1/T

    for mode in range(1,numEigen+1):
        idx = 0
        phi = np.zeros([N,1])                          # Eigen vector for the specified mode
        for node in used:
            for dof in used[node]:
                phi[idx,0]=op.nodeEigenvector(node,mode,dof)
                idx += 1
        
        for j in range(1,NDF+1):
            l = ldict['l'+str(j)]                      # Influence vector in specified dof
            Mtot = l.T@Mmatrix@l                       # Total mass in specified dof
            if not Mtot == 0:                          # Check if any mass is assigned in dof
                Mn = phi.T@Mmatrix@phi                 # Modal mass in specified direction
                Ln = phi.T@Mmatrix@l                   # Effective modal mass in specified dof
                Mnstar = (Ln**2/Mn/Mtot*100)[0,0]      # Normalised effective modal mass participating [%] in specified dof
            else: Mnstar = -1                          # No mass is assigned in specified dof, indicate it with -1.
            Mdict['M'+str(j)][mode-1] = Mnstar         # Save the modal mass participation ratio for the specified dof
    
    for j in range(1,7):
        try: 
            Mdict['M'+str(j)]
        except: 
            Mdict['M'+str(j)] = np.ones(numEigen)*(-1)
            
    if outname != None:
        M1 = Mdict['M1']; M2 = Mdict['M2']; M3 = Mdict['M3']; M4 = Mdict['M4']; M5 = Mdict['M5']; M6 = Mdict['M6']
        sM1 = np.cumsum(M1); sM2 = np.cumsum(M2); sM3 = np.cumsum(M3); sM4 = np.cumsum(M4); sM5 = np.cumsum(M5); sM6 = np.cumsum(M6)
        if not any(sM1>0): sM1 = np.ones(numEigen)*(-1)
        if not any(sM2>0): sM2 = np.ones(numEigen)*(-1)
        if not any(sM3>0): sM3 = np.ones(numEigen)*(-1)
        if not any(sM4>0): sM4 = np.ones(numEigen)*(-1)
        if not any(sM5>0): sM5 = np.ones(numEigen)*(-1)
        if not any(sM6>0): sM6 = np.ones(numEigen)*(-1)
        
        with open(outname+'.csv','w') as f:
            f.write('Mode,T [sec],f [Hz],\u03C9 [rad/sec],M1 [%],M2 [%],M3 [%],M4 [%],M5 [%],M6 [%],\u2211M1 [%],\u2211M2 [%],\u2211M3 [%],\u2211M4 [%],\u2211M5 [%],\u2211M6 [%]\n')
            for mode in range(1,numEigen+1):      
                f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
                      % ("{:.0f}".format(mode), "{:.3f}".format(T[mode-1]), "{:.3f}".format(frq[mode-1]), "{:.2f}".format(Omega[mode-1]), \
                          "{:.2f}".format(M1[mode-1]), "{:.2f}".format(M2[mode-1]), "{:.2f}".format(M3[mode-1]), \
                          "{:.2f}".format(M4[mode-1]), "{:.2f}".format(M5[mode-1]), "{:.2f}".format(M6[mode-1]), \
                          "{:.2f}".format(sM1[mode-1]), "{:.2f}".format(sM2[mode-1]), "{:.2f}".format(sM3[mode-1]), \
                          "{:.2f}".format(sM4[mode-1]), "{:.2f}".format(sM5[mode-1]), "{:.2f}".format(sM6[mode-1])))

    return T, Mdict
