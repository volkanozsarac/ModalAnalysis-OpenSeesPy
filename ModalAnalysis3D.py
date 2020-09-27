def ModalAnalysis3D(numEigen):
    """
    |-----------------------------------------------------------------------|
    |                                                                       |
    |    Modal Analysis of 3D systems                                       |
    |                                                                       |
    |    Author: Volkan Ozsarac                                             |
    |    Affiliation: University School for Advanced Studies IUSS Pavia     |
    |    Earthquake Engineering PhD Candidate                               |
    |                                                                       |
    |-----------------------------------------------------------------------|
    """
    import numpy as np
    import openseespy.opensees as op
    import time
    import sys

    print('Extracting the mass matrix, ignore the following warnings...\n')
    op.wipeAnalysis()
    op.system('FullGeneral')
    op.analysis('Transient')

    # Extract the Mass Matrix
    op.integrator('GimmeMCK',1.0,0.0,0.0)
    op.analyze(1,0.0) 
    time.sleep(0.5) 
    # Number of equations in the model
    N = op.systemSize() # Has to be done after analyze
    Mmatrix = op.printA('-ret') # Or use op.printA('-file','M.out')
    Mmatrix = np.array(Mmatrix) # Convert the list to an array
    Mmatrix.shape = (N,N) # Make the array an NxN matrix
    print('\nExtracted the mass matrix, ignore the previous warnings...')

    # Rerrange the mass matrix in accordance with nodelist order from getNodeTags()
    DOFs = [] # These are the idx of all the DOFs used in the extract mass matrix, order is rearranged
    used = {} # Save here the nodes and their associated dofs used in global mass matrix
    lx = np.zeros([N,1]) # influence vector (x)
    ly = np.zeros([N,1]) # influence vector (y)
    lz = np.zeros([N,1]) # influence vector (z)
    lrx = np.zeros([N,1]) # influence vector (rx)
    lry = np.zeros([N,1]) # influence vector (ry)
    lrz = np.zeros([N,1]) # influence vector (rz)

    idx = 0 # index
    NDF = 6 # NDF is number of DOFs/node
    for node in op.getNodeTags():
        used[node] = []
        for j in range(NDF): 
            temp = op.nodeDOFs(node)[j]
            if temp not in DOFs and temp >=0:
                DOFs.append(op.nodeDOFs(node)[j])
                used[node].append(j+1)
                if j == 0: lx[idx,0] = 1
                if j == 1: ly[idx,0] = 1
                if j == 2: lz[idx,0] = 1
                if j == 3: lrx[idx,0] = 1
                if j == 4: lry[idx,0] = 1
                if j == 5: lrz[idx,0] = 1
                idx += 1

    Mmatrix = Mmatrix[DOFs,:][:,DOFs] 

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
        print("Error on Modal Analysis...")
        sys.exit()
    else:  
        Lamda = np.asarray(eigenValues)
        Omega = Lamda**0.5
        T = 2*np.pi/Omega
        f = 1/T

    print('Modal properties for the first %d modes:' % numEigen)
    Mx = []; My = []; Mz = []; Mrx = []; Mry = []; Mrz = []
    dofs = ['x','y','z','rx','ry','rz']
    print('Mode| T [sec] | f [Hz] | \u03C9 [rad/sec] | Mx [%] | My [%] | Mz [%] | \u2211Mx [%] | \u2211My [%] | \u2211Mz [%]')

    for mode in range(1,numEigen+1): # Although Mrx, Mry and Mrz are calculated, I am not printing these
        idx = 0
        phi = np.zeros([N,1])
        for node in used:
            for dof in used[node]:
                phi[idx,0]=op.nodeEigenvector(node,mode,dof)
                idx += 1
        
        for dof in dofs:
            l = eval('l'+dof)
            Mtot = l.T@Mmatrix@l   # Total mass in specified global dof
            Mn = phi.T@Mmatrix@phi # Modal mass in specified global dof
            Ln = phi.T@Mmatrix@l   # Effective modal mass in specified global dof
            Mnstar = Ln**2/Mn/Mtot*100 # Normalised effective modal mass participating [%], in specified global dof
            eval('M'+dof).append(Mnstar[0,0]) # Save the modal mass for the specified dof
        
        print('%3s |%7s  |%6s  |%9s    |%6s  |%6s  |%6s  |%7s  |%7s  |%7s' \
              % ("{:.0f}".format(mode), "{:.3f}".format(T[mode-1]), "{:.3f}".format(f[mode-1]), "{:.2f}".format(Omega[mode-1]), \
                  "{:.2f}".format(Mx[mode-1]), "{:.2f}".format(My[mode-1]), "{:.2f}".format(Mz[mode-1]), \
                  "{:.2f}".format(sum(Mx)), "{:.2f}".format(sum(My)), "{:.2f}".format(sum(Mz))))

    Mtot = lx.T@Mmatrix@lx; Mtot = Mtot[0,0]