import scipy as sp
import numpy as np

def LU(A):
    '''
    LU(A)
    
    Function to perform LU Decomposition on a matrix (A) using the Crout algorithm.
    L[i][i] set to 1 (Doolittle).
    Outputs Lower (L), Upper (U) and result (Res) matrices.
    '''
    n = len(A) #the value of n for the nxn matrix A
    U = sp.zeros([n,n]) #creates the nxn matrix to hold the upper matrix
    L = sp.zeros([n,n]) #creates the nxn matrix to hold the lower matrix
    for j in range (n):
        L[j][j] = 1.0 #sets the diagonal in the lower matrix to have values of 1
        for i in range(j+1): #itterates through upper matrix in order
            U[i][j] = A[i][j] - sum(L[i][k]*U[k][j] for k in range(n)) #equation for calculating the upper matrix values (Crout's algorithm)
        for i in range(j+1,n): #itterates through lower matrix in order
            L[i][j] = (A[i][j] - sum(L[i][k]*U[k][j] for k in range(n)))/U[j][j] #equation for calculating the lower matrix values (Crout's algorithm)
    LUCheck(A,L,U) #runs the check function to make sure everything is working properly
    Res = U
    for i in range(n):
        for j in range(n):
            if j<i:
                Res[i][j] = L[i][j]
    return (L,U,Res)

def dt(A):
    '''
    dt(A)
    
    Function to find the determinant of a matrix (A) from LU decomposition.
    Outputs determinant (d).
    '''
    n = len(A) #the value of n for the nxn matrix A
    L,U, Res = LU(A) #runs LU decomposition
    d = 1 #creates variable for out determinant
    for i in range(n):
        d = d*U[i][i] #calculates product of diagonal elements of upper matrix
    dtCheck(A,d) #runs check function to make sure everything is working correctly
    return d

def solve(A,b):
    '''
    solve(A,b)
    
    Function to solve A.x=b for x given a matrix (A) and solution (b).
    Uses backwards and forwards substitution with LU decomposition.
    Outputs x.
    '''
    n = len(A) #the value of n for the nxn matrix A
    L,U, Res = LU(A) #runs LU decomposition
    #with L.U.x=b we can say L.y=b and U.x=y
    y = sp.zeros(n) #column vector to hold our x values
    x = sp.zeros(n) #column vector to hold our y values
    y[0]=b[0]/L[0][0] #equation for y0
    for i in range(1,n):
        y[i] = (b[i]-sum(L[i][j]*y[j] for j in range(i)))/L[i][i] #equation for solving the elements of y
    x[n-1]=y[n-1]/(U[n-1][n-1]) #equation for x0
    for i in range(n-2,-1,-1):
        x[i] = (y[i]-sum(U[i][j]*x[j] for j in range(i+1,n)))/U[i][i] #equation for solving the elements of x
    solveCheck(A,x,b) #runs check function to make sure everything is working correctly
    return x
    
def inv(A):
    '''
    inv(A)
    
    Function to find the inverse of a matrix (A) using LU decompostion.
    Outputs the inverse (inv)
    '''
    n = len(A) #the value of n for the nxn matrix A
    L,U, Res = LU(A) #runs LU decomposition
    I = sp.zeros([n,n]) #creates nxn matrix that will become the identity matrix
    inv = sp.zeros([n,n]) #creates nxn matrix that will become the inverse matrix
    for i in range(n):
        I[i][i] = 1.0 #makes I the identity matrix
    for j in range(n):
        Ii = sp.zeros(n) #column vector that will be column j of the identity matrix
        for i in range(n):
            Ii[i] = I[i][j] #sets Ii to be column j of the identity matrix
        #now for each column we can say A.inv(A)[j]=I[j]
        invi = solve(A,Ii) #solve the equation with our function
        for i in range(n):
            inv[i][j] = invi[i] #assembles the columns into our matrix
    invCheck(A,inv) #runs check function to make sure everything is working correctly
    return inv
        
def LUCheck(A,L,U):
    '''
    LUCheck(A,L,U)
    
    Checks the LU function is working correctly by confirming L.U=A.
    If not an Error is raised.
    '''
    n = len(A) #the value of n for the nxn matrix A
    A1 = sp.matmul(L,U) #multiplies L and U (should give A)
    for i in range (n):
        for j in range (n):
            if abs(A[i][j]-A1[i][j]) > 1e-14: #ensures the product of L and U is ~ A
                raise TypeError('LU Function Failed') #if not an error is raised
    return

def dtCheck(A,d):
    '''
    dtCheck(A,d)
    
    Checks the dt function is working correctly by comparing determinant 
    calculated determinant (d) to numpy function.
    If not an Error is raised.
    '''
    det = np.linalg.det(A) #calculates determinant with nump function
    if abs(det-d) > 1e-8: #ensures our determinant value matches the true value
        raise TypeError('Det Function Failed') #if not an error is raised
    return

def solveCheck(A,x,b):
    '''
    solveCheck(A,x,b)
    
    Checks the solve function is working correctly by confirming A.x=b and by
    comparing calculated x to numpy function
    If not an Error is raised.
    '''
    n = len(A) #the value of n for the nxn matrix A
    check = np.dot(A,x) - b  #calculates A.x - b (should be 0)
    true = np.linalg.solve(A, b) #calculates x with numpy function
    for i in range(n):
        if abs(check[i]) > 1e-14: #checks A.x - b = 0
            raise TypeError('Solve Function Failed') #if not an error is raised
        if abs(true[i]-x[i]) > 1e-14: #checks my x matches the numpy x
            raise TypeError('Solve Function Failed') #if not an error is raised
    return

def invCheck(A,inv):
    '''
    invCheck(A,inv)
    
    Checks the inv function is working correctly by making sure inv(A).A=I
    If not an Error is raised.
    '''
    n = len(A) #the value of n for the nxn matrix A
    check = sp.matmul(A,inv) #calculates A.inv(A) (should be identity matrix)
    for i in range(n):
        if (check[i][i] - 1) > 1e-14: #checks diagonals are 1
            raise TypeError('Inv Function Failed') #if not an error is raised
        for j in range(n):
            if i!=j:
                if check[i][j] > 1e-14: #checks non-diagonals are 0
                    raise TypeError('Inv Function Failed')#if not an error is raised
    return

'''
This section of the code produces the results for the question.
The checks are activated throughout the functions themselves.
'''
A = sp.array([[3,1,0,0,0],
              [3,9,4,0,0],
              [0,8,20,10,0],
              [0,0,-22,31,-25],
              [0,0,0,-35,61]]) #our input matrix
b = sp.array([2,5,-4,8,9]) #the 'b' for our A.x=b equation
L,U, Res = LU(A) #run LU decompetition
x = solve(A, b) #solve for x
Det = dt(A) #find the determinant
Inv = inv(A) #find the inverse matrix
