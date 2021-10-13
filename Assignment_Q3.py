import scipy as sp
import Assignment_Q2 as Q2
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

def LagInt(t,x,y):
    '''
    LagInt(t,x,y)
    
    This function performs Lagrangian interpolation on tabulated data set x and
    y.
    Outputs the functional value (summed) for a given x.
    '''
    n = len(x) #calculates the number of data points
    summed = 0 #creates variable to hold function
    for i in range(0,n):
        prod = 1 #creates variable to hold product part of equation
        for j in range(0,n):
            if j != i:
                prod = prod * ((t-x[j])/(x[i]-x[j])) #calculate product part of equations
        summed = summed + prod*y[i] #calculates function
    return summed

def CubSpline(t,x,y):
    '''
    CubSpline(t,x,y)
    
    This function performs cubic Spline interpolation on tabulated data set x and
    y.
    Outputs the functional value (f) for a given x.
    '''
    n = len(x) #calculates the number of data points
    A = sp.zeros([n-2,n-2]) #creates matrix to hold values for the tridiagonal matrix that multiples the column vector of second derivatives
    b = sp.zeros(n-2) #creates column vector to hold solutions to the second differential matrix equation
    f = sp.zeros(len(t)) #creates variablle to hold the function
    diffs = sp.zeros(n) #creates column vector to hold second differentials
    j = 0 #sets up initial j value (for the position in A)
    i = 1 #sets up initial i value (for the equation of x)
    while j < n-2 and i < n-1:
        A[j][j] = ((x[i+1]-x[i-1])/3)
        if j!=0:
            A[j][j-1] = ((x[i]-x[i-1])/6)
        if j!=n-3:
            A[j][j+1] = ((x[i+1]-x[i])/6)
        #the three equations above populate the matrix A with the appropriate values
        b[j] = ((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1])) #calculates the values of b
        j = j + 1 #iterates over j
        i = i + 1 #iterates over i at the same time as j
    diffcalc = Q2.solve(A,b) #calculates the values of the second differentials (without f0'' and fn'')
    for i in range(1,n-1):
        diffs[i] = diffcalc[i-1] #sets f0'' and fn'' to 0
    for j in range(len(t)):    
        i = 0
        while t[j] > x[i]:
            i = i+1 #finds the interval a value of x is in
        i = i-1
        A = (x[i+1]-t[j])/(x[i+1]-x[i])
        B = 1-A
        C = ((A**3-A)*(x[i+1]-x[i])**2)/6
        D = ((B**3-B)*(x[i+1]-x[i])**2)/6
        #the three equations above are the coefficients for the cubic spline funtion
        f[j] = A*y[i]+B*y[i+1]+C*diffs[i]+D*diffs[i+1] #calculates the cubic spline function
    return f

def CubSplineTest():
    '''
    CubSplineTest()
    
    This function tests the cubic spline function by putting in a third-order polynomial.
    The third order polynomial should be the same as the cubic spline outputted,
    except for effects at boundary conditions.
    A graph is displayed and the error is checked to be below a certain threshold.
    '''
    t = sp.linspace(-6,6,5001) #create a set of values t
    x = sp.linspace(-10,10,15) #these are the x values we will sample on our function
    y = 1/3 * x**3  + 3*x**2 #the y values we will sample for our function
    yt = 1/3 * t**3  + 3*t**2 #the y values of the the original function
    c = CubSpline(t,x,y) #run the cubic spline function
    plt.figure()
    plt.title('Cubic Spline Test')
    plt.plot(t,c, label = 'Cubic Spline', color = 'blue') #plot the cubic spline
    plt.plot(t,yt,'--',dashes = (10,10), color = 'red', label = 'Third Order Polynomial') #plot over the expected function
    plt.legend(loc='upper left', frameon=False)
    plt.xlim(-6,6)
    plt.show()
    for i in range(len(t)):
        if abs(c[i]-yt[i]) > 0.1:
            raise TypeError('Cubic Spline Failed') #reais error if difference is above threshold
    return
    
    

x = sp.array([-0.75,-0.5,-0.35,-0.1,0.05,0.1,0.23,0.29,0.48,0.6,0.92,1.05,1.5]) #our input x
y = sp.array([0.1,0.3,0.47,0.66,0.6,0.54,0.3,0.15,-0.32,-0.54,-0.6,-0.47,-0.08]) #our input y
t = sp.linspace(-.75,1.5,1000) #the ts to plot against

fig, axs = plt.subplots(2) #create subplots
fig.tight_layout(pad=3.0) 
fig.suptitle('Interpolation Plots')

axs[0].set_ylim(-120,10)
axs[0].set_xlabel('x')
axs[0].set_ylabel('y')
axs[0].plot(x,y,'x')
axs[0].plot(t,LagInt(t,x,y), label = 'Lagrangian Polynomial', color = 'blue') #plot the langrangian
axs[0].plot(x,y,'x')
axs[0].plot(t,CubSpline(t,x,y),label = 'Cubic Spline', color = 'red') #plot the cubic spline
axs[0].legend(loc='lower left', frameon=False)

axs[1].set_ylim(-1,1)
axs[1].set_xlabel('x')
axs[1].set_ylabel('y')
c = CubicSpline(x,y,bc_type='natural') #the built in cubic spline function with natural BCs
axs[1].plot(x,y,'x')
axs[1].plot(t,LagInt(t,x,y), label = 'Lagrangian Polynomial', color = 'blue') #plot the langrangian
axs[1].plot(x,y,'x')
axs[1].plot(t,CubSpline(t,x,y),label = 'Cubic Spline', color = 'red') #plot the cubic spline
axs[1].plot(t,c(t),'--',dashes=(5, 15),label = 'Scipy Cubic Spline', color = 'black')
axs[1].legend(loc='upper right', frameon=False)

CubSplineTest() #run the test function