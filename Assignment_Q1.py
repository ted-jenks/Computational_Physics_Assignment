import scipy as sp
import matplotlib.pyplot as plt

def near(r):
    '''
    near(r)
    
    Function to find the nearest representable real numbers higher and lower than a 
    given floating point value (r).
    Also finds the upper and lower difference within which real numbers are rounded to
    the given floating point value and expresses it as a fractional range
    Outputs the upper number, lower number and the fractional rounding range.
    '''
    x = r + r/4 #selects a value slightly higher than input number
    upper = float('NaN') #create variable 'upper' to hold nearest representable number higher than r
    #perform a binary search style algorithm until we 'run out' of numbers
    while x != r and x != upper: #condition only met when there are no floats between our values
        upper = x #stores the previous value of x
        x = (x+r)/2 #sets the next value of x to the midpoint
    y = r - r/4 #selects a value slightly lower than input number
    lower = float('NaN') #create variable 'lower' to hold nearest representable number lower than r
    #perform a binary search style algorithm until we 'run out' of numbers
    while y != r and y != lower: #condition only met when there are no floats between our values
        lower = y #stores the previous value of y
        y = (y+r)/2 #sets the next value of x to the midpoint
    difU = (upper - r)/2 #finds the range above r that rounds to r
    difL = (r-lower)/2 #finds the range below r that rounds to r
    fracRange = (difU+difL)/r #finds the range that rounds to r as a fraction of r
    return (upper, lower, fracRange)

def nearCheck(r):
    '''
    nearCheck(r)
    
    Checks the near function by ensuring the upper and lower of a given floating
    point number in the near function return the given floating point number as one
    of their upper/lower accordingly.
    If not an Error is raised.
    '''
    upper, lower, fracRange = near(r) #performs our 'near' function on a value r
    upper1, lower1, fracRange1 = near(upper) #performs our 'near' function on a the 'upper' value for r
    upper2, lower2, fracRange2 = near(lower) #performs our 'near' function on a the 'lower' value for r
    if lower1 != r or upper2 != r: #confirms the lower of upper is r (and visa versa)
        raise TypeError('Near Function Failed') #if not an error is raised
    return

'''
This section of code returns the results for the question and runs a check
to validate the program.
'''
nearCheck(0.25) #this performs the check that my function is working properly
upper, lower, fracRange = near(0.25) #this Performs the function
print('Upper from 0.25:',upper) #prints the upper from 0.25
print('Lower from 0.25:',lower) #prints the lower from 0.25
print('Fractional Range Rounding to 0.25:\n', fracRange) #prints the fractional range rounding to 0.25
upper1, lower1, fracRange1 = near(upper) #performs the function on the number higher than 0.25
print('Upper from the upper from 0.25:',upper1) #prints the upper the upper from 0.25
print('Lower from the lower from 0.25:',lower1) #prints the lower the upper from 0.25
print('Fractional Range Rounding to next representable real number higher than 0.25:\n',fracRange1) #prints the fractional range rounding to the number higher than 0.25
upper2, lower2, fracRange2 = near(lower) #performs the function on the number lower than 0.25
print('Upper from the lower from 0.25:',upper2) #prints the upper the lower from 0.25
print('Lower from the lower from 0.25:',lower2) #prints the lower the lower from 0.25
print('Fractional Range Rounding to next representable real number lower than 0.25:\n',fracRange2) #prints the fractional range rounding to the number lower than 0.25

'''
This section of code produces the graph.
The aim of this is to show how fractional error range changes with increasing numbers.
'''
x = sp.linspace(0.06251,1.999999,10000) #create a big set of numbers to put through the function
x2 = [2,1,2**-1,2**-2,2**-3,2**-4] #create some powers of 2s to go through the function
FracRanges = [] #create a list to store the fractional ranges
FracRanges2 = [] #create a list to store the fractional ranges of the powers of 2s
for i in range(len(x)):
    upper, lower, fracRange = near(x[i]) #runs the function on big set of numbers
    FracRanges.append(fracRange) #stores the fractional ranges
for i in range(len(x2)):
    upper, lower, fracRange = near(x2[i])#runs the function on the powers of two
    FracRanges2.append(fracRange) #stores the fractional ranges
plt.figure()
plt.plot(x,FracRanges,'.') #plots all the fractional ranges
plt.plot(x2,FracRanges2,'x', label = 'Powers of Two') #plots all the fractional ranges for the powers of two
plt.grid()
plt.legend(frameon=False)
plt.xlabel('Number')
plt.ylabel('Fractional Error')
plt.title('Fractional Error for Rounding to Representable Numbers')
plt.show()
