import scipy as sp
import matplotlib.pyplot as plt

params = {
   'axes.labelsize': 12,
   'font.size': 10,
   'legend.fontsize': 15,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'figure.figsize': [15, 15]
   } 
plt.rcParams.update(params)

def ImpDiv(t, Vin, Vinh, method):
    '''
    ImpDiv(t, Vin, Vinh, method)
    
    This function calculate the voltage output of an impedance divider for a given
    input voltage function (Vin).
    It can use two different ODE methods, fourth-order Adams-Bashforth (AB4) and
    fourth-orde Runge-Kutta (RB4) which need to be chosen between when running the function.
    The inupt Vinh is the function of Vin(t+h/2) where h is time spacing to allow the 
    Runge-Kutta method to work.
    '''
    n = len(Vin) #finds the number of samples of the input function
    Vout = sp.zeros(n) #create variable to hold the output voltage
    h = t[1]-t[0] #finds the time spacing
    Vout[0] = 1 #for the initial value Vout=Vin as time isn't changing in this period
    for i in range(0,3):
        fa = Vin[i] - Vout[i]
        fb = Vinh[i] - (Vout[i] + h*fa/2)
        fc = Vinh[i] - (Vout[i] + h*fb/2)
        fd = Vin[i+1] - (Vout[i] + h*fc)
        Vout[i+1] = (Vout[i] + (h/6) * (fa + 2*fb + 2*fc +fd)) #Use Runge-Kutta to get first few values for AB4
    if method == 'AB4': #for the fourth-order Adams-Bashforth method
        for i in range(3,n-1):
            fa = (Vin[i] - Vout[i])
            fb = (Vin[i-1]-Vout[i-1])
            fc = (Vin[i-2]-Vout[i-2])
            fd = (Vin[i-3]-Vout[i-3])
            Vout[i+1] = (Vout[i] + (h/24) * (55*fa-59*fb+37*fc-9*fd)) #these equations itteratively excecute the AB4 method
        return Vout
    if method == 'RK4':#for the fourth-order Runge-Kutta method
        for i in range(3,n-1):
            fa = Vin[i] - Vout[i]
            fb = Vinh[i] - (Vout[i] + h*fa/2)
            fc = Vinh[i] - (Vout[i] + h*fb/2)
            fd = Vin[i+1] - (Vout[i] + h*fc)
            Vout[i+1] = (Vout[i] + (h/6) * (fa + 2*fb + 2*fc +fd)) #these equations itteratively excecute the RB4 method
        return Vout
    else:
        raise TypeError('Invalid Method') #if the function is run with an invalid method name an error is raised
    return

def Vin1(t):
    '''
    Vin1(t)
    
    This function creates the Vin values for question 5(c) on the assignment.
    '''
    n = len(t) #finds the number of samples
    y = sp.zeros(n) #creates variable to hold signal
    for i in range(n):
        if t[i]<0:
            y[i] = 1 #sets V to 1 in negative t
        else:
            y[i] = 0 #sets V to 0 in positivee t
    return y

def Vin2(t,T):
    '''
    Vin1(t,T)
    
    This function creates the Vin values for question 5(e) on the assignment.
    It is a square function.
    It can be assigned a specific period, T.
    '''
    n = len(t) #finds the number of samples
    y = sp.zeros(n) #creates variable to hold signal
    for i in range(n):
        if t[i]<0:
            y[i] = 1 #in negative time, V is 1
        else:
            if int(2*(t[i]/T))%(2) == 0: #creates periodicity of square function
                y[i]=0
            else:
                y[i]=1
    return y
     
     
def ErrorV1(Vout, t):
    '''
    ErrorV1(Vout,t)
    
    This function finds the mean error of the method when it solves Vin1, the signal from
    5(c).
    It does this by comparing values to the analytically calculated output expected.
    It will output the mean error.
    '''
    errors = [] #create variable to hold the errors
    for i in range(1,len(t)):
        expect = sp.exp(-t) #the expected exponential decay function
        err = abs(Vout[i] - expect[i])/expect[i] #the difference between the expected and calculated values
        errors.append(err) #stores error values
    MeanErr = sp.average(errors) #calculates the mean of the errors
    return (MeanErr, errors)
  
       
t = sp.linspace(0,40,401) #the selected amount of sampling
h = t[1]-t[0] #find the timespacing

Vinf = Vin1(t) #the first input voltage
Vinfh = Vin1(t+h/2) #the first input voltage at points Vin(t+h/2) for the Runge-Kutta method
Vout1R = ImpDiv(t,Vinf,Vinfh, method='RK4') #find Vout with RK4
MeanErr, errors = ErrorV1(Vout1R,t)
print('Mean Relative Error from RB4 h = 0.1s:', MeanErr) #prints the mean error

#This code is identical and repeats the above for h = 0.05s
t = sp.linspace(0,40,801)
h = t[1]-t[0]

Vinf = Vin1(t)
Vinfh = Vin1(t+h/2)
Vout1R = ImpDiv(t,Vinf,Vinfh, method='RK4')

MeanErr, errors = ErrorV1(Vout1R,t)
print('Mean Relative Error from RB4 h = 0.05s:',MeanErr)

#This code is identical and repeats the above for h = 0.025s
fig, axs = plt.subplots(2,2)
fig.tight_layout(pad=9.0)
   
t = sp.linspace(0,40,1601)
tplot = sp.linspace(-10,0,801)
h = t[1]-t[0]

Vinf = Vin1(t)
Vinfh = Vin1(t+h/2)
Vout1A = ImpDiv(t,Vinf,Vinfh, method='AB4')
Vout1R = ImpDiv(t,Vinf,Vinfh, method='RK4')
te1 = sp.linspace(-4,0,1000)
te2 = sp.linspace(0,10,1000)
exp1 = te1**0
exp2 = sp.exp(-te2)

Vins = Vin2(t, 2)
Vinsh = Vin2(t+h/2, 2)
Vout2 = ImpDiv(t,Vins,Vinsh,method='RK4')

MeanErr, errors = ErrorV1(Vout1R,t)
print('Mean Relative Error from RB4 h = 0.025s:',MeanErr)

axs[0,0].plot(t,Vinf, color = 'cornflowerblue', label = 'Vin')
axs[0,0].plot(tplot,Vin1(tplot), color = 'cornflowerblue')
axs[0,0].plot(t,Vout1A , label = 'Vout by AB4',color = 'darkorange')
axs[0,0].plot(te1,exp1,'--',dashes=(5, 15), label = 'Expected',color='red')
axs[0,0].plot(te2,exp2,'--',dashes=(5, 15), color='red')
axs[0,0].legend(loc='upper right', frameon=False)
axs[0,0].set_xlim(-1,6)
axs[0,0].set_ylim(0,1.1)
axs[0,0].set_xlabel('t̂ (s)')
axs[0,0].set_ylabel('Voltage (V0)')
axs[0,0].set_title('Voltage Input 1 using Fourth-Order Adams-Bashforth Method')

axs[0,1].plot(t,Vinf, color = 'cornflowerblue', label = 'Vin')
axs[0,1].plot(tplot,Vin1(tplot), color = 'cornflowerblue')
axs[0,1].plot(t,Vout1R, label = 'Vout by RK4', color = 'green')
axs[0,1].plot(te1,exp1,'--',dashes=(5, 15), label = 'Expected',color='red')
axs[0,1].plot(te2,exp2,'--',dashes=(5, 15), color='red')
axs[0,1].legend(loc='upper right', frameon=False)
axs[0,1].set_xlim(-1,6)
axs[0,1].set_ylim(0,1.1)
axs[0,1].set_xlabel('t̂ (s)')
axs[0,1].set_ylabel('Voltage (V0)')
axs[0,1].set_title('Voltage Input 1 using Fourth-Order Runge-Kutta Method')

axs[1,0].plot(t,Vins, label = 'Vin', color = 'cornflowerblue')
axs[1,0].plot(tplot,Vin2(tplot,2), color = 'cornflowerblue')
axs[1,0].plot(t,Vout2 , label = 'Vout by RB4',color = 'green')
axs[1,0].legend(loc='upper right', frameon=False)
axs[1,0].set_xlim(-1,10)
axs[1,0].set_ylim(0,1.5)
axs[1,0].set_xlabel('t̂ (s)')
axs[1,0].set_ylabel('Voltage (V0)')
axs[1,0].set_title('Voltage Input 2 with T = 2RC using Fourth-Order Runge-Kutta Method')

Vins = Vin2(t, 0.5)
Vinsh = Vin2(t+h/2, 0.5)
Vout2 = ImpDiv(t,Vins,Vinsh,method='RK4')

axs[1,1].plot(t,Vins, label = 'Vin', color = 'cornflowerblue')
axs[1,1].plot(tplot,Vin2(tplot,0.5), color = 'cornflowerblue')
axs[1,1].plot(t,Vout2 , label = 'Vout by RB4',color = 'green')
axs[1,1].legend(loc='upper right', frameon=False)
axs[1,1].set_xlim(-1,10)
axs[1,1].set_ylim(0,1.5)
axs[1,1].set_xlabel('t̂ (s)')
axs[1,1].set_ylabel('Voltage (V0)')
axs[1,1].set_title('Voltage Input 2 with T = RC/2 using Fourth-Order Runge-Kutta Method')

plt.figure()
plt.plot(t[1:],errors)
plt.xlim(0,40)
plt.grid()
plt.ylim(0,1.4e-7)
plt.xlabel('t̂ (s)')
plt.ylabel('Relative Error')
plt.title('Relative Error with Time')

