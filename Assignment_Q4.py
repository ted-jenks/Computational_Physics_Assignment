import scipy as sp
import scipy.signal as spf
import numpy as np
import matplotlib.pyplot as plt


def tophat(t):
    '''
    tophat(t)
    
    This is the rect function returning the values (h) at points (t).
    '''
    n=len(t) #find the the number of data points
    h = sp.zeros(n) #variable to hold values of function
    for i in range(n):
        if t[i]<=7 and t[i]>=5:
            h[i] = 4 #function returns 4 in non-zero region
        else:
            h[i] = 0 #function returns 0 elsewhere
    return h

def gauss(t):
    '''
    gauss(t)
    
    This is the gaussian function returning the values (g) at points (t).
    '''
    n=len(t)  #find the the number of data points
    g = sp.zeros(n) #variable to hold values of function
    for i in range(n):
        g[i] = (1/(sp.sqrt(2*sp.pi)))*sp.exp((-1*t[i]**2)/4) #outputs values of the gaussian
    return g

def conv(t):
    '''
    conv(t)
    
    This is the function to convolve the gaussian and rect functions ussing FFTs.
    It outputs data for the fourier transfroms as well as the convolution.
    '''
    n = len(t) #find the the number of data points
    timestep = (t[1]-t[0]) #find the timestep between data points
    
    s1 = gauss(t) #our signal one is the gaussian
    np.pad(s1,(0,n-1),'constant') #pad signal one with n-1 zeroes to the right
    s2 = tophat(t) #our signal two is the rect function
    np.pad(s2,(0,n-1),'constant') #pad signal two with n-1 zeroes to the right
    Fg = np.fft.fft(s1, norm = 'ortho') #find the FFT os signal one
    Ft = np.fft.fft(s2, norm = 'ortho') #find the FFT os signal two
    
    xg = np.fft.fftshift(np.fft.fftfreq(n, d = timestep)) #finds the frequancies for the fourier transform plot
    yg = np.abs(np.fft.fftshift(Fg)) #this shifts the data into the order we'd expect/want for graphing
    
    xt = np.fft.fftshift(np.fft.fftfreq(n, d = timestep)) #finds the frequancies for the fourier transform plot
    yt = np.abs(np.fft.fftshift(Ft)) #this shifts the data into the order we'd expect/want for graphing
    
    fc = Ft*Fg #applies convolution theorem
    c1 = np.fft.ifft(fc) #takes the inverse FFT
    c = np.fft.fftshift(c1).real #this shifts the data into the order we'd expect/want for graphing

    return (xg,yg,xt,yt, c)

'''
The following code plots the results for the questions as well as the true convolution
to compare my results against.
'''

t = sp.linspace(-20,20,41) #the selected amount of sampling
xg,yg,xt,yt,c = conv(t) #run the convolution function
fig, axs = plt.subplots(4) #create subplots
fig.tight_layout(pad=3.0)


axs[2].plot(t,c, color = 'forestgreen', label = 'Convolution (g*h)(t)') #plot the convolution I calculated
axs[2].set_xlim(-12,12)
axs[2].set_xlabel('t (s)')
axs[2].set_title('Convolution with Time Step = 0.1s')

t = sp.linspace(-20,20,4001) #the selected amount of sampling
xg,yg,xt,yt,c = conv(t) #run the convolution function

axs[0].plot(t,gauss(t), color = 'cornflowerblue', label = 'Gauss Function') #plot the gauss
axs[0].plot(t,tophat(t), color = 'darkorange', label = 'Rect Function') #plot the tophat function
axs[0].set_xlabel('t (s)')
axs[0].legend(loc='upper left', frameon=False)
axs[0].set_xlim(-12,12)
axs[0].set_title('Input Functions')

axs[1].plot(xg,yg, color = 'cornflowerblue', label = 'Modulus of Fourier Transfrom of g') #plot the fourier transform of the gaussian
axs[1].plot(xt,yt, color = 'darkorange', label = 'Modulus of Fourier Transfrom of h') #plot the fourier transform of the tophat
axs[1].set_xlim(-3.6,3.6)
axs[1].legend(loc='upper left', frameon=False)
axs[1].set_xlabel('ω (s$^-$$^1$)')
axs[1].set_title('Modulus of Fourier Transforms')

axs[3].plot(t,c, color = 'forestgreen', label = 'Convolution (g*h)(t)') #plot the convolution I calculated
axs[3].set_xlim(-12,12)
axs[3].set_xlabel('t (s)')
axs[3].set_title('Convolution with Time Step = 0.01s')
#I now plot the built in convolution with very high sampling rate so I can assume it is the true value
t = sp.linspace(-20,20,40001) #very high sampling rate so we can assume this convolution is the true value
t1 = sp.linspace(-40,40,80001) #t values to plot the built in convolution against
y = spf.fftconvolve(gauss(t),tophat(t))/len(t) #find the built in convolution 
axs[3].plot(t1,y, '--', dashes = (5,15), color = 'red', label = 'True Convolution') #plot the convolution built in
axs[3].legend(loc='upper left', frameon=False)
axs[2].plot(t1,y, '--', dashes = (5,15), color = 'red', label = 'True Convolution') #plot the convolution built in
axs[2].legend(loc='upper left', frameon=False)
plt.show()



