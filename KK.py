#!/home/coffee/anaconda3/bin/python
import numpy as np;
import matplotlib.pyplot as plt;
import scipy
from scipy import special
import scipy.special as special
import scipy.integrate as integrate

#laser pulse
def gauss(x,x0,xw):
    return np.exp(-((x-x0)/(xw))**2);
# Define time domain vector for intensity and one for phase

nsamples=10000#4000
t0=200;
wt=20;#5
centerfreq=.5;
#i=np.arange(nsamples);
t = np.linspace(0,1000,nsamples);
aoft = np.zeros(nsamples);
poft = np.cos(centerfreq*2.*np.pi*t);

def addquadratic(value,freq,f0):
    return (value * np.sign(freq) * (abs(freq)-f0)**2);

aoft = aoft + gauss(t,t0,wt);
zoft = aoft*poft;
zofw = np.fft.fft(zoft.real); # electric field is the real part... 
f = np.fft.fftfreq(t.shape[-1],0.1);#frequency domain 0.25
aofw = np.abs(zofw);
pofw = np.angle(zofw);
pofw = pofw + addquadratic(500,f,centerfreq);#1000

#ifft into the time domain, after adding chirp to the laser pulse
chirped = aofw*np.exp(1j*pofw);
chirpedtime = np.fft.ifft(chirped);

# Plot the unchirped and chirped pulse
fig1=plt.figure();
plt.plot(t,zoft.real,color='blue',label='zoft.real');
plt.plot(t,chirpedtime.real,color='green',label='chirpedtime.real');
plt.show();
chirpedw = np.fft.fft(chirpedtime);
#plot the argument of unchirped and chirped pulse(in frequency domain)
fig2=plt.figure();
plt.plot(f,np.abs(zofw),color='black',label='zofw amp)');
plt.scatter(f,np.abs(chirpedw),color='blue',label='zofw amp)');
plt.show();
#plot the amplitude of unchirped and chirped pulse(in frequency domain)
fig3=plt.figure();
plt.scatter(f,np.unwrap(np.angle(zofw)),color='red',label='zofw arg');
plt.scatter(f,np.unwrap(np.angle(chirpedw)),color='green',label='zofw arg');
plt.show();


#def integrand(x,w,w0,c):
   # return (absorption(x,w0,c)*x/((x**2-w**2)+0.0000000000001))#
# refractive(imag_epsilon,f):
 #   return np.fft.fft(-1j*np.sign(f)/np.pi*np.fft.ifft(imag_epsilon))

gamma = 0.05
f0 = 3#resonance frequency
t_rise = 1
t_fall = 100
win=50.;
tau = 1;#2*pi*d/c

absrp_spec = None
refr_spec = None
density = None
ndelays = 1000;
f_ = list(f)
f_modified = f_[:]
for i in range(len(f_modified)):
    if f_modified[i]<0:
        f_modified[i] = f_modified[i]+4


def absorption(w,w0,g):
    #global gamma
    #g = gamma + 0.3*n**2
    return 0.01/(g**2+((w-w0)**2))
    '''if w > 0:
        return 1./(g**2+((w-w0)**2))
    elif w < 0:
        return (-1.)/(g**2+((w-w0)**2))
    else:
        return 0'''

def carrier_density(t0_laser,t0_xray):
    global t_rise,t_fall
    y1 = (1+special.erf((t0_laser-t0_xray)/t_rise))/2.0
    y2 = (1-special.erf((t0_laser-t0_xray)/t_fall))/2.0
    return y1*y2

#def switch_reso_freq(t):
#    global f0
#    return f0 - (carrier_density()

t0_vec=np.linspace(t0-win,t0+win,ndelays);
for i in range(len(t0_vec)):
    density = carrier_density(t0,t0_vec[i])
    f0_dyn = f0 - 4*density**2
    g = gamma + 0.3*density**2
    absrp_spec = np.array([absorption(fx,f0_dyn,g) for fx in f_])
    #absrp_spec = [absorption(fx,f0) for fx in f_]
    #print "i"
    refr_spec = np.fft.fft(-1j*np.sign(f)/np.pi*np.fft.ifft(absrp_spec))
    #refr_spec = [refractive(fx,f0,c) for fx in f_]
    #print "r"
    aofw = aofw*np.exp(-np.array(f_modified)*tau*absrp_spec); #np.abs(chirpedtime)*(1-.5*gauss(t,t0_vec[i],wt));
    pofw = pofw + np.array(f_modified)*tau*refr_spec; #np.angle(chirpedtime)-erf(t,t0_vec[i],wt,np.pi/20);#gauss(t,t0,wt)*h_phase
    #pofc=np.angle(chirpedtime)+erf(t,delay,wt,np.pi/20);#gauss(t,t0,wt)*h_phase
    #aofc=np.abs(chirpedtime)*(1-.5*gauss(t,delay,wt));

    new_F = aofw*np.exp(1j*pofw)
    #new_T = np.fft.ifft(new_F)
    if i%200==0:
        #print zip(absrp_spec,f)
        #print zip(refr_spec,f)
        fig4=plt.figure();
        plt.plot(f,absrp_spec)
        plt.show()
        fig5=plt.figure();
        plt.plot(f,refr_spec)
        plt.show()
        fig6=plt.figure();
        plt.plot(f,aofw)
        plt.show()
        fig7=plt.figure();
        plt.plot(f,pofw)
        plt.show()
        new_T = np.fft.ifft(new_F)
        fig8=plt.figure();
        plt.plot(f,np.abs(new_F))
        plt.show()

    part = new_F[0:len(new_F)//2]
    if i==0:
        R=np.array(np.abs(part))
    else:
        R=np.row_stack((np.abs(part),R))


plt.imshow(R)
plt.hot()
plt.colorbar()
plt.show()

