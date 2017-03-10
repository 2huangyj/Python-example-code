
# coding: utf-8

# In[1]:

import csv
import numpy as np 
import matplotlib.pyplot as plt  
import scipy
from scipy.signal import find_peaks_cwt
import peakutils
from scipy.optimize import curve_fit


# In[26]:

lower_lim_1 = 1.7
upper_lim_1 = 2.0
lower_lim_2 = 3.5
upper_lim_2 = 4.0
time_cutoff = 80
R=20
power = "200mw"


# In[56]:

f = open(power+".dat")
for row in csv.reader(f):#this is picking out each row 
	try:
		numbers = row[0].split('\t')
		new_numbers = [float(number) for number in numbers]
		data_array = np.row_stack((data_array,np.array(new_numbers)))
 	except ValueError:#I did tihs for the command "new_numbers"
	 	print row
	 	continue
	except NameError:#I did this for the command"data_array"
	#, but both exception handling are for the first row!
 		data_array = np.array(new_numbers)


# In[57]:

time = data_array[:,0]#np.shape(time)
reflectivity = data_array[:,1]


# In[58]:

fig1 = plt.figure()
plt.plot(time, reflectivity/R, color = 'red', label = 'reflectivity versus time')
plt.savefig("trace for "+power+".png")

time_sub = time[time_cutoff:500]
reflectivity_sub = reflectivity[time_cutoff:500]


# In[45]:

R_fft = np.fft.fft(reflectivity_sub)
space = (time_sub[0]-time_sub[-1])/float(time_sub.shape[-1]-1)#time.shape[-1] = 500L,type is long
f = np.fft.fftfreq(time_sub.shape[-1],space)

fig2 = plt.figure()
plt.plot(f[250:500],np.abs(R_fft)[250:500],color = 'red', label = 'FFT of reflectivity')
plt.savefig("FFT for "+power+".png")


# In[46]:

cb_w = np.abs(R_fft)
indexes_w = find_peaks_cwt(cb_w,np.arange(5,15))
l = list(f[indexes_w])
l.sort(reverse = False)
print l

mode1 = 0
mode2 = 0
for peak in l:
    if (peak<upper_lim_1 and peak>lower_lim_1):
        mode1 = peak
        print "mode1 is"
        print mode1
    elif (peak<upper_lim_2 and peak>lower_lim_2):
        mode2 = peak
        print "mode2 is"
        print mode2


# In[47]:

def func(x, A,B,C,tau0,tau1,offset):
    return A*np.exp(-x/tau0)+(B*np.sin(mode1*x)+C*np.cos(mode1*x))* np.exp(-x/tau1) + offset


# In[48]:

popt,pcov = curve_fit(func,time_sub,reflectivity_sub)


# In[50]:

(A, B, C, tau0,tau1,offset)=popt
delta_0 = (np.average(reflectivity[0:15])-A-offset)/R 



# In[ ]:
with open("result.txt",'a') as f:
    f.write('power:{}\n'.format(power))
    f.write('mode1:{}\n'.format(mode1))
    f.write('mode2:{}\n'.format(mode2))
    f.write('delta_0:{}\n'.format(delta_0))
    f.write('time_cutoff:{}\n'.format(time_cutoff))
