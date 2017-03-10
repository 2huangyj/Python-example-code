import numpy as np
import matplotlib.pyplot as plt 
import random
from scipy.optimize import curve_fit
import csv

#cluster Monte Carlo
D = 3
l = 20
t_step = 7
J = 1
T_list = [0.5,1.0,2.0,3.0,3.5,4.0,4.2,4.4,4.5,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0]

#T_list = [1.0,1.25,1.5,1.75,2.0,2.1,2.15,2.20,2.24,2.25,2.26,2.27,2.28,2.29,2.30,2.4,2.5,2.75,3.0,3.25,3.75,4.0,4.25,4.5,4.75]
#T = 2.5
M_list = []
C_list = []
chi_list = []
ksai_list = []
#create a stack
class Stack:
	def __init__(self):
		self.items = []

	def isEmpty(self):
		return self.items == []

	def push(self,item):
		self.items.append(item)

	def pop(self):
		return self.items.pop(0)
def func(x,ksai,A,Const):
	return A*np.exp(-x/ksai)+Const

for T in T_list:
	print 'T:', T
	t = 0
	t_max = 980#100000
	beta = 1/T
	p = 1-np.exp(-2*beta*J)
	grid = 2*np.random.randint(2, size=(l,l,l))-1

	s = Stack()
	sum_M = 0.0
	sum_E = 0.0
	sum_M_2 = 0.0
	sum_E_2 = 0.0
	correlation = np.zeros(shape = (l,))

	#flip the spin according to the algorithm
	while t < t_max:
		flipped = np.zeros(shape=(l,l,l))
		i = np.random.randint(l - 1)
		j = np.random.randint(l - 1)
		k = np.random.randint(l - 1)
		s.push((i,j,k))
		flipped[i,j,k] = 1
		
		M = 0.0
		E = 0.0

		while not s.isEmpty():
			x,y,z = s.pop()
			for (x_,y_,z_) in [((x-1)%l,y,z),(x,(y-1)%l,z),(x,y,(z-1)%l),((x+1)%l,y,z),(x,(y+1)%l,z),(x,y,(z+1)%l)]:
				if (grid[x,y,z] == grid[x_,y_,z_]) and (flipped[x_,y_,z_]==0) and (p>random.random()):
					s.push((x_,y_,z_))
					flipped[x_,y_,z_]=1
			grid[x,y,z] = -grid[x,y,z] 
			#inner while loop is over


		if (t+1)%t_step == 0:
			print t
			#calculate M
			M = sum(sum(sum(grid)))
			M = M/(l**D*1.0)
			sum_M = sum_M + np.abs(M)
			sum_M_2 = sum_M_2 + M * M

			#calculate E		
			for i in range(l):
				for j in range(l):
					for k in range(l):
						E = E - J*(grid[i,j,k]*grid[i,(j+1)%l,k] + grid[i,j,k]*grid[(i+1)%l,j,k] + grid[i,j,k] * grid[i,j,(k+1)%l])
						
			E = E/(l**D*1.0)
			sum_E = sum_E + np.abs(E)
			sum_E_2 = sum_E_2 + E * E

			#calculate the correlation length
			cross  = np.zeros(shape = (l,))
			for i in range(l):
				for j in range(l):
					grid_piece_k = (np.abs(np.fft.fft(grid[:,i,j])))**2
					cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))
					grid_piece_k = (np.abs(np.fft.fft(grid[i,:,j])))**2
					cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))
					grid_piece_k = (np.abs(np.fft.fft(grid[i,j,:])))**2
					cross = np.vstack((cross,np.real(np.fft.ifft(grid_piece_k))))

			cross = np.sum(cross,axis = 0)/(1.0*D*l**D) 
			correlation = correlation + cross
		

		#print t
		t = t + 1
		#outer while loop is over

	ave_M = sum_M/((t_max//t_step)*1.0)
	ave_M_2 = sum_M_2/((t_max//t_step)*1.0)

	ave_E = sum_E/((t_max//t_step)*1.0)
	ave_E_2 = sum_E_2/((t_max//t_step)*1.0)

	correlation = (correlation/((t_max//t_step)*1.0))[0:10]
	popt,pcov = curve_fit(func,np.arange(0,10),correlation)
	ksai = popt[0]
	offset = popt[2]

	C = ave_E_2-ave_E**2
	chi = ave_M_2-ave_M**2

	if T == 4.5:
		fig = plt.figure()
		plt.plot(np.arange(len(correlation)),correlation-offset)
		plt.savefig("correlation @ " + str(T) +".png")
		fig.clear()
#T_list
	M_list.append(ave_M)
	C_list.append(C)
	chi_list.append(chi)
	ksai_list.append(ksai)
	
with open('result.csv','w') as csvfile:
	fieldnames = ['T','M','C','chi','ksai']
	writer = csv.DictWriter(csvfile,fieldnames = fieldnames)

	writer.writeheader()
	for i in range(len(T_list)):
		writer.writerow({'T':T_list[i],'M':M_list[i],'C':C_list[i],'chi':chi_list[i],'ksai':ksai_list[i]})
