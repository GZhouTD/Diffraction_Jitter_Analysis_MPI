#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:32:28 2020

@author: gzhou
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close()
#basic parameters
ph_E = 7.5 #keV
c_const = 12.39842017144820
lr = c_const/ph_E #0.1 nm
lr = lr /1e10
c_speed = 299792458
c_f = c_speed/lr
r_sig_f = 5e-4
sig_f = c_f*r_sig_f


# transverse profile
sig_x = 20 # um
m_x = 5
x_l = -m_x*sig_x# um
x_r =  m_x*sig_x # um

n_x = 100 # x samples
x = np.linspace(x_l,x_r,n_x)

trans = np.exp(-x**2/2/sig_x**2)
plt.figure()
plt.subplot(1,2,1)
plt.plot(x,trans)


# longitudinal Guassian

n_f = 8192*2
m_f = 3 # 3 sigma range
f_l = c_f-m_f*sig_f
f_r = c_f+m_f*sig_f
freq = np.linspace(f_l,f_r,n_f)

ltd = np.exp(-(freq-c_f)**2/2/sig_f**2)
plt.subplot(1,2,2)

plt.plot(freq,ltd)

[f_mesh, x_mesh]=np.meshgrid(freq,x)

nor_profile = np.exp(-(f_mesh-c_f)**2/2/sig_f**2)*np.exp(-x_mesh**2/2/sig_x**2)
#plt.figure()
#plt.contourf(x_mesh,f_mesh,nor_profile)
plt.figure()
plt.contourf(f_mesh,x_mesh,nor_profile)

e_field = np.sqrt(nor_profile)
efield_r = np.real(e_field)
efield_i = np.imag(e_field)

np.savetxt('laser_real.in',efield_r)
np.savetxt('laser_imag.in',efield_i)
np.savetxt('x.in',x)
np.savetxt('freq.in',freq)

#r =np.loadtxt('out.dat')
#
eph = 12.386/(c_speed/freq*1e10)
#plt.plot(eph*1e3-7500,r[:,1])

plt.plot(eph*1e3-7500,ltd)




##shape profile
gl =-1000 #um
gr = 1000 #um
n_g = 100000
g_x = np.linspace(gl,gr,n_g)
sig_gx1 = 200/2.355
sig_gx2 = 1600/2.355


g_h = 0#0.0012*np.exp(-g_x**2/2/sig_gx1**2)+0.003*np.exp(-g_x**2/2/sig_gx2**2)

plt.plot(g_x,g_h)

shape_data = np.zeros([n_g,4])
shape_data[:,0]=g_x
shape_data[:,1]=g_h
#np.ones([n_g])*0
#g_h
shape_data[:,2]= np.ones([n_g])*0
shape_data[:,3]=  np.ones([n_g])*0#*2.66e-5

np.savetxt('shape.in',shape_data)






#out=np.abs(np.sum(e_field,axis=0)*(r[:,5]+1j*r[:,6]))**2
#plt.plot(out)
#out2 = np.zeros([500,8192])*1j
#for i in range(500):
#    out2[i,:] = e_field[i,:]*(r[:,5]+np.random.randn(8192)*1+1j*(r[:,6]+np.random.randn(8192)*1))
#
#plt.plot(np.abs(np.sum(out2,axis=0))**2)






















