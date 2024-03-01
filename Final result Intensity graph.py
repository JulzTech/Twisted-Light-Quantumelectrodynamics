# -*- coding: utf-8 -*-
"""
Created on Thu May 20 13:11:37 2021

@author: JulzTech
"""

# Import math functions

import numpy as np

import matplotlib.pyplot as plt



# Complex conjugate using Numpy function

def bar(num):

    return np.conj(num)


# Defined constant Clp that we will need for later

def clp(p,l):

    return np.sqrt((2*(np.math.factorial(p)/(np.pi*np.math.factorial((p+abs(l)))))))


# Levi Civita

T = np.array([[[0,0,0],[0,0,1],[0,-1,0]],[[0,0,-1],[0,0,0],[1,0,0]],[[0,1,0],[1,0,0],[0,0,0]]])

newarr1 = T.reshape(3, 3, 3)

E=newarr1          


#Making a grid in x,y (projecting at defied z)

x = np.linspace(-3.0e-6,3.0e-6,1001)

x[500] = 1e-16 #Remove r=0

y = np.linspace(-3.0e-6,3.0e-6,1001)  

y[500] = 1e-16 #Remove r=0  

X, Y = np.meshgrid(x,y)


#Z = np.linspace(-3.0e-6,3.0e-6,1001)

 

# Input parameters to be altered as necessary

k=((2*np.pi)/(729e-9))


r=np.sqrt(x**2+y**2)

R=np.sqrt(X**2+Y**2)

w0 = 2*(729e-9)


z = k*w0


print (z)


iunit = (0+1j)

hbar = 1.054e-34 #SI

ke = 8.99e9 #SI

c = 299792458


V=1

w=1

n=1

N=1

rho = 1

r=1


mu_base = 3.33e-29 #Cm


def I(w):

    return (rho*2*((np.pi)*n*c*hbar*w)/V)

   

mu_f0 = [mu_base, 0.0, 0.0]


#Same thing for e???


e0 =  8.854187817e-12

km = [0,0,1]

ell = -1


#def phi[i,k,m]:

#    angle=np.arctan(Y/X)

#    return [-np.sin(phi),np.cos(phi),0]


xx = 0.5 

yy = 0.5

zz = -0.5 


xy=yx=0.001

yz=zy=0.001

xz=zx=0.001


arr = np.array([[xx,xy,xz],[yx,yy,yz],[zx,zy,zz]])

newarr = arr.reshape(3, 3)

Q=newarr

   


#other way around

#quadrupole symmetric and traceless


print (Q)

# Define our equation (f)**2

def f(X,Y):

    return ((clp(0,1)/w0))**2*((np.sqrt(2*R)/w0))**2\
        *(np.exp(-(R**2/(w0**2)))**2)*((((2*R)**2)/(w0**2))**2)


#Plot Fermi

for i in range(3):

        for k in range(3):

            for m in range(3) :

                print(i,k,m)

                print(E[i,k,m])



# First Fermi contribution  

def Fermi_1(X,Y):

    i=0

    k=1

    m=2

    phi=np.arctan2(Y,X) 

    return (1/R)*((f(X,Y))*E[i,k,m]*km[m]*(-np.sin(phi))*ell*((mu_f0[k]*(Q[i,i]))-(mu_f0[i]*Q[k,i])))

# Second Fermi contribution with phi j as y    
def Fermi_2(X,Y):

    i=1

    k=0

    m=2

    phi=np.arctan2(Y,X) 

    return (1/R)*((f(X,Y))*E[i,k,m]*km[m]*(np.cos(phi))*ell*((mu_f0[k]*(Q[i,k]))-(mu_f0[k]*Q[k,k])))







# Phi?
# Order of E from screenshot for i and k - 12 more equations with the variables in different order?


# plot



plt1=plt


#plt1.plot(R*1e6,Fermi_2(X,Y))



contour_levels = np.linspace(0,1,1001)

clrs = 'viridis' #Colormap 'viridis' or 'plasma'


plt2=plt

plt2.figure(figsize=(8,6))

plt2.imshow(Fermi_1(X,Y),interpolation='nearest',extent=[-2.0,2.0,-2.0,2.0])

plt2.xlabel(r'$\mathrm{X\,/\mathrm{{\mu}m}}$')

plt2.ylabel(r'$\mathrm{Y\,/\mathrm{{\mu}m}}$')

plt2.colorbar()


plt.show()
