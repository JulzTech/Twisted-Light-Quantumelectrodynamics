# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 23:31:50 2021

@author: JulzTech
"""

# Imports, initialisations, and physical constants



import numpy as np

import matplotlib.pyplot as plt



iunit = (0+1j)

hbar = 1.054e-34 #SI

ke = 8.99e9 #SI

lightspeed = 299792458 #SI





# Functions



def delta(i,j):

    if i==j:

        return 1.0

    else:

        return 0.0



def alpha_D(i,j):

    return mu_00D[i]*mu_0bD[j]/(e_gap+iunit*damping_D0)\
        +mu_00D[j]*mu_0bD[i]/(e_gap-e_photon+iunit*damping_D0)\
            +mu_0bD[i]*mu_bbD[j]/(iunit*damping_Db)\
                +mu_0bD[j]*mu_bbD[i]/(-e_photon+iunit*damping_Db)



def alpha_A(i,j):

    return mu_b0A[i]*mu_00A[j]/(e_photon+iunit*damping_A0)\
        +mu_b0A[j]*mu_00A[i]/(iunit*damping_A0)\
            +mu_bbA[i]*mu_b0A[j]/(-e_gap+e_photon+iunit*damping_Ab)\
                +mu_bbA[j]*mu_b0A[i]/(-e_gap+iunit*damping_Ab)



def alpha_T(i,j):

    return mu_00T[i]*mu_00T[j]/(e_photon+iunit*damping_T0)\
        +mu_00T[j]*mu_00T[i]/(-e_photon+iunit*damping_T0)\
            +mu_0bT[i]*mu_b0T[j]/(-e_gap+e_photon+iunit*damping_Tb)\
                +mu_0bT[j]*mu_b0T[i]/(-e_gap-e_photon+iunit*damping_Tb)



def eta_DA(m,i,j):

    return delta(i,j)-m*R_DA[i]*R_DA[j]*R_DAs**-2

def eta_DT(m,i,j):

    return delta(i,j)-m*R_DT[i]*R_DT[j]*R_DTs**-2

def eta_AT(m,i,j):

    return delta(i,j)-m*R_AT[i]*R_AT[j]*R_ATs**-2



def V_kDA(i,j):

    return ((ref_index**2+2)/(3*ref_index))**2*R_DAs**-3*ke**-1\
        *np.exp(iunit*kw*R_DAs*ref_index.real-kw*R_DAs*ref_index.imag)\
            *(eta_DA(3,i,j)+kw*R_DAs*ref_index.imag*eta_DA(3,i,j)\
              +kw**2*R_DAs**2*(ref_index.imag**2-ref_index.real**2)*eta_DA(1,i,j)\
                  -iunit*kw*R_DAs*ref_index.real*eta_DA(3,i,j)-2*iunit*kw**2*R_DAs**2*ref_index.real*ref_index.imag*eta_DA(1,i,j))



def V_0DT(i,j):

    return ((ref_index**2+2)/(3*ref_index))**2*R_DAs**-3*ke**-1*eta_DT(3,i,j)



def V_kDT(i,j):

    return ((ref_index**2+2)/(3*ref_index))**2*R_DTs**-3*ke**-1\
        *np.exp(iunit*kw*R_DTs*ref_index.real-kw*R_DTs*ref_index.imag)\
            *(eta_DT(3,i,j)+kw*R_DTs*ref_index.imag*eta_DT(3,i,j)\
              +kw**2*R_DTs**2*(ref_index.imag**2-ref_index.real**2)*eta_DT(1,i,j)\
                  -iunit*kw*R_DTs*ref_index.real*eta_DT(3,i,j)-2*iunit*kw**2*R_DTs**2*ref_index.real*ref_index.imag*eta_DT(1,i,j))



def V_0AT(i,j):

    return ((ref_index**2+2)/(3*ref_index))**2*R_ATs**-3*ke**-1*eta_AT(3,i,j)



def V_kAT(i,j):

    return ((ref_index**2+2)/(3*ref_index))**2*R_ATs**-3*ke**-1\
        *np.exp(iunit*kw*R_ATs*ref_index.real-kw*R_ATs*ref_index.imag)\
            *(eta_AT(3,i,j)+kw*R_ATs*ref_index.imag*eta_AT(3,i,j)\
              +kw**2*R_ATs**2*(ref_index.imag**2-ref_index.real**2)*eta_AT(1,i,j)\
                  -iunit*kw*R_ATs*ref_index.real*eta_AT(3,i,j)-2*iunit*kw**2*R_ATs**2*ref_index.real*ref_index.imag*eta_AT(1,i,j))





# Parameters



fermi_pre = 1.0  #2*pi*rho/hbar

e_gap = 6.32e-19 #J transition energy

ref_index = 1.1+iunit*0.1



mu_base = 3.33e-29 #Cm

mu_00D = [mu_base, 0.0, 0.0]

mu_0bD = [mu_base, 0.0, 0.0]

mu_b0D = [mu_base, 0.0, 0.0]

mu_bbD = [mu_base, 0.0, 0.0]



mu_00A = [mu_base*0.707, mu_base*0.707, 0.0]

mu_0bA = [mu_base*0.707, mu_base*0.707, 0.0]

mu_b0A = [mu_base*0.707, mu_base*0.707, 0.0]

mu_bbA = [mu_base*0.707, mu_base*0.707, 0.0]



damping_D0 = 1.5e-20*e_gap

damping_Db = 1.5e-20*e_gap

damping_A0 = 1.5e-20*e_gap

damping_Ab = 1.5e-20*e_gap

damping_T0 = 1.5e-20*e_gap

damping_Tb = 1.5e-20*e_gap



R_DA = [0.0, 0.0, 1.4e-9]

R_DT = [0.0, 0.0, 0.7e-9]





# Calculations



R_AT = [R_DA[0]-R_DT[0], R_DA[1]-R_DT[1], R_DA[2]-R_DT[2]]

R_DAs = np.sqrt(R_DA[0]**2+R_DA[1]**2+R_DA[2]**2)

R_DTs = np.sqrt(R_DT[0]**2+R_DT[1]**2+R_DT[2]**2)

R_ATs = np.sqrt(R_AT[0]**2+R_AT[1]**2+R_AT[2]**2)



e_photon = e_gap #Exact resonance

kw = e_photon/(hbar*lightspeed)



abscissa_list = []

da_da = []

da_tda = []

da_dat = []

da_dta = []

tda_tda = []

tda_dat = []

tda_dta = []

dat_dat = []

dat_dta = []

dta_dta = []



for iteration in range(0,200):



    theta = 6.2824*iteration/200.0

    mu_00T = [mu_base*np.sin(theta), mu_base*np.cos(theta), 0.0]

    mu_0bT = [mu_base*np.sin(theta), mu_base*np.cos(theta), 0.0]

    mu_b0T = [mu_base*np.sin(theta), mu_base*np.cos(theta), 0.0]

    mu_bbT = [mu_base*np.sin(theta), mu_base*np.cos(theta), 0.0]



    mDA_list = []

    for i in range(3):

        for j in range(3):

            mDA_list.append(mu_0bD[i] * V_kDA(i,j) * mu_b0A[j])

    M_DA = sum(mDA_list)



    mTDA_list = []

    for i in range(3):

        for j in range(3):

            for k in range(3):

                for l in range(3):

                    mTDA_list.append(mu_00T[i] * V_0DT(i,j) * alpha_D(j,k) * V_kDA(k,l) * mu_b0A[l])

    M_TDA = sum(mTDA_list)



    mDAT_list = []

    for i in range(3):

        for j in range(3):

            for k in range(3):

                for l in range(3):

                    mDAT_list.append(mu_0bD[i] * V_kDA(i,j) * alpha_A(j,k) * V_0AT(k,l) * mu_00T[l])

    M_DAT = sum(mDAT_list)



    mDTA_list = []

    for i in range(3):

        for j in range(3):

            for k in range(3):

                for l in range(3):

                    mDTA_list.append(mu_0bD[i] * V_kDT(i,j) * alpha_T(j,k) * V_kAT(k,l) * mu_b0A[l])

    M_DTA = sum(mDTA_list)



    abscissa_list.append(theta)

    da_da.append(1.0)

    da_tda.append(2*(M_DA*(M_TDA.real-iunit*M_TDA.imag)).real/abs(M_DA)**2)

    da_dat.append(2*(M_DA*(M_DAT.real-iunit*M_DAT.imag)).real/abs(M_DA)**2)

    da_dta.append(2*(M_DA*(M_DTA.real-iunit*M_DTA.imag)).real/abs(M_DA)**2)

    tda_tda.append(abs(M_TDA)**2/abs(M_DA)**2)

    tda_dat.append(2*(M_TDA*(M_DAT.real-iunit*M_DAT.imag)).real/abs(M_DA)**2)

    tda_dta.append(2*(M_TDA*(M_DTA.real-iunit*M_DTA.imag)).real/abs(M_DA)**2)

    dat_dat.append(abs(M_DAT)**2/abs(M_DA)**2)

    dat_dta.append(2*(M_DAT*(M_DTA.real-iunit*M_DTA.imag)).real/abs(M_DA)**2)

    dta_dta.append(abs(M_DTA)**2/abs(M_DA)**2)





plt.rc('font', size=12)          # controls default text sizes

plt.rc('axes', titlesize=12,labelsize=8,linewidth=0.5)     # fontsize of the axes title,labels and linewidth for subplots

plt.rc('xtick', labelsize=10)    # fontsize of the tick labels

plt.rc('ytick', labelsize=10)    # fontsize of the tick labels

plt.rc('legend', fontsize=10)

plt.rc('lines', markersize=2,markeredgewidth=0.2,linewidth=1.0)

plt.xlabel(r'$\theta$/radians',fontsize=12)

plt.ylabel('Total Rate',fontsize=12)







#plt.plot(abscissa_list,da_da,label='DA pure')

#plt.plot(abscissa_list,da_tda,label='DA:TDA')

#plt.plot(abscissa_list,da_dat,label='DA:DAT')

#plt.plot(abscissa_list,da_dta,label='DA:DTA')

#plt.plot(abscissa_list,tda_tda,label='TDA pure')

#plt.plot(abscissa_list,tda_dat,label='TDA:DAT')

#plt.plot(abscissa_list,tda_dta,label='TDA:DTA')

#plt.plot(abscissa_list,dat_dat,label='DAT pure')

#plt.plot(abscissa_list,dat_dta,label='DAT:DTA')

#plt.plot(abscissa_list,dta_dta,label='DTA pure')



total_rate = []

for i in range(0,200):

    total_rate.append(da_da[i]+da_tda[i]+da_dat[i]+da_dta[i]+tda_tda[i]+tda_dat[i]+tda_dta[i]+dat_dat[i]+dat_dta[i]+dta_dta[i])

plt.plot(abscissa_list,total_rate,label='Total rate')



#plt.legend()

plt.show()



plt.savefig('./fig_1b.png')
