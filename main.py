# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

g = 9.81
dm = 0.5 * 1e-6 #+/-10mg
dI = 0.01
dH = 2 * 1e-3

def quad_mono(x, a):
    return a * x**2

def curie_fit(T, acm, theta):
    return acm / (T-theta)

# a very unclean way to to do a linear interpolation from I to H from the "Linearity between H and I" measurement:
def I_to_H(I):
    N = len(I)
    n = len(I_lin)
    H = np.zeros(N)-1
    for i in range(N):
        for j in range(n):
            if(I[i] < I_lin[j]):
                dist = I_lin[j] - I_lin[j-1]
                part = I[i] - I_lin[j-1]
                pd = part/dist
                if(H[i] == -1):
                    H[i] = pd * H_lin[j] + (1-pd)* H_lin[j-1]
            elif(I[i] == 10):
                H[i] = H_lin[n-1]
    return H

#plt.errorbar(fmt="none", color='b', alpha=0.3, capsize=2)


#Maximize Force as function of height, I = 6 A
l_height = np.array([3, 5, 5.3, 7.8, 10.8, 12.8])
K_height = np.array([63110-63090, 63155-63118, 63379-63335,  63724-63495, 64397-63650, 64462-63880])*1e-6*g #mg->kg->N
dl_height = 0.1
dK_height = dm*g

plt.errorbar(l_height, K_height*1e3, xerr=dl_height, yerr=dK_height, fmt='x', capsize=2)
plt.grid()
plt.title("Force depending on the height of the sample")
plt.xlabel("Height (cm)")
plt.ylabel("Downforce (mN)")
plt.savefig('height.pdf', format='pdf')
plt.show()


#Linearity between H and I
H_lin = np.array([0.8, 17.96, 35.9, 53.8, 71.4, 88.8, 105.9, 122.4, 137.7, 150.7, 162.5])*1e-3 #Tesla
I_lin = np.array([0, 1, 2, 3, 4.01, 5, 5.99, 7, 8, 9, 10])

num=100000
I_test = np.linspace(0,10,num=num)
H_test = I_to_H(I_test)

plt.plot(I_test, H_test*1e3, color='g', label="linear interpolation", linewidth=0.7)
plt.errorbar(I_lin, H_lin*1e3, xerr=dI, yerr=dH, fmt='x', capsize=2, label="measurement")
plt.grid()
plt.title("Magnetic field strength depending on the current")
plt.xlabel("Current through the electromagnet (A)")
plt.ylabel("Magnetic field strength (mT)")
plt.legend()
plt.savefig('HI.pdf', format='pdf')
plt.show()


#Current at Room Temperature (26.5 Celsius)
I_RT = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.51, 4.01, 4.51, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
m0_RT = 63350 #Leergewicht in mg
m_RT = np.array([63351, 63370, 63395, 63428, 63474, 63529, 63590, 63661, 63742, 63829, 63927, 64034, 64159, 64267, 64386, 64500, 64618, 64735, 64851, 64965])
K_RT = (m_RT - m0_RT)*1e-6*g
H_RT = I_to_H(I_RT)
dK_RT = g * np.sqrt((dm*m0_RT)**2 + (dm*m_RT)**2)

fit_RT, covfit_RT = curve_fit(quad_mono, H_RT, K_RT)
dfit_RT = np.sqrt(covfit_RT[0])
print("fit_RT: " + str(fit_RT[0]) + " +/- " + str(dfit_RT[0]))

plt.errorbar(H_RT*1e3, K_RT*1e3, xerr=dH, yerr=dK_RT, fmt='x', capsize=2, label="measurement")
plt.plot(H_RT*1e3, quad_mono(H_RT, fit_RT)*1e3, label="quadratic fit", color='r')
plt.grid()
plt.title("Force depending on the magnetic field (T=$26.5^\circ$C)")
plt.legend()
plt.xlabel("Magnetic field strength (mT)")
plt.ylabel("Downforce (mN)")
plt.savefig('RT.pdf', format='pdf')
plt.show()


#Current at Room Temperature (20.5 Celsius) (no heater)
I_RT2 = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.51, 4.01, 4.51, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
m0_RT2 = 63350 #Leergewicht in mg
m_RT2 = np.array([63356, 63370, 63395, 63431, 63474, 63528, 63592, 63665, 63747, 63839, 63937, 64044, 64161, 64289, 64400, 64521, 64648, 64764, 64880, 64995])
K_RT2 = (m_RT2 - m0_RT2)*1e-6*g
H_RT2 = I_to_H(I_RT2)
dK_RT2 = g * np.sqrt((dm*m0_RT2)**2 + (dm*m_RT2)**2)

fit_RT2, covfit_RT2 = curve_fit(quad_mono, H_RT2, K_RT2)
dfit_RT2 = np.sqrt(covfit_RT2[0])
print("fit_RT2: " + str(fit_RT2[0]) + " +/- " + str(dfit_RT2[0]))

plt.errorbar(H_RT2*1e3, K_RT2*1e3, xerr=dH, yerr=dK_RT2, fmt='x', capsize=2, label="measurement")
plt.plot(H_RT2*1e3, quad_mono(H_RT2, fit_RT2)*1e3, label="quadratic fit", color='r')
plt.grid()
plt.legend()
plt.title("Force depending on the magnetic field (T=$20.5^\circ$C, heater removed)")
plt.xlabel("Magnetic field strength (mT)")
plt.ylabel("Downforce (mN)")
plt.savefig('RT2.pdf', format='pdf')
plt.show()


#Current at High Temperature (525 Celsius pm 2 Grad)
I_HT = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.51, 4.01, 4.51, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
m0_HT = 63346 #Leergewicht in mg
m_HT = np.array([63346, 63354, 63363, 63377, 63394, 63414, 63438, 63465, 63496, 63530, 63568, 63612, 63653, 63698, 63742, 63788, 63833, 63878, 63922, 63966])
K_HT = (m_HT - m0_HT)*1e-6*g
H_HT = I_to_H(I_HT)
dK_HT = g * np.sqrt((dm*m0_HT)**2 + (dm*m_HT)**2)

fit_HT, covfit_HT = curve_fit(quad_mono, H_HT, K_HT)
dfit_HT = np.sqrt(covfit_HT[0])
print("fit_HT: " + str(fit_HT[0]) + " +/- " + str(dfit_HT[0]))

plt.errorbar(H_HT*1e3, K_HT*1e3, xerr=dH, yerr=dK_HT, fmt='x', capsize=2, label="measurement")
plt.plot(H_HT*1e3, quad_mono(H_HT, fit_HT)*1e3, label="quadratic fit", color='r')
plt.grid()
plt.legend()
plt.title("Force depending on the magnetic field (T=$525^\circ$C)")
plt.xlabel("Magnetic field strength (mT)")
plt.ylabel("Downforce (mN)")
plt.savefig('HT.pdf', format='pdf')
plt.show()


#K as function of Temperature 
#I_T = 6
Temp_T = np.array([627, 586, 549, 507, 470, 431, 387, 350, 315, 269, 229, 190, 152, 112, 71, 40.2, 27.5]) + 273.15
m0_T = np.array([63346, 63346, 63346, 63347, 63347, 63347, 63347, 63347, 63348, 63348, 63348, 63349, 63349, 63349, 63350, 63350, 63350])
m_T = np.array([63581, 63592, 63601, 63615, 63628, 63642, 63661, 63680, 63699, 63729, 63758, 63793, 63839, 63880, 63944, 64005, 64027])
K_T = (m_T-m0_T)*1e-6*g
dK_T = g * np.sqrt((dm*m0_T)**2 + (dm*m_T)**2)
dT = 5

fit_T, covfit_T = curve_fit(curie_fit, Temp_T, K_T)
acm = fit_T[0]
theta = fit_T[1]
dacm = np.sqrt(np.diag(covfit_T)[0])
dtheta = np.sqrt(np.diag(covfit_T)[1])
print("curie_fit:")
print("acm = " + str(acm) + " +/- " + str(dacm))
print("theta = " + str(theta) + " +/- " + str(dtheta))

plt.errorbar(Temp_T, K_T*1e3, xerr=dT, yerr=dK_T, fmt='x', capsize=2, label="measurement")
plt.plot(Temp_T, curie_fit(Temp_T, acm, theta)*1e3, color='r', label="curie law")
plt.grid()
plt.legend()
plt.title("Force depending on the Temperature (I=6A / H=106mT)")
plt.xlabel("Temperature (K)")
plt.ylabel("Downforce (mN)")
plt.savefig('T.pdf', format='pdf')
plt.show()








