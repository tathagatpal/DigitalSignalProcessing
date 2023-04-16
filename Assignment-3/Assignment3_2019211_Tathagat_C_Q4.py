#Q4
import numpy as np
import scipy as sc
from scipy import signal
import math
import cmath
import matplotlib.pyplot as plt

Omega_s= 800
Omega_p = np.array([100, 200])
Omega_stop = np.array([110, 190])

#Calculating desired critical frequencies of the analog filter
omega_p = np.multiply(Omega_p,2*math.pi/Omega_s)
omega_stop = np.multiply(Omega_stop,2*math.pi/Omega_s)
Omega_p_analog  = np.multiply(np.tan(np.multiply(omega_p,0.5)),2)
Omega_stop_analog  = np.multiply(np.tan(np.multiply(omega_stop,0.5)),2)
B = Omega_stop_analog[1]-Omega_stop_analog[0]
Omega_0_analog = math.sqrt(Omega_stop_analog[1]*Omega_stop_analog[0])


#Calculating ripples 
stop_ripple_db = 40
pass_ripple_db = 30
stop_ripple = 10**(-stop_ripple_db/20)
pass_ripple = 10**(-pass_ripple_db/20)

#Calculating pass band and stop band frequencies of the analog LPF
Omega_p_lpf = 1
Omega_stop_lpf = min((Omega_p_analog[1]**2-Omega_0_analog**2)/Omega_p_analog[1]/B, (Omega_0_analog**2-Omega_p_analog[0]**2)/Omega_p_analog[0]/B)

#Chebyshev I filter to obtain the LPF
rp = -20*math.log10(1-pass_ripple)
epsilon = math.sqrt(1/(1-pass_ripple)**2-1)
N = math.ceil(math.log((math.sqrt(1-stop_ripple**2)+math.sqrt(1-stop_ripple**2*(1+epsilon**2)))/stop_ripple/epsilon)/math.log(Omega_stop_lpf/Omega_p_lpf+math.sqrt(Omega_stop_lpf**2/Omega_p_lpf**2-1)))
b_lpf, a_lpf = signal.cheby1(N,rp,Omega_p_lpf,analog=True)

#Plotting the initial low pass filter
plt.subplot(1,2,1)
w_lpf, h_lpf = signal.freqs(b_lpf, a_lpf, worN=1000)
plt.plot(w_lpf,abs(h_lpf))
plt.hlines(10**(-stop_ripple_db/20), 0, 10, 'r', 'dotted')
plt.hlines(1+10**(-pass_ripple_db/20), 0, 10, 'r', 'dotted')
plt.hlines(1-10**(-pass_ripple_db/20), 0, 10, 'r', 'dotted')
plt.xlabel(r'$\Omega$')
plt.ylabel(r'$|H_{LP}(\Omega)|$')
plt.title('Magnitude response (linear scale) of the initial low-pass filter')
plt.grid(linestyle='dashed')

plt.subplot(1,2,2)
plt.plot(w_lpf,20*np.log10(abs(h_lpf)))
plt.hlines(-stop_ripple_db, 0, 10, 'r', 'dotted')
plt.hlines(20*np.log10(1+10**(-pass_ripple_db/20)), 0, 10, 'r', 'dotted')
plt.hlines(20*np.log10(1-10**(-pass_ripple_db/20)), 0, 10, 'r', 'dotted')
plt.xlabel(r'$\Omega$')
plt.ylabel(r'$|H_{LP}(\Omega)|$')
plt.title('Magnitude response (dB) of the initial low-pass filter')
plt.grid(linestyle='dashed')
plt.show()

#Obtaining the analog BSF using frequency transformation
b_bsf, a_bsf = signal.lp2bs(b_lpf, a_lpf, Omega_0_analog, B)

w_bsf, h_bsf = signal.freqs(b_bsf, a_bsf, worN=1000)

plt.subplot(1,2,1)
plt.plot(w_bsf[:600],abs(h_bsf)[:600])
plt.xlabel(r'$\Omega$')
plt.ylabel(r'$|H_{BS}(\Omega)|$')
plt.title('Magnitude response (linear scale) of the initial bandstop filter')
plt.grid(linestyle='dashed')

plt.subplot(1,2,2)
plt.plot(w_bsf[:600],20*np.log10(abs(h_bsf)[:600]))
plt.xlabel(r'$\Omega$')
plt.ylabel(r'$|H_{BS}(\Omega)|$')
plt.title('Magnitude response (dB) of the initial bandstop filter')
plt.grid(linestyle='dashed')

plt.show()

#Applying Bilinear transform
b, a = signal.bilinear(b_bsf, a_bsf)

w, h =signal.freqz(b, a, worN=10000)

w *=Omega_s/2/math.pi

mag_response = 20*np.log10(abs(h))

plt.subplot(1,2,1)
plt.plot(w,abs(h))
plt.axis([0, Omega_s/2, 0, 1.1])
plt.xlabel(r'$\omega$')
plt.ylabel(r'$|H(\omega)|$')
plt.hlines(10**(-stop_ripple_db/20), 0, Omega_s/2, 'r', 'dotted')
plt.hlines(1+10**(-pass_ripple_db/20), 0, Omega_s/2, 'r', 'dotted')
plt.hlines(1-10**(-pass_ripple_db/20), 0, Omega_s/2, 'r', 'dotted')
plt.title('Magnitude response (linear scale)')
plt.grid(linestyle='dashed')

plt.subplot(1,2,2)
plt.plot(w,mag_response)
plt.xlabel(r'$\omega$')
plt.ylabel(r'$|H(\omega)|$ (dB)')
plt.hlines(-stop_ripple_db, 0, Omega_s/2, 'r', 'dotted')
plt.hlines(20*math.log10(1+10**(-pass_ripple_db/20)), 0, Omega_s/2, 'r', 'dotted')
plt.hlines(20*math.log10(1-10**(-pass_ripple_db/20)), 0, Omega_s/2, 'r', 'dotted')
plt.title('Magnitude response (in dB)')
plt.grid(linestyle='dashed')

plt.suptitle('IIR BandStop Filter using Bilinear Transform')
plt.show()