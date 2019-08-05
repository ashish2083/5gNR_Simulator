import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
import pucch

delta_fmax = 480e+3
N_f = 4096
delta_fref = 15e+3
N_fref = 2048

Ts = 1/(delta_fref*N_fref)
Tc = 1/(delta_fmax*N_f)
T_f = ((delta_fmax*N_f)/100)*Tc
k = Ts/Tc

numerology = [[0, 15e+3, 14, 10, 1], [1, 30e+3, 14, 20, 2], [2, 60e+3, 14, 40, 4], [3, 120e+3, 14, 80, 8],
              [4, 240e+3, 160, 16]]

transition_time = [[1, 25600], [2, 13792]] # In the multiple of Tc

# print(delta_fmax, Ts, Tc, k, T_f, numerology[0][1], transition_time[0][1]*Tc)

fft_size = 4096       #one symbol
num_slot_sym = 14     #one_slot

nTxGrid = [[complex(0, 0) for x in range(fft_size)] for x in range(num_slot_sym)]

p1 = pucch.pucch()

p1.pucch_format0_param["pucchGroupHopping"] = 'neither'
p1.pucch_format0_param["pucchFrequencyHopping"] = 'disable'
p1.pucch_format0_param["initialCyclicShift"] = 5
p1.pucch_format0_param["nrOfSymbols"] = 2
p1.pucch_format0_param["startSymbolIndex"] = 10
p1.pucch_format0_param["nHarqBit"] = 2
p1.pucch_format0_param["harqBit0"] = 0
p1.pucch_format0_param["harqBit1"] = 1
p1.pucch_format0_param["startPRB"] = 10
p1.pucch_format0_param["nHopPRB"] = 20


# Create Gird #

txGrid = p1.pucch_format_0(nTxGrid, 0, 400, 0, p1.pucch_format0_param)
#txGrid = p1.pucch_format_0(nTxGrid, 0, 400, 0, p1.pucch_format0_param)

txVector = [[complex(0, 0) for x in range(fft_size)] for x in range(num_slot_sym)]

for n in range(0, num_slot_sym, 1):
    txVector[n] = np.fft.ifft(txGrid[n])*np.sqrt(fft_size)

snr_db = 10 #dB
snr = np.sqrt(10**(-snr_db/10))
for n in range(0, num_slot_sym, 1):
    noise_real = snr*np.random.normal(0, 1, 4096)
    noise_imag = snr*np.random.normal(0, 1, 4096)
    txVector[n].real = txVector[n].real + noise_real
    txVector[n].imag = txVector[n].imag + noise_imag


rxVector = [[complex(0, 0) for x in range(fft_size)] for x in range(num_slot_sym)]

for n in range(0, num_slot_sym, 1):
    rxVector[n] = np.fft.fft(txVector[n])/np.sqrt(fft_size)

# Receiver
harq_bit = p1.pucch_format_0_rec(rxVector, 0, 400, 0, p1.pucch_format0_param)

print(harq_bit)





# X = [x.real for x in p2[10]]
# Y = [x.imag for x in p2[10]]
# f1 = plt.figure()
#plt.scatter(X, Y, color='red')

# for n in range(0, 14, 1):
#    for m in range(0, 4096,1):
#        X[n][m] = abs(p2[n][m])

#f2 = plt.figure()
plt.plot(abs(txVector[10]), color='blue')

# plt.show()

# t = np.arange(400)
# n = np.zeros((400,), dtype=complex)
# n[40:60] = np.exp(1j*np.random.uniform(0, 2*np.pi, (20,)))
# s = np.fft.ifft(n)
# plt.plot(t, s.real, 'b-', t, s.imag, 'r--')
# plt.legend(('real', 'imaginary'))
# plt.show()
