import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
import pucch


# Some system parameters #
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

transition_time = [[1, 25600], [2, 13792]]  # In the multiple of Tc

# print(delta_fmax, Ts, Tc, k, T_f, numerology[0][1], transition_time[0][1]*Tc)

fft_size = 4096       # One symbol
num_slot_sym = 14     # One slot
nFrame = 500      # Frame to simulate

nFrameGrid = []

p4 = pucch.pucch() # Format 3

p4.pucch_format3_param["pucchGroupHopping"] = 'neither'
p4.pucch_format3_param["pucchFrequencyHopping"] = 'disable'
p4.pucch_format3_param["startSymbolIndex"] = 0
p4.pucch_format3_param["nrOfSymbols"] = 5 # 1,2 Only
p4.pucch_format3_param["nPRB"] = 1 # RM code only
p4.pucch_format3_param["startPRB"] = 0 # RM code only
p4.pucch_format3_param["n_rnti"] = 100
p4.pucch_format3_param["cqi_bit_len"] = 6
p4.pucch_format3_param["cqi_bit"] = 0
p4.pucch_format3_param["modBPSK"] = 0

nTxBits = np.random.randint(0, 2, nFrame*p4.pucch_format3_param["cqi_bit_len"])

snr_sweep = []
ber_sweep = []

for snr_db in range(-4, 10, 1):
    nRxBits = []
    print(snr_db)
    for frame in range(0, nFrame, 1):
        nTxGrid = [[complex(0, 0) for x in range(fft_size)] for x in range(num_slot_sym)]

        cqi_bit_start = frame*p4.pucch_format3_param["cqi_bit_len"]
        cqi_bit_stop = (frame+1)*p4.pucch_format3_param["cqi_bit_len"]
        p4.pucch_format3_param["cqi_bit"] = nTxBits[cqi_bit_start:cqi_bit_stop]

        # Create Gird #
        txGrid = p4.pucch_format_3(nTxGrid, 0, 400, 100, 0, p4.pucch_format3_param)
        txVector = [[complex(0, 0) for x in range(fft_size)] for x in range(num_slot_sym)]
        for n in range(0, num_slot_sym, 1):
            txVector[n] = np.fft.ifft(txGrid[n])*np.sqrt(fft_size)

        # Add Cyclic Prefix
        # TBD

        # Add Channel
        snr = 10**(snr_db/10)
        for n in range(0, num_slot_sym, 1):
            sig_power = snr  # Noise Power  == 1
            noise_power = 1
            noise_real = np.sqrt(1/2)*np.random.normal(0, 1, fft_size)
            noise_imag = np.sqrt(1/2)*np.random.normal(0, 1, fft_size)

            txVector[n].real = np.sqrt(sig_power)*txVector[n].real + noise_real
            txVector[n].imag = np.sqrt(sig_power)*txVector[n].imag + noise_imag

        rxVector = [[complex(0, 0) for x in range(fft_size)] for x in range(num_slot_sym)]

        # Remove Cyclic Prefix
        #TBD

        for n in range(0, num_slot_sym, 1):
            rxVector[n] = np.fft.fft(txVector[n])/np.sqrt(fft_size)

        # Receiver
        cqi_bit, snr = p4.pucch_format_3_rec(rxVector, 0, 400, 100, 0, p4.pucch_format3_param, noise_power)

        for n in range(0, p4.pucch_format3_param["cqi_bit_len"]): # Serializing
            nRxBits.append(cqi_bit[n])


    bit_error = np.sum(np.abs(nRxBits-nTxBits))/len(nTxBits)
    snr_sweep.append(snr_db)
    ber_sweep.append(bit_error)
    print(snr_sweep, ber_sweep)

print(snr_sweep, ber_sweep)

X = [x for x in snr_sweep]
Y = [x for x in ber_sweep]

f1 = plt.figure()
plt.semilogy(X, Y, color='red')
plt.grid
plt.show()

