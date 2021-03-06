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
nFrame = 16000      # Frame to simulate

nFrameGrid = []

p2 = pucch.pucch() # Format 1

p2.pucch_format1_param["pucchGroupHopping"] = 'neither'
p2.pucch_format1_param["pucchFrequencyHopping"] = 'disable'
p2.pucch_format1_param["initialCyclicShift"] = 5
p2.pucch_format1_param["nrOfSymbols"] = 6
p2.pucch_format1_param["startSymbolIndex"] = 0
p2.pucch_format1_param["nHarqBit"] = 2   # SR Can be triggered by keeping the HARQ bit as 1 and HARQ bit as 0
p2.pucch_format1_param["harqBit0"] = 0
p2.pucch_format1_param["harqBit1"] = 1
p2.pucch_format1_param["startPRB"] = 10
p2.pucch_format1_param["nHopPRB"] = 20
p2.pucch_format1_param["spread_seq_idx"] = 1

nTxBits = np.random.randint(0, 2, nFrame*p2.pucch_format1_param["nHarqBit"])

snr_sweep = []
ber_sweep = []

for snr_db in range(-11, 2, 1):
    nRxBits = []
    print(snr_db)
    for frame in range(0, nFrame, 1):
        nTxGrid = [[complex(0, 0) for x in range(fft_size)] for x in range(num_slot_sym)]

        if p2.pucch_format1_param["nHarqBit"] == 2:
            p2.pucch_format1_param["harqBit0"] = nTxBits[2 * frame]
            p2.pucch_format1_param["harqBit1"] = nTxBits[2 * frame + 1]
        else:
            p2.pucch_format1_param["harqBit0"] = nTxBits[frame]

        # Create Gird #
        txGrid = p2.pucch_format_1(nTxGrid, 0, 400, 0, p2.pucch_format1_param)

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
        harq_bit = p2.pucch_format_1_rec(rxVector, 0, 400, 0, p2.pucch_format1_param, noise_power)

        harqBits = harq_bit[0].tolist()
        harqBits.reverse()
        nRxBits  = nRxBits + harqBits

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

