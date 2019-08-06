import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
import referenceSignal
import cSequence

# Creating PUCCH Channel Format 1,2,3,4


class pucch(referenceSignal.ReferenceSignal, cSequence.CSequence):
    pucch_resource_set =  [[0, 12, 2, 0, [0, 3]],
                         [0, 12, 2, 0, [0, 4, 8]],
                         [0, 12, 2, 3, [0, 4, 8]],
                         [1, 10, 4, 0, [0, 6]],
                         [1, 10, 4, 0, [0, 3, 6, 9]],
                         [1, 10, 4, 2, [0, 3, 6, 9]],
                         [1, 10, 4, 4, [0, 3, 6, 9]],
                         [1, 4, 10, 0, [0, 6]],
                         [1, 4, 10, 0, [0, 3, 6, 9]],
                         [1, 4, 10, 2, [0, 3, 6, 9]],
                         [1, 4, 10, 4, [0, 3, 6, 9]],
                         [1, 0, 14, 0, [0, 6]],
                         [1, 0, 14, 0, [0, 3, 6, 9]],
                         [1, 0, 14, 2, [0, 3, 6, 9]],
                         [1, 0, 14, 4, [0, 3, 6, 9]],
                         [1, 0, 14, 123, [0, 3, 6, 9]]]
    mcs_one_bit = [0, 6]
    mcs_two_bits = [0, 3, 9, 6]
    N_sc_rb = 12

    pucch_format0_param = {"pucchGroupHopping": 'neither',
                           "initialCyclicShift": -1,
                           "nrOfSymbols": -1,
                           "startSymbolIndex": -1,
                           "nHarqBit": -1,
                           "harqBit0": -1,
                           "harqBit1": -1,
                           "startPRB": -1,
                           "nHopPRB": -1
                           }

    def generate_u_v(self, pucchGroupHopping, pucchFrequencyHopping, n_id, n_sf_u, n_hop):
        # Generating u,v
        if pucchGroupHopping == 'neither':
            fgh = 0
            v = 0
            fss = (n_id % 30)
        elif pucchGroupHopping == 'enable':
            c_init = n_id//30
            c_len = 8*(2*n_sf_u+n_hop)+8+1
            c_seq = self.generate_c_sequence(c_init, c_len)
            c_seq = c_seq[-8:]
            fgh = cinit_calc(c_seq, 8)
            v = 0
            fss = (n_id % 30)
        else:
            c_init = (2**5)*(n_id // 30) + n_id % 30
            c_len = (2 * n_sf_u + n_hop + 1)
            c_seq = self.generate_c_sequence(c_init, c_len)
            fgh = 0
            fss = (n_id % 30)
            v = c_seq[-1]

        u = (fgh+fss) % 30
        return [u, v]

    def cyclic_shift_ncs(self, n_sf_u, l, l_dash, m_o, m_cs, n_id):

        c_init = n_id
        c_len = 14*8*n_sf_u+8*(l+l_dash)+8+1
        c_seq = self.generate_c_sequence(c_init, c_len)
        c_seq = c_seq[-8:]
        ncs_n_sf_u_l = self.cinit_calc(c_seq, 8)
        cyclic_shift = 2*cmath.pi/12*((m_o + m_cs + ncs_n_sf_u_l) % 12)
        return cyclic_shift

    def pucch_format_0(self, nGrid, n_sf_u, n_id, n_hop, pucch_format0_param):

        r_u_v_c_alpha_sym_0 = []
        r_u_v_c_alpha_sym_1 = []
        l = pucch_format0_param["startSymbolIndex"]

        # Generate Cyclic shift for symbol 0
        m_o = pucch_format0_param["initialCyclicShift"]
        m_cs = 0

        if pucch_format0_param["nHarqBit"] == 1:

            m_cs = self.mcs_one_bit[pucch_format0_param["harqBit0"]]

        elif pucch_format0_param["nHarqBit"] == 2:

            index = pucch_format0_param["harqBit0"]+2*pucch_format0_param["harqBit1"]
            m_cs = self.mcs_two_bits[index]

        else:
            print("PUCCH / pucch_format_0: Error with no of Harq bits")
        if pucch_format0_param["nrOfSymbols"] == 1:

            [u, v] = self.generate_u_v(pucch_format0_param["pucchGroupHopping"], pucch_format0_param["pucchFrequencyHopping"], n_id, n_sf_u, n_hop)
            cyclic_shift_sym_0 = self.cyclic_shift_ncs(n_sf_u, l, 0, m_o, m_cs, n_id)
            r_u_v_c_alpha_sym_0 = self.reference_seq(u, v, cyclic_shift_sym_0, 1, 0)

        elif pucch_format0_param["nrOfSymbols"] == 2:

            [u, v] = self.generate_u_v(pucch_format0_param["pucchGroupHopping"], pucch_format0_param["pucchFrequencyHopping"], n_id, n_sf_u, n_hop)
            cyclic_shift_sym_0 = self.cyclic_shift_ncs(n_sf_u, l, 0, m_o, m_cs, n_id)
            cyclic_shift_sym_1 = self.cyclic_shift_ncs(n_sf_u, l+1, 0, m_o, m_cs, n_id)
            r_u_v_c_alpha_sym_0 = self.reference_seq(u, v, cyclic_shift_sym_0, 1, 0)
            r_u_v_c_alpha_sym_1 = self.reference_seq(u, v, cyclic_shift_sym_1, 1, 0)

        else:
            print("PUCCH / pucch_format_0: Error with no of pucch symbols")

        # Mapping to PRB
        startSymbolIndex = pucch_format0_param["startSymbolIndex"]
        startPRB = pucch_format0_param["startPRB"]

        if pucch_format0_param["nrOfSymbols"] == 1:

            for n in range(0, 12, 1):
               nGrid[startSymbolIndex][startPRB*12+n] = r_u_v_c_alpha_sym_0[n]

        elif pucch_format0_param["nrOfSymbols"] == 2:

            for n in range(0, 12, 1):

                nGrid[startSymbolIndex][startPRB*12 + n] = r_u_v_c_alpha_sym_0[n]
                nGrid[startSymbolIndex+1][startPRB*12 + n] = r_u_v_c_alpha_sym_1[n]

        return nGrid

    def pucch_format_0_rec(self, nGrid, n_sf_u, n_id, n_hop, pucch_format0_param, noise_power):

        r_u_v_c_alpha_sym_0 = []
        r_u_v_c_alpha_sym_1 = []

        l = pucch_format0_param["startSymbolIndex"]

        # Generate Cyclic shift for symbol 0

        m_o = pucch_format0_param["initialCyclicShift"]

        if pucch_format0_param["nrOfSymbols"] == 1:

            [u, v] = self.generate_u_v(pucch_format0_param["pucchGroupHopping"], pucch_format0_param["pucchFrequencyHopping"], n_id, n_sf_u, n_hop)
            cyclic_shift_sym_0 = self.cyclic_shift_ncs(n_sf_u, l, 0, m_o, 0, n_id)
            r_u_v_c_alpha_sym_0 = self.reference_seq(u, v, cyclic_shift_sym_0, 1, 0)

        elif pucch_format0_param["nrOfSymbols"] == 2:

            [u, v] = self.generate_u_v(pucch_format0_param["pucchGroupHopping"], pucch_format0_param["pucchFrequencyHopping"], n_id, n_sf_u, n_hop)
            cyclic_shift_sym_0 = self.cyclic_shift_ncs(n_sf_u, l, 0, m_o, 0, n_id)
            cyclic_shift_sym_1 = self.cyclic_shift_ncs(n_sf_u, l+1, 0, m_o, 0, n_id)
            r_u_v_c_alpha_sym_0 = self.reference_seq(u, v, cyclic_shift_sym_0, 1, 0)
            r_u_v_c_alpha_sym_1 = self.reference_seq(u, v, cyclic_shift_sym_1, 1, 0)

        else:
            print("PUCCH / pucch_format_0: Error with no of pucch symbols")

        # Mapping to PRB
        startSymbolIndex = pucch_format0_param["startSymbolIndex"]
        startPRB = pucch_format0_param["startPRB"]

        rec_symbol = []
        if pucch_format0_param["nrOfSymbols"] == 1:

            # Get PUCCH format 0 Symbol and derotate with reference symbol
            for n in range(0, 12, 1):
                rec_symbol.append(nGrid[startSymbolIndex][startPRB*12 + n]*np.conj(r_u_v_c_alpha_sym_0[n]))
            noisePower = noise_power[startSymbolIndex]

        elif pucch_format0_param["nrOfSymbols"] == 2:  # Coherent Combining

            for n in range(0, 12, 1):
                rec_symbol.append(nGrid[startSymbolIndex][startPRB * 12 + n] * np.conj(r_u_v_c_alpha_sym_0[n]))
                rec_symbol[n] += nGrid[startSymbolIndex+1][startPRB * 12 + n] * np.conj(r_u_v_c_alpha_sym_1[n])

            noisePower = (noise_power[startSymbolIndex] + noise_power[startSymbolIndex+1])/2


        # if Harq Bit is 1
        if pucch_format0_param["nHarqBit"] == 1:

            corr_bit_1 = 0

            corr_bit_0 = np.sum(rec_symbol)

            for n in range(0, 12, 1):
                corr_bit_1 += rec_symbol[n]*cmath.exp(complex(0, -math.pi*n))

            if corr_bit_1.real > corr_bit_0.real:
                harq_bit = 1
                sig_power = np.absolute(corr_bit_1)**2
                if sig_power > noisePower:
                    dtx = 0
                else:
                    dtx = 1
            else:
                harq_bit = 0
                sig_power = np.absolute(corr_bit_0)**2
                if sig_power > noisePower:
                    dtx = 0
                else:
                    dtx = 1
        elif pucch_format0_param["nHarqBit"] == 2:

            corr_bit_0_1 = 0
            corr_bit_1_0 = 0
            corr_bit_1_1 = 0

            corr_bit_0_0 = np.sum(rec_symbol)

            for n in range(0, 12, 1):
                corr_bit_0_1 += rec_symbol[n]*cmath.exp(complex(0, -math.pi/2 * n))
                corr_bit_1_0 += rec_symbol[n]*cmath.exp(complex(0, -3*math.pi/2 * n))
                corr_bit_1_1 += rec_symbol[n]*cmath.exp(complex(0, -math.pi * n))

            dec_mat = [corr_bit_0_0.real, corr_bit_0_1.real, corr_bit_1_0.real, corr_bit_1_1.real]

            harq_bit = dec_mat.index(max(dec_mat))

            if harq_bit == 0:
                harq_bit = [0, 0]
                sig_power = np.absolute(corr_bit_0_0 / 12) ** 2
                if sig_power > noisePower:
                    dtx = 0
                else:
                    dtx = 1

            if harq_bit == 1:
                harq_bit = [0, 1]
                sig_power = np.absolute(corr_bit_0_1 / 12) ** 2
                if sig_power > noisePower:
                    dtx = 0
                else:
                    dtx = 1

            if harq_bit == 2:
                harq_bit = [1, 0]
                sig_power = np.absolute(corr_bit_1_0 / 12) ** 2
                if sig_power > noisePower:
                    dtx = 0
                else:
                    dtx = 1

            if harq_bit == 3:
                harq_bit = [1, 1]
                sig_power = np.absolute(corr_bit_1_1 / 12) ** 2
                if sig_power > noisePower:
                    dtx = 0
                else:
                    dtx = 1
        return [harq_bit, dtx]


