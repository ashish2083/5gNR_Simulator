import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
import referenceSignal
import cSequence

# Creating PUCCH Channel Format 1,2,3,4


class pucch(referenceSignal.ReferenceSignal, cSequence.CSequence):

    pucch_resource_set = [[0, 12, 2, 0, [0, 3]],
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

    pucch_format_1_table_data = [[4, 2, 1, 1],
                                 [5, 2, 1, 1],
                                 [6, 3, 1, 2],
                                 [7, 3, 1, 2],
                                 [8, 4, 2, 2],
                                 [9, 4, 2, 2],
                                 [10, 5, 2, 3],
                                 [11, 5, 2, 3],
                                 [12, 6, 3, 3],
                                 [13, 6, 3, 3],
                                 [14, 7, 3, 4]]

    pucch_format_1_table_dmrs = [[4, 2, 1, 1],
                                 [5, 3, 1, 2],
                                 [6, 3, 2, 1],
                                 [7, 4, 2, 2],
                                 [8, 4, 2, 2],
                                 [9, 5, 2, 3],
                                 [10, 5, 3, 2],
                                 [11, 6, 3, 3],
                                 [12, 6, 3, 3],
                                 [13, 7, 3, 4],
                                 [14, 7, 4, 3]]

    pucch_format_1_i_0_table = [[0], [0, 0], [0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0]]

    pucch_format_1_i_1_table = [[0], [0, 1], [0, 1, 2], [0, 2, 0, 2], [0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5],
                                [0, 1, 2, 3, 4, 5, 6]]

    pucch_format_1_i_2_table = [[0], [0, 0], [0, 2, 1], [0, 0, 2, 2], [0, 2, 4, 1, 3], [0, 2, 4, 0, 2, 4],
                                [0, 2, 4, 6, 1, 3, 5]]

    pucch_format_1_i_3_table = [[0], [0, 0], [0, 0, 0], [0, 2, 2, 0], [0, 3, 1, 4, 2], [0, 3, 0, 3, 0, 3],
                                [0, 3, 6, 2, 5, 1, 4]]

    pucch_format_1_i_4_table = [[0], [0, 0], [0, 0, 0], [0, 0, 0, 0], [0, 4, 3, 2, 1], [0, 4, 2, 0, 4, 2],
                                [0, 4, 1, 5, 2, 6, 3]]

    pucch_format_1_i_5_table = [[0], [0, 0], [0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 5, 4, 3, 2, 1],
                                [0, 5, 3, 1, 6, 4, 2]]

    pucch_format_1_i_6_table = [[0], [0, 0], [0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0],
                                [0, 6, 5, 4, 3, 2, 1]]
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
    pucch_format1_param = {"pucchGroupHopping": 'neither',
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

    def pucch_format_1(self, nGrid, n_sf_u, n_id, n_hop, pucch_format1_param):

        # Modulating bits
        mod_sym = complex(0,0)
        if pucch_format1_param["nHarqBit"] == 1:
            mod_sym = 1 / np.sqrt(2) * complex(1 - 2 * pucch_format1_param["harqBit0"],
                                               1 - 2 * pucch_format1_param["harqBit0"])  # BPSK
        else:
            mod_sym = 1 / np.sqrt(2) * complex(1 - 2 * pucch_format1_param["harqBit0"],
                                               1 - 2 * pucch_format1_param["harqBit1"])  # QPSK

        # Generate DMRS Sequence

        l = pucch_format1_param["startSymbolIndex"]
        n_symbol = pucch_format1_param["nrOfSymbols"]
        n_pucch_1_sf_m_dash = self.pucch_format_1_table_dmrs[n_symbol-4][1]  #No support for pucch hopping for now

        n_dmrs_symbol = (n_symbol + 1) // 2

        i = pucch_format1_param["spread_seq_idx"]

        if i == 0:
            w_i = self.pucch_format_1_i_0_table[n_pucch_1_sf_m_dash-1]

        if i == 1:
            w_i = self.pucch_format_1_i_1_table[n_pucch_1_sf_m_dash - 1]

        if i == 2:
            w_i = self.pucch_format_1_i_2_table[n_pucch_1_sf_m_dash - 1]

        if i == 3:
            w_i = self.pucch_format_1_i_3_table[n_pucch_1_sf_m_dash - 1]

        if i == 4:
            w_i = self.pucch_format_1_i_4_table[n_pucch_1_sf_m_dash - 1]

        if i == 5:
            w_i = self.pucch_format_1_i_5_table[n_pucch_1_sf_m_dash - 1]

        if i == 6:
            w_i = self.pucch_format_1_i_6_table[n_pucch_1_sf_m_dash - 1]

        if i == 7:
            w_i = self.pucch_format_1_i_7_table[n_pucch_1_sf_m_dash - 1]

        [u, v] = self.generate_u_v(pucch_format1_param["pucchGroupHopping"],
                                   pucch_format1_param["pucchFrequencyHopping"], n_id, n_sf_u, n_hop)

        dmrs_sequence = []
        for n in range(0, n_pucch_1_sf_m_dash, 1):
            m_o = pucch_format1_param["initialCyclicShift"]
            cyclic_shift_sym = self.cyclic_shift_ncs(n_sf_u, (l+n*2), 0, m_o, 0, n_id)
            r_u_v_c_alpha_seq = self.reference_seq(u, v, cyclic_shift_sym, 1, 0)
            seq_mult = cmath.exp(complex(0, 2*math.pi * w_i[n]/n_dmrs_symbol))
            dmrs_sequence = np.r_[dmrs_sequence, np.dot(seq_mult, r_u_v_c_alpha_seq)]

        # Generate Data Sequence

        l = pucch_format1_param["startSymbolIndex"]
        n_symbol = pucch_format1_param["nrOfSymbols"]
        n_pucch_1_sf_m_dash = self.pucch_format_1_table_data[n_symbol - 4][1]  # No support for pucch hopping for now

        n_data_symbol = n_symbol - n_dmrs_symbol

        if i == 0:
            w_i = self.pucch_format_1_i_0_table[n_pucch_1_sf_m_dash - 1]

        if i == 1:
            w_i = self.pucch_format_1_i_1_table[n_pucch_1_sf_m_dash - 1]

        if i == 2:
            w_i = self.pucch_format_1_i_2_table[n_pucch_1_sf_m_dash - 1]

        if i == 3:
            w_i = self.pucch_format_1_i_3_table[n_pucch_1_sf_m_dash - 1]

        if i == 4:
            w_i = self.pucch_format_1_i_4_table[n_pucch_1_sf_m_dash - 1]

        if i == 5:
            w_i = self.pucch_format_1_i_5_table[n_pucch_1_sf_m_dash - 1]

        if i == 6:
            w_i = self.pucch_format_1_i_6_table[n_pucch_1_sf_m_dash - 1]

        if i == 7:
            w_i = self.pucch_format_1_i_7_table[n_pucch_1_sf_m_dash - 1]

        [u, v] = self.generate_u_v(pucch_format1_param["pucchGroupHopping"],
                                   pucch_format1_param["pucchFrequencyHopping"], n_id, n_sf_u, n_hop)

        data_sequence = []

        for n in range(0, n_pucch_1_sf_m_dash, 1):
            m_o = pucch_format1_param["initialCyclicShift"]
            cyclic_shift_sym = self.cyclic_shift_ncs(n_sf_u, (l + n * 2 + 1), 0, m_o, 0, n_id)
            r_u_v_c_alpha_seq = self.reference_seq(u, v, cyclic_shift_sym, 1, 0)
            spread_coefficient = mod_sym*cmath.exp(complex(0, 2 * math.pi * w_i[n] / n_data_symbol))
            data_sequence = np.r_[data_sequence, np.dot(spread_coefficient, r_u_v_c_alpha_seq)]

        # Mapping on to the Grid
        startPRB = pucch_format1_param["startPRB"]

        # DMRS
        for n in range(0, n_dmrs_symbol, 1):
            for m in range(0, self.N_sc_rb):
                nGrid[l+2*n][startPRB*12 + m] = dmrs_sequence[n*self.N_sc_rb +m]

        # Data

        for n in range(0, n_data_symbol, 1):
            for m in range(0, self.N_sc_rb):
                nGrid[l + (2*n)+1][startPRB * 12 + m] = data_sequence[n * self.N_sc_rb + m]

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

        elif pucch_format0_param["nrOfSymbols"] == 2:  # Coherent Combining

            for n in range(0, 12, 1):
                rec_symbol.append(nGrid[startSymbolIndex][startPRB * 12 + n] * np.conj(r_u_v_c_alpha_sym_0[n])/2)
                rec_symbol[n] += nGrid[startSymbolIndex+1][startPRB * 12 + n] * np.conj(r_u_v_c_alpha_sym_1[n])/2

        # if Harq Bit is 1
        if pucch_format0_param["nHarqBit"] == 1:

            corr_bit_1 = 0

            corr_bit_0 = np.sum(rec_symbol)

            for n in range(0, 12, 1):
                corr_bit_1 += rec_symbol[n]*cmath.exp(complex(0, -math.pi*n))

            if corr_bit_1.real > corr_bit_0.real:
                harq_bit = np.array(1)
                sig_power = np.absolute(corr_bit_1/12)**2

                if sig_power > noise_power:
                    dtx = 0
                else:
                    dtx = 1
            else:
                harq_bit = np.array(0)
                sig_power = np.absolute(corr_bit_0/12)**2

                if sig_power > noise_power:
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
                harq_bit = np.array([0, 0])
                sig_power = np.absolute(corr_bit_0_0/12) ** 2

                if sig_power > noise_power:
                    dtx = 0
                else:
                    dtx = 1

            if harq_bit == 1:
                harq_bit = np.array([0, 1])
                sig_power = np.absolute(corr_bit_0_1/12) ** 2
                if sig_power > noise_power:
                    dtx = 0
                else:
                    dtx = 1

            if harq_bit == 2:
                harq_bit = np.array([1, 0])
                sig_power = np.absolute(corr_bit_1_0/12) ** 2
                if sig_power > noise_power:
                    dtx = 0
                else:
                    dtx = 1

            if harq_bit == 3:
                harq_bit = np.array([1, 1])
                sig_power = np.absolute(corr_bit_1_1/12) ** 2
                if sig_power > noise_power:
                    dtx = 0
                else:
                    dtx = 1
        return [harq_bit, dtx]

    def pucch_format_1_rec(self, nGrid, n_sf_u, n_id, n_hop, pucch_format1_param,noise_power):

        # Generate DMRS Sequence

        l = pucch_format1_param["startSymbolIndex"]
        n_symbol = pucch_format1_param["nrOfSymbols"]
        n_pucch_1_sf_m_dash = self.pucch_format_1_table_dmrs[n_symbol - 4][1]  # No support for pucch hopping for now

        n_dmrs_symbol = (n_symbol + 1) // 2

        i = pucch_format1_param["spread_seq_idx"]

        if i == 0:
            w_i = self.pucch_format_1_i_0_table[n_pucch_1_sf_m_dash - 1]

        if i == 1:
            w_i = self.pucch_format_1_i_1_table[n_pucch_1_sf_m_dash - 1]

        if i == 2:
            w_i = self.pucch_format_1_i_2_table[n_pucch_1_sf_m_dash - 1]

        if i == 3:
            w_i = self.pucch_format_1_i_3_table[n_pucch_1_sf_m_dash - 1]

        if i == 4:
            w_i = self.pucch_format_1_i_4_table[n_pucch_1_sf_m_dash - 1]

        if i == 5:
            w_i = self.pucch_format_1_i_5_table[n_pucch_1_sf_m_dash - 1]

        if i == 6:
            w_i = self.pucch_format_1_i_6_table[n_pucch_1_sf_m_dash - 1]

        if i == 7:
            w_i = self.pucch_format_1_i_7_table[n_pucch_1_sf_m_dash - 1]

        [u, v] = self.generate_u_v(pucch_format1_param["pucchGroupHopping"],
                                   pucch_format1_param["pucchFrequencyHopping"], n_id, n_sf_u, n_hop)

        dmrs_sequence = []
        for n in range(0, n_pucch_1_sf_m_dash, 1):
            m_o = pucch_format1_param["initialCyclicShift"]
            cyclic_shift_sym = self.cyclic_shift_ncs(n_sf_u, (l + n * 2), 0, m_o, 0, n_id)
            r_u_v_c_alpha_seq = self.reference_seq(u, v, cyclic_shift_sym, 1, 0)
            seq_mult = cmath.exp(complex(0, 2 * math.pi * w_i[n] / n_dmrs_symbol))
            dmrs_sequence = np.r_[dmrs_sequence, np.dot(seq_mult, r_u_v_c_alpha_seq)]

        # Generate Data Sequence

        l = pucch_format1_param["startSymbolIndex"]
        n_symbol = pucch_format1_param["nrOfSymbols"]
        n_pucch_1_sf_m_dash = self.pucch_format_1_table_data[n_symbol - 4][1]  # No support for pucch hopping for now

        n_data_symbol = n_symbol - n_dmrs_symbol

        if i == 0:
            w_i = self.pucch_format_1_i_0_table[n_pucch_1_sf_m_dash - 1]

        if i == 1:
            w_i = self.pucch_format_1_i_1_table[n_pucch_1_sf_m_dash - 1]

        if i == 2:
            w_i = self.pucch_format_1_i_2_table[n_pucch_1_sf_m_dash - 1]

        if i == 3:
            w_i = self.pucch_format_1_i_3_table[n_pucch_1_sf_m_dash - 1]

        if i == 4:
            w_i = self.pucch_format_1_i_4_table[n_pucch_1_sf_m_dash - 1]

        if i == 5:
            w_i = self.pucch_format_1_i_5_table[n_pucch_1_sf_m_dash - 1]

        if i == 6:
            w_i = self.pucch_format_1_i_6_table[n_pucch_1_sf_m_dash - 1]

        if i == 7:
            w_i = self.pucch_format_1_i_7_table[n_pucch_1_sf_m_dash - 1]

        [u, v] = self.generate_u_v(pucch_format1_param["pucchGroupHopping"],
                                   pucch_format1_param["pucchFrequencyHopping"], n_id, n_sf_u, n_hop)

        data_sequence = []
        for n in range(0, n_pucch_1_sf_m_dash, 1):
            m_o = pucch_format1_param["initialCyclicShift"]
            cyclic_shift_sym = self.cyclic_shift_ncs(n_sf_u, (l + n * 2 + 1), 0, m_o, 0, n_id)
            r_u_v_c_alpha_seq = self.reference_seq(u, v, cyclic_shift_sym, 1, 0)
            spread_coefficient = cmath.exp(complex(0, 2 * math.pi * w_i[n] / n_data_symbol))
            data_sequence = np.r_[data_sequence, np.dot(spread_coefficient, r_u_v_c_alpha_seq)]

        # Getting Start PRB
        startPRB = pucch_format1_param["startPRB"]

        # Channel Estimate - Using simple averaging in as the PRB are less and channel is assumed to coherent during the
        # duration of the PUCCH transmission. Can be changed based on use cases.
        channel_estimate = 0
        for n in range(0, n_dmrs_symbol, 1):
            for m in range(0, self.N_sc_rb):
                channel_estimate += nGrid[l + 2 * n][startPRB * 12 + m]*np.conj(dmrs_sequence[n * self.N_sc_rb + m])

        channel_estimate = channel_estimate/n_dmrs_symbol/self.N_sc_rb

        # Equalizing Data / MMSE receiever
        data_estimate = 0
        for n in range(0, n_data_symbol, 1):
            for m in range(0, self.N_sc_rb):
                data_estimate += nGrid[l + 2 * n + 1][startPRB * 12 + m]*np.conj(data_sequence[n * self.N_sc_rb + m])

        data_estimate = data_estimate/n_data_symbol/self.N_sc_rb

        equalized_data = data_estimate*np.conj(channel_estimate)/(np.absolute(channel_estimate)**2 + noise_power) # MMSE Equalizer

        harq_bit = np.array([0])

        sig_power = np.absolute(data_estimate)**2

        # Demodulation
        if pucch_format1_param["nHarqBit"] == 1:
           if equalized_data.real < 0 and equalized_data.imag < 0:
                harq_bit = np.array([1])
        else:
            if equalized_data.real > 0:
                if equalized_data.imag > 0:
                    harq_bit = np.array([0, 0])
                else:
                    harq_bit = np.array([1, 0])
            else:
                if equalized_data.imag < 0:
                    harq_bit = np.array([1, 1])
                else:
                    harq_bit = np.array([0, 1])

        dtx = 0

        if sig_power < noise_power:
            dtx = 1

       # print(sig_power, harq_bit, dtx)

        return [harq_bit, dtx]



