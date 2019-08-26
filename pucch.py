import cmath
import math
import matplotlib.pyplot as plt
import numpy as np
import referenceSignal
import cSequence
import encoder
# Creating PUCCH Channel Format 1,2,3,4


class pucch(referenceSignal.ReferenceSignal, cSequence.CSequence, encoder.encoder):

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
    mcs_one_bit_sr = [3, 9]
    mcs_two_bits = [0, 3, 9, 6]
    mcs_two_bits_sr = [1, 4, 7, 10]
    N_sc_rb = 12

    pucch_format0_param = {"pucchGroupHopping": 'neither',
                           "initialCyclicShift": -1,
                           "nrOfSymbols": -1,
                           "startSymbolIndex": -1,
                           "nHarqBit": -1,
                           "harqBit0": -1,
                           "harqBit1": -1,
                           "sr": -1,
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
    pucch_format2_param = {"startSymbolIndex": -1,
                           "nrOfSymbols": -1,
                           "nPRB": -1,
                           "startPRB": -1,
                           "n_rnti": -1,
                           "cqi_bit_len": -1,
                           "cqi_bit": -1
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
            if pucch_format0_param["sr"] == 1:
                m_cs = self.mcs_one_bit_sr[pucch_format0_param["harqBit0"]]
            else:
                m_cs = self.mcs_one_bit[pucch_format0_param["harqBit0"]]

        elif pucch_format0_param["nHarqBit"] == 2:

            index = pucch_format0_param["harqBit0"]+2*pucch_format0_param["harqBit1"]
            if pucch_format0_param["sr"] == 1:
                m_cs = self.mcs_two_bits_sr[index]
            else:
                m_cs = self.mcs_two_bits[index]

        else:
            if pucch_format0_param["sr"] == 1:
                if pucch_format0_param["sr"] == 1:
                    m_cs = 0
                else:
                    print("PUCCH / pucch_format_0: Invalid format")

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

    def pucch_format_2(self, nGrid, n_sf_u, n_id, n_id_0, n_hop, pucch_format2_param):
        # Generate Scrambling cSequence
        if pucch_format2_param["cqi_bit_len"] < 3:
            print("PUCCH:pucch_format_2 Invalid no of CQI bits")

        if pucch_format2_param["cqi_bit_len"] <= 11:
            coded_bit = self.reed_muller_encoder(pucch_format2_param["cqi_bit"], pucch_format2_param["cqi_bit_len"])
            c_len = 32

        # TBD Polar Code here

        cinit = pucch_format2_param["n_rnti"]*(2**15) + n_id
        c_seq = self.generate_c_sequence(cinit, c_len)

        # Scramble
        scram_bit = (coded_bit + c_seq) % 2

        data_sym_array = []

        # Modulation QPSK
        for n in range(0, c_len, 2):
            mod_sym = 1 / np.sqrt(2) * complex(1 - 2 * scram_bit[n],
                                               1 - 2 * scram_bit[n+1])  # QPSK
            data_sym_array.append(mod_sym)

        # Generating Reference Signal
        l = pucch_format2_param["startSymbolIndex"]
        n_symbol = pucch_format2_param["nrOfSymbols"]
        n_resource_block = pucch_format2_param["nPRB"]
        n_resource_element = n_resource_block*self.N_sc_rb
        n_pilot_bits = int(n_resource_block*4)*2 # QPSK Bits
        n_pilot_sym  = (n_resource_block*4)*n_symbol
        n_data_sym   = (n_resource_block*8)*n_symbol
        ref_sym_array = []

        if n_symbol == 1:
            c_init = ((2**17)*(14*n_sf_u + l + 1)*(2*n_id_0+1) + 2*n_id_0) % (2**31)
            c_len = n_pilot_bits
            c_seq = self.generate_c_sequence(c_init, c_len)

            for n in range(0, n_pilot_bits, 2):
                mod_sym = 1 / np.sqrt(2) * complex(1 - 2 * c_seq[n],
                                                   1 - 2 * c_seq[n+1])  # QPSK
                ref_sym_array.append(mod_sym)
        else: # 2 Symbol

            # First Symbol
            c_init = ((2**17)*(14*n_sf_u + l + 1)*(2*n_id_0+1) + 2*n_id_0) % (2**31)
            c_len = n_pilot_bits
            c_seq = self.generate_c_sequence(c_init, c_len)

            for n in range(0, n_pilot_bits, 2):
                mod_sym = 1 / np.sqrt(2) * complex(1 - 2 * c_seq[n],
                                                   1 - 2 * c_seq[n+1])  # QPSK
                ref_sym_array.append(mod_sym)
            # Second Symbol
            c_init = ((2**17)*(14*n_sf_u + l+1 + 1)*(2*n_id_0+1) + 2*n_id_0) % (2**31)
            c_len = n_pilot_bits
            c_seq = self.generate_c_sequence(c_init, c_len)
            for n in range(0, n_pilot_bits, 2):
                mod_sym = 1 / np.sqrt(2) * complex(1 - 2 * c_seq[n],
                                                   1 - 2 * c_seq[n+1])  # QPSK
                ref_sym_array.append(mod_sym)

        # Mapping on to the resource Grid
        startPRB = pucch_format2_param["startPRB"]

        # DMRS
        if n_symbol == 1:
            dmrs_loc = 0
            data_loc = 0

            for m in range(0, n_resource_block*self.N_sc_rb,1):
                if(m % (3*dmrs_loc+1)) == 0 and m !=0:
                    nGrid[l][startPRB*12 + m] = ref_sym_array[dmrs_loc]
                    dmrs_loc += 1
                else:
                    nGrid[l][startPRB*12 + m] = data_sym_array[data_loc]
                    data_loc += 1


        else:
            dmrs_loc = 0
            data_loc = 0

            for m in range(0, n_resource_block*self.N_sc_rb, 1):
                if(m % (3*dmrs_loc+1)) == 0 and m != 0:
                    nGrid[l][startPRB*12 + m] = 0 #ref_sym_array[dmrs_loc]
                    nGrid[l+1][startPRB*12 + m] = 0 #*ref_sym_array[dmrs_loc + int(n_pilot_sym/2)]
                    dmrs_loc += 1
                else:
                    nGrid[l][startPRB*12 + m] = data_sym_array[data_loc]
                    nGrid[l+1][startPRB*12 + m] = data_sym_array[data_loc +int(n_data_sym/2)]
                    data_loc += 1
            return nGrid

    def pucch_format_0_rec(self, nGrid, n_sf_u, n_id, n_hop, pucch_format0_param, noise_power):

        expn = np.vectorize(cmath.exp)
        cmplx = np.vectorize(complex)

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
            rec_symbol = nGrid[startSymbolIndex][startPRB*12 + np.r_[:12]]*np.conj(r_u_v_c_alpha_sym_0)

        elif pucch_format0_param["nrOfSymbols"] == 2:  # Coherent Combining

            rec_symbol = nGrid[startSymbolIndex][startPRB*12 + np.r_[:12]]*np.conj(r_u_v_c_alpha_sym_0)
            rec_symbol += nGrid[startSymbolIndex+1][startPRB*12 + np.r_[:12]]*np.conj(r_u_v_c_alpha_sym_1)

        # if Harq Bit is 1
        sr = 0
        harq_bit = np.array([0])
        if pucch_format0_param["nHarqBit"] == 1:
            corr_bit_0 = np.sum(rec_symbol)
            corr_bit_1 = np.sum(rec_symbol*expn(cmplx(0, -math.pi*np.r_[:12])))
            corr_bit_0_sr = np.sum(rec_symbol*expn(cmplx(0, -math.pi/2 *np.r_[:12])))
            corr_bit_1_sr = np.sum(rec_symbol*expn(cmplx(0, -3*math.pi/2 *np.r_[:12])))

            dec_mat = [corr_bit_0.real,corr_bit_1,corr_bit_0_sr.real, corr_bit_1_sr.real]
            harq_sr_dec = dec_mat.index(max(dec_mat))

            sig_power = 0
            harq_bit  = np.array(0)
            if harq_sr_dec == 0:
                harq_bit = np.array(0)
                sig_power = np.absolute(corr_bit_0)**2
                sr = 0

            if harq_sr_dec == 1:
                harq_bit = np.array(1)
                sig_power = np.absolute(corr_bit_0_sr)**2
                sr = 0

            if harq_sr_dec == 2:
                harq_bit = np.array(0)
                sig_power = np.absolute(corr_bit_0_sr)**2
                sr = 1

            if harq_sr_dec == 3:
                harq_bit = np.array(1)
                sig_power = np.absolute(corr_bit_1_sr)**2
                sr = 1

            if sig_power >= noise_power:
                dtx = 0
            else:
                dtx = 1

        elif pucch_format0_param["nHarqBit"] == 2:
            corr_bit_0_0 = np.sum(rec_symbol*expn(cmplx(0, -math.pi/6 *np.r_[:12])))
            corr_bit_0_1 = np.sum(rec_symbol*expn(cmplx(0, -4*math.pi/6 *np.r_[:12])))
            corr_bit_1_0 = np.sum(rec_symbol*expn(cmplx(0, -7*math.pi/6 *np.r_[:12])))
            corr_bit_1_1 = np.sum(rec_symbol*expn(cmplx(0, -10*math.pi/6 *np.r_[:12])))

            # Decision Matrix Without SR
            dec_mat_sr = [corr_bit_0_0.real, corr_bit_0_1.real, corr_bit_1_0.real, corr_bit_1_1.real]
            dec_mat_sr_raw = [corr_bit_0_0, corr_bit_0_1, corr_bit_1_0, corr_bit_1_1]

            harq_bit_dec_sr = dec_mat_sr.index(max((dec_mat_sr)))

            # Decision Matrix With SR
            corr_bit_0_0 = np.sum(rec_symbol)
            corr_bit_0_1 = np.sum(rec_symbol*expn(cmplx(0, -math.pi/2 *np.r_[:12])))
            corr_bit_1_0 = np.sum(rec_symbol*expn(cmplx(0, -3*math.pi/2 *np.r_[:12])))
            corr_bit_1_1 = np.sum(rec_symbol*expn(cmplx(0, -math.pi*np.r_[:12])))

            dec_mat = [corr_bit_0_0.real, corr_bit_0_1.real, corr_bit_1_0.real, corr_bit_1_1.real]
            dec_mat_raw = [corr_bit_0_0, corr_bit_0_1, corr_bit_1_0, corr_bit_1_1]

            harq_bit_dec = dec_mat.index(max((dec_mat)))

            if dec_mat_sr[harq_bit_dec_sr] >= dec_mat[harq_bit_dec]:
                sr = 1
                dec_mat = dec_mat_sr
                dec_mat_raw = dec_mat_sr_raw
                harq_bit_dec = harq_bit_dec_sr
            else:
                sr = 0

            if harq_bit_dec == 0:
                harq_bit = np.array([0, 0])
                sig_power = np.absolute(dec_mat[0]) ** 2

            if harq_bit_dec == 1:
                harq_bit = np.array([0, 1])
                sig_power = np.absolute(dec_mat_raw[1]) ** 2

            if harq_bit_dec == 2:
                harq_bit = np.array([1, 0])
                sig_power = np.absolute(dec_mat_raw[2]) ** 2

            if harq_bit_dec == 3:
                harq_bit = np.array([1, 1])
                sig_power = np.absolute(dec_mat_raw[3]) ** 2

            if sig_power > noise_power:
                dtx = 0
            else:
                dtx = 1
        else:
            if pucch_format0_param["sr"] == 1:
                corr_bit = np.sum(rec_symbol)
                sig_power = np.absolute(corr_bit) ** 2
                if sig_power > noise_power:
                    dtx = 0
                    sr = 1
                else:
                    dtx = 1
                    sr = 0
        snr = 10*np.log10(sig_power/noise_power)

        return [harq_bit, sr, snr, dtx]

    def pucch_format_1_rec(self, nGrid, n_sf_u, n_id, n_hop, pucch_format1_param, noise_power):

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
            channel_estimate += np.sum(nGrid[l + 2 * n][startPRB * 12 + np.r_[:12]]*np.conj(dmrs_sequence[n * self.N_sc_rb + np.r_[:12]]))

        channel_estimate = channel_estimate/n_dmrs_symbol/self.N_sc_rb

        # Equalizing Data / MMSE receiever

        data_estimate = 0
        for n in range(0, n_data_symbol, 1):
            data_estimate += np.sum(nGrid[l + 2 * n + 1][startPRB * 12 + np.r_[:12]]*np.conj(data_sequence[n * self.N_sc_rb + np.r_[:12]]))

        #data_estimate = data_estimate/n_data_symbol/self.N_sc_rb

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

        snr = 10*np.log10(sig_power/noise_power)

        # print(sig_power, harq_bit, dtx)

        return [harq_bit, snr, dtx]



