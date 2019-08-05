import numpy as npy


class CSequence:

    def generate_c_sequence(self, cinit, c_len):
        c_seq = []
        x2 = []
        x2_temp = [int(x) for x in bin(cinit)[2:]]

        if (len(x2_temp)) > 31:
           print("Error: Cinit sequence is greater than 31 bits")
        else:
            x2 = [0]*(31-len(x2_temp))
            for n in range(0, len(x2_temp), 1):
                x2.append(x2_temp[n])

            x1 = [0]*30+[1]
            for n in range(0, 1600, 1):
                x2_update = (x2[27] + x2[28] + x2[29] + x2[30]) % 2
                x2 = [x2_update] + x2[:-1]
                x1_update = (x1[27]+x1[30]) % 2
                x1 = [x1_update] + x1[:-1]

            for n in range(0, c_len, 1):
                c_seq.append((x1[30]+x2[30]) % 2)
                x2_update = (x2[27] + x2[28] + x2[29] + x2[30]) % 2
                x2 = [x2_update] + x2[:-1]
                x1_update = (x1[27]+x1[30]) % 2
                x1 = [x1_update] + x1[:-1]
            return c_seq
        return [0]

    def cinit_calc(self, c_seq, c_seq_len):
        c_init = 0
        for n in range(0, c_seq_len, 1):
            c_init = c_init+(2**n)*c_seq[n]
        return c_init







