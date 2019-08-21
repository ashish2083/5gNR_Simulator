import numpy as np


class encoder():
    rm_table = np.array([[1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                         [1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1],
                         [1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1],
                         [1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1],
                         [1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1],
                         [1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1],
                         [1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
                         [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1],
                         [1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1],
                         [1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1],
                         [1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1],
                         [1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1],
                         [1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
                         [1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
                         [1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1],
                         [1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1],
                         [1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0],
                         [1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0],
                         [1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0],
                         [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                         [1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
                         [1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1],
                         [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1],
                         [1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1],
                         [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0],
                         [1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1],
                         [1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0],
                         [1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0],
                         [1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0],
                         [1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                         [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    def reed_muller_encoder(self, bits, len):

        bits = np.array(bits)
        if len < 3 and len > 11:
            print("Encoder:ReedMuller, Ivnvalid bit length")

        # Appending Zeros and repeating 32 times
        zer_t_app = 11-len

        if zer_t_app > 0:
            bits = np.concatenate((bits, np.zeros(zer_t_app)))

        bits = np.tile(bits, (32, 1))
        output = np.sum((bits*self.rm_table), 1) % 2
        return output

    def reed_muller_decoder_soft(self, soft_bits, len):
        # Soft bits bit=0 -> + ; bit ==1 -> -1
        # Correlate with RM_Table to get bits
        dec_mat = []
        for i in range(0, 2**len, 1):
            bits = np.binary_repr(i, len)
            bits = np.array([int(x) for x in bits])
            enc_bits = (1-2*self.reed_muller_encoder(bits, len))
            corr_bits = np.dot(enc_bits, soft_bits)
            dec_mat.append(corr_bits)

        dec_bit = dec_mat.index(max(dec_mat))
        bits = np.binary_repr(dec_bit, len)
        bits = np.array([int(x) for x in bits])
        return bits




