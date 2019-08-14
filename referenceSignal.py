import cmath
import math
import numpy as np

class ReferenceSignal:

    N_sc_rb =12
    lut_phi_mzc_six = [[-3, -1,  3,  3, -1, -3],
                       [3,  3, -1, -1,  3, -3],
                       [-3, -3, -3,  3,  1, -3],
                       [1,  1,  1,  3, -1, -3],
                       [1,  1,  1, -3, -1,  3],
                       [-3,  1, -1, -3, -3, -3],
                       [3,  1,  3, -3, -3, -3],
                       [-3, -1,  1, -3,  1, -1],
                       [-3, -1, -3,  1, -3, -3],
                       [-3, -3,  1, -3,  3, -3],
                       [-3,  1,  3,  1, -3, -3],
                       [-3, -1, -3,  1,  1, -3],
                       [1,  1,  3, -1, -3,  3],
                       [1,  1,  3,  3, -1,  3],
                       [1,  1,  1, -3,  3, -1],
                       [1,  1,  1, -1,  3, -3],
                       [-3, -1, -1, -1,  3, -1],
                       [-3, -3, -1,  1, -1, -3],
                       [-3, -3, -3,  1, -3, -1],
                       [-3,  1,  1, -3, -1, -3],
                       [-3,  3, -3,  1,  1, -3],
                       [-3,  1, -3, -3, -3, -1],
                       [1,  1, -3,  3,  1,  3],
                       [1,  1, -3, -3,  1, -3],
                       [1,  1,  3, -1,  3,  3],
                       [1,  1, -3,  1,  3,  3],
                       [1,  1, -1, -1,  3, -1],
                       [1,  1, -1,  3, -1, -1],
                       [1,  1, -1,  3, -3, -1],
                       [1,  1, -3,  1, -1, -1]]

    lut_phi_mzc_tweleve = [[-3, 1, -3, -3, -3, 3, -3, -1, 1, 1, 1, -3],
                           [-3, 3, 1, -3, 1, 3, -1, -1, 1, 3, 3, 3],
                           [-3, 3, 3, 1, -3, 3, -1, 1, 3, -3, 3, -3],
                           [-3, -3, -1, 3, 3, 3, -3, 3, -3, 1, -1, -3],
                           [-3, -1, -1, 1, 3, 1, 1, -1, 1, -1, -3, 1],
                           [-3, -3, 3, 1, -3, -3, -3, -1, 3, -1, 1, 3],
                           [1, -1, 3, -1, -1, -1, -3, -1, 1, 1, 1, -3],
                           [-1, -3, 3, -1, -3, -3, -3, -1, 1, -1, 1, -3],
                           [-3, -1, 3, 1, -3, -1, -3, 3, 1, 3, 3, 1],
                           [-3, -1, -1, -3, -3, -1, -3, 3, 1, 3, -1, -3],
                           [-3, 3, -3, 3, 3, -3, -1, -1, 3, 3, 1, -3],
                           [-3, -1, -3, -1, -1, -3, 3, 3, -1, -1, 1, -3],
                           [-3, -1, 3, -3, -3, -1, -3, 1, -1, -3, 3, 3],
                           [-3, 1, -1, -1, 3, 3, -3, -1, -1, -3, -1, -3],
                           [1, 3, -3, 1, 3, 3, 3, 1, -1, 1, -1, 3],
                           [-3, 1, 3, -1, -1, -3, -3, -1, -1, 3, 1, -3],
                           [-1, -1, -1, -1, 1, -3, -1, 3, 3, -1, -3, 1],
                           [-1, 1, 1, -1, 1, 3, 3, -1, -1, -3, 1, -3],
                           [-3, 1, 3, 3, -1, -1, -3, 3, 3, -3, 3, -3],
                           [-3, -3, 3, -3, -1, 3, 3, 3, -1, -3, 1, -3],
                           [3, 1, 3, 1, 3, -3, -1, 1, 3, 1, -1, -3],
                           [-3, 3, 1, 3, -3, 1, 1, 1, 1, 3, -3, 3],
                           [-3, 3, 3, 3, -1, -3, -3, -1, -3, 1, 3, -3],
                           [3, -1, -3, 3, -3, -1, 3, 3, 3, -3, -1, -3],
                           [-3, -1, 1, -3, 1, 3, 3, 3, -1, -3, 3, 3],
                           [-3, 3, 1, -1, 3, 3, -3, 1, -1, 1, -1, 1],
                           [-1, 1, 3, -3, 1, -1, 1, -1, -1, -3, 1, -1],
                           [-3, -3, 3, 3, 3, -3, -1, 1, -3, 3, 1, -3],
                           [1, -1, 3, 1, 1, -1, -1, -1, 1, 3, -3, 1],
                           [-3, 3, -3, 3, -3, -3, 3, -1, -1, 1, 3, -3]]

    def prime_search(self, m):
        # Initialize a list
        result = 1

        # Assume number is prime until shown it is not.

        for possiblePrime in range(2, m):

            IsPrime = True

            for num in range(2, possiblePrime):
                if possiblePrime % num == 0:
                    IsPrime = False

                if IsPrime:
                    result = possiblePrime

        return result

    def base_seq_gr_thirty_six(self, mzc, u, v):
        r_u_v = []
        nzc = self.prime_search(mzc)

        q_bar = nzc*(u+1)/31
        q = math.floor(q_bar+0.5) + v*(-1)**(2*q_bar)

        # Creating the base Sequence
        for n in range(0, mzc, 1):
            m = n % nzc
            r_u_v.append(cmath.exp(complex(0, -math.pi*q*m*(m+1)/nzc)))

        return r_u_v

    def base_seq_lt_thirty_six(self, u, mzc):
        r_u_v = []

        if mzc == 6:
            for n in range(0, 6, 1):
                r_u_v.append(cmath.exp(complex(0, -math.pi*self.lut_phi_mzc_six[u][n]/4)))
        elif mzc == 12:
            for n in range(0, 12, 1):
                r_u_v.append(cmath.exp(complex(0, -math.pi*self.lut_phi_mzc_tweleve[u][n]/4)))
        else:
            print("Reference Sequence Not Supported")

        return r_u_v

    def base_seq_eq_thirty(self, u):
        r_u_v = []

        for n in range(0, 30, 1):
            r_u_v.append(cmath.exp(complex(0, -math.pi*(u+1)*(n+1)*(n+2)/31)))

        return r_u_v

    def reference_seq(self, u, v, alpha, m, delta):

        mzc = int(m*self.N_sc_rb/(2**delta))

        r_u_v_alpha_delta = []

        if mzc >= 36:
            r_u_v = self.base_seq_gr_thirty_six(mzc, u, v)
        elif mzc < 36:
            r_u_v = self.base_seq_lt_thirty_six(u, mzc)
        elif mzc == 30:
            r_u_v = self.base_seq_lt_thirty_six(u, mzc)
        else:
            print("Unknown Reference Sequence Length")

        for n in range(0, mzc, 1):
            r_u_v_alpha_delta.append(cmath.exp(complex(0, alpha*n))*r_u_v[n])

        return r_u_v_alpha_delta

