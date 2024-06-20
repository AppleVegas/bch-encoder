import numpy as np
import matplotlib as plt
import galois

# Generator for generating galois field. Equals to 2, because we're generating galois field for base 2.
# Actually is p in GF(p**m)
GENERATOR = 0b10

def printbin(n, l):
    print(bin(n)[2:].zfill(l))

# Input of primitive polynomial
input_pm = int("0x" + '83', 16)#input("Primitive polynomial (hex): "), 16)
input_t = 5 #input("Correctable error count t (int): ")
input_encodetype = 2
while input_encodetype is None:
    try:
        t = int(input("Encoding type\n1 - Polynomial (non-systematic)\n2 - Matrix (systematic)\n"))
        if t < 1 or t > 2:
            raise Exception()
        input_encodetype = t
    except:
        print("ffs...")

poly_hex = hex(input_pm)
poly_bin = bin(input_pm)
print("Power (m): %d\nInt: %d\nHex: %s\nBinary: %s\n" % (input_pm.bit_length() - 1, input_pm, poly_hex, poly_bin))


class GF():
    def __init__(self, poly: int) -> 'GF':
        self.p = poly # Primitive polynomial used in generating the field and multiplication
        self.m = poly.bit_length() - 1 # Power m in GF(p**m)
        self.field = [] # Array for storing the field
        self.lookup = {} # Dictionary for reverse lookup of index of field item
        self.generate()
    
    def generate(self) -> None:
        a = 1
        for i in range(0, 2**self.m):
            self.field.append(a)
            self.lookup[a] = i
            a = self._mul(a, GENERATOR)
    
    def a(self, index) -> 'int':
        return self.field[index % (2**self.m - 1)]
    
    def look(self, a) -> 'int':
        return self.lookup[a]

    def add(self, a: int, b: int) -> 'int':
        return a ^ b

    def pow(self, a: int, pow: int) -> 'int':
        return self.field[(self.lookup[a] * pow) % (2**self.m - 1)]
    
    def mul(self, a: int, b: int) -> 'int':
        return self.field[(self.lookup[a] + self.lookup[b]) % (2**self.m - 1)] if a != 0 and b != 0 else 0
    
    def div(self, a: int, b: int) -> 'int':
        if b == 0:
            raise ZeroDivisionError()
        return self.field[(self.lookup[a] - self.lookup[b]) % (2**self.m - 1)] if a != 0 else 0
    
    def _mod(self, a: int, b: int):
        while a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        return a
    
    def _mul(self, a: int, b: int) -> 'int':
        r = 0
        while b != 0:
            if (b & 1):
                r ^= a
            if (a & 2**(self.m - 1)):
                a = (a << 1) ^ self.p
            else:
                a <<= 1
            b >>= 1
        return r
    
    def poly_mul(self, a: list, b: list) -> 'list':
        r = [0] * (len(a) + len(b) - 1)
        for i, coef_a in enumerate(a):
            for j, coef_b in enumerate(b):
                r[i + j] ^= self._mul(coef_a, coef_b)
        return r
    
    def poly_eval(self, a: list, x: int) -> 'int':
        r = a[0]
        for c in a[1:]:
            r = self._mul(r, x) ^ c
        return r
    
    def poly_add(self, a: list, b: list):
        l = len(a) - len(b)
        if l > 0:
            b = [0]*l + b
        elif l < 0:
            a = [0]*(-l) + a
        return [(c1 ^ c2) for c1, c2 in zip(a, b)]

    def poly_scale(self, a: list, s: int) -> 'int':
        return [self._mul(c, s) for c in a]
    
    def poly_to_int(self, poly: list) -> 'int':
        #poly = [i % 2 for i in poly]
        r = 0
        for coef in poly:
            r = (r << 1) | coef
        return r
    
    def int_to_poly(self, integer: int) -> 'list':
        l = integer.bit_length()
        r = [0] * l
        for i in range(0, l):
            if integer & (1 << i):
                r[(l - 1) - i] = 1
                continue
            r[(l - 1) - i] = 0 
        return r # [int(x) for x in np.binary_repr(integer)]
    
    def min_poly(self, root):
        min_poly = [1]  # Let f(x) be 1

        seen_roots = set()  # Separating cyclotomic classes 

        for _ in range(1, 2**self.m): # i <= 2^m
            if root in seen_roots:
                break

            seen_roots.add(root)
            min_poly = self.poly_mul(min_poly, [1, root])  # Multiply by (x - root)
            root = self.pow(root, 2)
    
        return min_poly

    # String representaion of the field for printing
    def __str__(self) -> str:
        return ('GF(2^%d) = \n' % self.m) + '\n'.join([('a^%d = ' % i) + bin(s)[2:].zfill(self.m) for i,s in enumerate(self.field)])

class BCHEncoder():
    def __init__(self) -> 'BCHEncoder':
        # Creating field
        self.field = GF(input_pm) # ideally polynomials should be generated ig, using input polys from table rn
        self.m = input_pm.bit_length() - 1
        self.n = 2**self.m - 1
        self.t = input_t
        self.d = 2*self.t + 1
        self.generator_poly = self.get_generator()
        self.generator = self.field.poly_to_int(self.generator_poly)
        self.k = self.n - (self.generator.bit_length() - 1)
        self.r = self.n - self.k
        self.G = None # Generator matrix
        self.H = None # Parity-check matrix
        self.HT = None # Parity-check matrix transposed
        self.build_matrix()
    
    def get_generator(self) -> 'list':
        g_x = [1]

        seen_roots = set()  # Separating cyclotomic classes 

        for i in range(1, 2*self.t + 1): # i <= 2t
            root = self.field.a(i)

            if root in seen_roots:
                continue

            f_x = [1]

            for _ in range(1, 2**self.m): # i <= 2^m
                if root in seen_roots:
                    break

                seen_roots.add(root)
                f_x = self.field.poly_mul(f_x, [1, root])  # Multiply by (x - root)
                root = self.field.pow(root, 2)

            g_x = self.field.poly_mul(g_x, f_x)
        return g_x

    def build_matrix(self):
        R = np.zeros(self.k, np.int64)
        R[self.k - 1] = 1 << self.r # Beginning from XORing x^n
        for k in range(self.k - 1, -1, -1): # Reversibly filling a matrix so we don't have to reverse it later
            if (R[k] >> self.r) & 1:
                R[k] ^= self.generator
            if k:
                R[k - 1] = R[k] << 1
        # this^ could have been replaced with self.field._mod(1 << k, self.generator)
        
        R = (((R[:,None] & (1 << np.arange(self.r))[::-1])) > 0).astype(int)
        self.G = np.hstack((np.identity(self.k, int),R))
        self.H = np.hstack((R.transpose(),np.identity(self.r, int)))
        self.HT = np.vstack((R,np.identity(self.r, int)))
        
        
    def encode_systematic(self, data: int) -> 'int':
        t = (data << self.r) 
        return self.field.int_to_poly( t ^ ( self.field._mod(t, self.generator)) )

    def encode_non_systematic(self, data: int) -> 'list':
        encoded = self.field.poly_mul(self.field.int_to_poly(data), self.generator_poly)
        return encoded
    
    def get_syndromes(self, data_v: list) -> 'list':
        S = [self.field.poly_eval(data_v, self.field.a(i)) for i in range(1, 2*self.t + 1)]

        ''' this also can be used but data_v has to be reversed
        data_v = data_v[::-1]
        S = [0] * (2*self.t)
        for i in range(1, 2*self.t + 1): # i <= 2t
            for j in range(0, len(data_v)): # i <= v
                if data_v[j] != 0:
                    S[i - 1] ^= self.field.a(i * j)
        '''
        return S

    def find_error_locator(self, synd): # Berlekamp-Massey algorithm
        err_loc = [1]
        old_loc = [1]

        for i in range(len(synd)):
            old_loc.append(0)
            delta = synd[i]
            for j in range(1, len(err_loc)):
                delta ^= self.field.mul(err_loc[- (j + 1) ], synd[i - j])
            if delta != 0:
                if len(old_loc) > len(err_loc):
                    new_loc = self.field.poly_scale(old_loc, delta)
                    old_loc = self.field.poly_scale(err_loc, self.field.div(1, delta))
                    err_loc = new_loc
                err_loc = self.field.poly_add(err_loc, self.field.poly_scale(old_loc, delta))
        return err_loc

    def find_errors(self, locator_poly, length):
        err_pos = []
        for i in range(length):
            if self.field.poly_eval(locator_poly, self.field.a(i)) == 0:
                err_pos.append((length - i) % length)
        return err_pos
    
    def decode(self, data_v: list) -> 'list':
        fixed_data = data_v
        length = len(data_v)
        syndromes = self.get_syndromes(data_v)

        for v in range(self.t, 0, -1): # initially v = t 
            s_matrix = np.zeros((v, v), dtype=int)
            for i in range(v):
                for j in range(v):
                    s_matrix[i, j] = syndromes[i + j]

            det = np.linalg.det(s_matrix)

            if det == 0: # if matrix determinator is 0, then v = v - 1
                continue
            
            loc = self.find_error_locator(syndromes)
            errors = self.find_errors(loc, length)
            print(errors)
            for error in errors:
                fixed_data[(length - 1) - error] ^= 1
            break
        
        syndromes = self.get_syndromes(fixed_data)

        if max(syndromes) != 0:
            raise Exception("Too many errors.")
        
        return fixed_data


encoder = BCHEncoder()

# Printing out field
print(encoder.field)
print(bin(encoder.generator)) # Printing generator polynomial
print("Generator matrix\nG = \n%s\n\nParity-check matrix\nH = \n%s\n\nParity-check matrix transposed\nH^T = \n%s\n" % (encoder.G, encoder.H, encoder.HT))

word = 0b1011101
encoder1 = encoder.encode_non_systematic(word)
encoder2 = encoder.encode_systematic(word)

#encoder1 ^= (0b001110000000000)
#encoder2 ^= (0b1001010)

print(encoder1, "encode non sys")
print(encoder2, "encode sys\n")
print(encoder.decode([1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]), "decode non sys")
print(encoder.decode([1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0]), "decode sys")

#TODO: cleanup, decoding