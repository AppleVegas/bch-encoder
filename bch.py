import numpy as np
import matplotlib as plt

# Generator for generating galois field. Equals to 2, because we're generating galois field for base 2.
# Actually is p in GF(p**m)
GENERATOR = 0b10

def printbin(n, l):
    print(bin(n)[2:].zfill(l))

# Input of primitive polynomial
input_pm = int("0x" + '13', 16)#input("Primitive polynomial (hex): "), 16)
input_t = 2 #input("Correctable error count t (int): ")
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

def _mod(a: int, b: int):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a

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
        return self.field[(self.lookup[a] + self.lookup[b]) % (2**self.m - 1)]
    
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
        #r = [i % 2 for i in r]
        return r
    
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
            r[i] = 0 
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
        self.build_matrix()
    
    def get_generator(self) -> 'list':
        g_x = [1]

        seen_roots = set()  # Separating cyclotomic classes 

        for i in range(1, 2*self.t + 1): # i <= 2t
            root = self.field.a(i)

            if root in seen_roots:
                continue

            #print(f"Visiting a^{i} = {bin(field.a(i))[2:].zfill(m)}")

            f_x = [1]

            for _ in range(1, 2**self.m): # i <= 2^m
                if root in seen_roots:
                    break

                seen_roots.add(root)
                f_x = self.field.poly_mul(f_x, [1, root])  # Multiply by (x - root)
                root = self.field.pow(root, 2)
            
            #print(bin(self.poly_to_int(f_x)))

            g_x = self.field.poly_mul(g_x, f_x)
        return g_x

    def build_matrix(self) -> 'np.ndarray':
        R = np.zeros(self.k, np.int64)
        R[self.k - 1] = 1 << self.r # Beginning from XORing x^n
        for k in range(self.k - 1, -1, -1): # Reversibly filling a matrix so we don't have to reverse it later
            if (R[k] >> self.r) & 1:
                R[k] ^= self.generator
            if k:
                R[k - 1] = R[k] << 1
        # ^ could have been replaced with _mod(1 << k, self.generator)

        for k in R:
            print(np.binary_repr(k).zfill(self.r))
        
        R = (((R[:,None] & (1 << np.arange(self.r))[::-1])) > 0).astype(int)
        self.G = np.hstack((np.identity(self.k, int),R))
        self.H = np.hstack((R.transpose(),np.identity(self.r, int)))
        
        
    def encode_systematic(self, data: int) -> 'int':
        t = self.field.poly_to_int(self.field.poly_mul(self.field.int_to_poly(2 << self.r - 1), self.field.int_to_poly(data)))
        return int( t ^ ( _mod(t, self.generator)) )

    def encode_non_systematic(self, data: int) -> 'int':
        encoded = self.field.poly_mul(self.field.int_to_poly(data), self.generator_poly)
        return self.field.poly_to_int(encoded)

    def decode(self, data: int) -> 'int':
        pass

encoder = BCHEncoder()

# Printing out field
print(encoder.field)
print(bin(encoder.generator)) # Printing generator polynomial
print("Generator matrix\nG = \n%s\n\nParity-check matrix\nH = \n%s\n" % (encoder.G, encoder.H))
print(bin(encoder.encode_non_systematic(0b1110001)))
print(bin(encoder.encode_systematic(0b1110001)))
#TODO: cleanup, decoding