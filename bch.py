import numpy as np
import matplotlib as plt

# Generator for generating galois field. Equals to 2, because we're generating galois field for base 2.
# Actually is p in GF(p**m)
GENERATOR = 0b10

def printbin(n, l):
    print(bin(n)[2:].zfill(l))

# Input of primitive polynome
input_pm = int("0x" + '83', 16)#input("Primitive polynome (hex): "), 16)
poly_hex = hex(input_pm)
poly_bin = bin(input_pm)
m = input_pm.bit_length() - 1
n = 2**m - 1
t = 1
d = 2*t + 1
print("Power (m): %d\nInt: %d\nHex: %s\nBinary: %s\n" % (m, input_pm, poly_hex, poly_bin))

class GF():
    def __init__(self, poly: int) -> 'GF':
        self.p = poly # Primitive polynome used in generating the field and multiplication
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
    
    def poly_to_int(self, poly):
        r = 0
        for coef in poly:
            r = (r << 1) | coef
        return r
    
    def min_poly(self, root):
        min_poly = [1]  # Let f(x) be 1

        seen_roots = set()  # Separating cyclotomic classes 

        for _ in range(1, 2**m): # i <= 2^m
            if root in seen_roots:
                break

            seen_roots.add(root)
            min_poly = self.poly_mul(min_poly, [1, root])  # Multiply by (x - root)
            root = self.pow(root, 2)
    
        return min_poly
    
    def get_generator(self):
        g_x = [1]

        seen_roots = set()  # Separating cyclotomic classes 

        for i in range(1, 2*t + 1): # i <= 2t
            root = field.a(i)

            if root in seen_roots:
                continue

            #print(f"Visiting a^{i} = {bin(field.a(i))[2:].zfill(m)}")

            f_x = [1]

            for _ in range(1, 2**m): # i <= 2^m
                if root in seen_roots:
                    break

                seen_roots.add(root)
                f_x = self.poly_mul(f_x, [1, root])  # Multiply by (x - root)
                root = field.pow(root, 2)

            #print(bin(self.poly_to_int(f_x)))

            g_x = self.poly_mul(g_x, f_x)
        return self.poly_to_int(g_x)

    # String representaion of the field for printing
    def __str__(self) -> str:
        return ('GF(2^%d) = \n' % self.m) + '\n'.join([('a^%d = ' % i) + bin(s)[2:].zfill(self.m) for i,s in enumerate(self.field)])

# Creating field
field = GF(input_pm)
# Printing out field
print(field)

print(bin(field.get_generator())) # Printing generator polynome
#TODO: cleanup, encoding decoding