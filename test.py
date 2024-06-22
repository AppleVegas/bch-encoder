from bch import BCHEncoder
import numpy as np
import matplotlib as plt
import random
#import sys
#np.set_printoptions(threshold=sys.maxsize)

if __name__ == "__main__":
    # Input of primitive polynomial
    input_pm = int("0x" + '83', 16)#input("Primitive polynomial (hex): "), 16)
    input_t = 2 #input("Correctable error count t (int): ")
    input_encodetype = 2
    while input_encodetype is None:
        try:
            t = int(input("Encoding type\n1 - Polynomial (non-systematic)\n2 - Matrix (systematic)\n"))
            if t < 1 or t > 2:
                raise Exception()
            input_encodetype = t
        except:
            print("Wrong.")

    poly_hex = hex(input_pm)
    poly_bin = bin(input_pm)
    print("Power (m): %d\nInt: %d\nHex: %s\nBinary: %s\n" % (input_pm.bit_length() - 1, input_pm, poly_hex, poly_bin))

    encoder = BCHEncoder(input_pm, input_t) # Initialize encoder

    print("BCH Characteristics:\n(n, k, t) = (%d, %d, %d)\nd = %d\n" % (encoder.n, encoder.k, encoder.t, encoder.d))
    print("Galois Field:\n%s\n" % encoder.field) # Printing out field
    print("g(x) = %s" % (np.binary_repr(encoder.generator))) # Printing generator polynomial
    print("Generator matrix\nG = \n%s\n\nParity-check matrix\nH = \n%s\n\nParity-check matrix transposed\nH^T = \n%s\n" % (encoder.G, encoder.H, encoder.HT))

    word = 0b101101#(1 << encoder.k - 15) - (1 << encoder.k - random.randint(1, encoder.k - 20))
    print("Data to encode: %s\n" % np.binary_repr(word))

    encoder1 = encoder.encode_non_systematic(word)
    encoder2 = encoder.encode_systematic(word)
    
    print("Encoding result (non-systemical): %s\n" % encoder1,
          "Encoding result (systemical):     %s\n" % encoder2, sep="")

    errors = 2 # Number of errors to be added
    ers = random.sample(range(0, len(encoder1)), errors)
    for i in ers:
        encoder1[i] ^= 1
        encoder2[i] ^= 1
    
    print("Decoding result (non-systemical): %s, Errored bits: %s\n" % encoder.decode(encoder1, False),
          "Decoding result (systemical):     %s, Errored Bits: %s\n" % encoder.decode(encoder2, True), sep="")
