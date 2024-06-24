from bch import BCHEncoder
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import random
#import sys
#np.set_printoptions(threshold=sys.maxsize)

def int_input(question, check_lambda):
    while True:
        try:
            t = int(input(question))
            if not check_lambda(t):
                raise Exception()
            return t
        except:
            print("Wrong input.")

if __name__ == "__main__":
    # Input of primitive polynomial
    input_pm = int("0x" + input("Primitive polynomial (hex): "), 16)
    input_t = int_input("Correctable error amount t (int): ", lambda x: x >= 1)

    poly_hex = hex(input_pm)
    poly_bin = bin(input_pm)
    print("Power (m): %d\nInt: %d\nHex: %s\nBinary: %s\n" % (input_pm.bit_length() - 1, input_pm, poly_hex, poly_bin))

    encoder = BCHEncoder(input_pm, input_t) # Initialize encoder

    print("BCH Characteristics:\n(n, k, t) = (%d, %d, %d)\nd = %d\n" % (encoder.n, encoder.k, encoder.t, encoder.d))
    print("Galois Field:\n%s\n" % encoder.field) # Printing out field
    print("g(x) = %s\n" % (np.binary_repr(encoder.generator))) # Printing generator polynomial
    print("Generator matrix\nG = \n%s\n\nParity-check matrix\nH = \n%s\n\nParity-check matrix transposed\nH^T = \n%s\n" % (encoder.G, encoder.H, encoder.HT))

    word = 0b10100111001 # 1337
    maxword = ((1 << encoder.k) - 1)
    word = int_input("Input word to incode (int, %i max)\n" % maxword, lambda x: x >= 1 and x <= maxword)
    print("Data to encode: %s\n" % np.binary_repr(word))

    encoder1 = encoder.encode_non_systematic(word)
    encoder2 = encoder.encode_systematic(word)
    
    print("Encoding result (non-systemical): %s\n" % encoder1,
          "Encoding result (systemical):     %s\n" % encoder2, sep="")

    input_errors = int_input("Input number of errors to be added (int)\n", lambda x: x >= 1)

    ers = random.sample(range(0, len(encoder1)), input_errors)
    for i in ers:
        encoder1[i] ^= 1
        encoder2[i] ^= 1

    print("Decoding result (non-systemical): %s, Errored bits: %s\n" % encoder.decode(encoder1, False),
          "Decoding result (systemical):     %s, Errored Bits: %s\n" % encoder.decode(encoder2, True), sep="")

    frames_amount = 64 # number of frames
    sample_data = np.array([encoder.field.int_to_poly((1 << encoder.n) - 1)] * frames_amount) 
    test_range = np.arange(0, 50)
    BER = []
    FER = []
    for probability in test_range:
        error_probability = (probability/100)*1
        bit_errors = 0
        frame_errors = 0
        for frame in sample_data:
            frame_bit_errors = 0
            for bit in frame:
                is_error = np.random.choice(np.arange(0, 2), p=[1 - error_probability, error_probability])
                if is_error:
                    frame_bit_errors += 1
            if frame_bit_errors > encoder.t:
                frame_errors += 1
            bit_errors += frame_bit_errors
            
        FER.append(frame_errors/frames_amount)
        BER.append(bit_errors/(encoder.n*frames_amount))

    test_range_lin = np.linspace(test_range.min(), test_range.max(), 200) 
    
    BER = savgol_filter(BER, 10, 1)
    FER = savgol_filter(FER, 10, 1)

    ber_smooth = np.poly1d(np.polyfit(test_range, BER, 2))
    fer_smooth = np.poly1d(np.polyfit(test_range, FER, 2))

    plt.yscale("log")
       
    plt.plot(test_range_lin, ber_smooth(test_range_lin), 'red')
    plt.plot(test_range_lin, fer_smooth(test_range_lin), 'blue')

    plt.legend (('Bit Error Rate (BER)', 'Frame Error Rate (FER)'))
    plt.xlabel('Error Probability, %')
    plt.ylabel('Error Rate')
    plt.gca().invert_xaxis()
    plt.grid(True, which="both")
    plt.show()
