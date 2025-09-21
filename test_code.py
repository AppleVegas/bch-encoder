from bch import BCHEncoder
import numpy as np
import matplotlib.pyplot as plt
import random

def int_input(question, check_lambda):
    while True:
        try:
            t = int(input(question + "\n>>> "))
            if not check_lambda(t):
                raise Exception()
            return t
        except:
            print("Неверный ввод!")

def count_errors(a1: list, a2: list) -> int:
    c = 0
    for i in range(0, len(a1)):
        if a1[i] != a2[i]:
            c += 1
    return c


if __name__ == "__main__":

    input_pm = int("0x" + input("Примитивный полином (в шестнадцатеричном формате):\n>>> "), 16) # Примитивный полином = 0x5B
    input_t = int_input("Число исправляемых ошибок t (целое, десятичное):", lambda x: x >= 1) # Число исправляемых ошибок t = 4

    poly_hex = hex(input_pm)
    poly_bin = bin(input_pm)
    print("Степень (m): %d\nДесятичное представление: %d\nШестнадцатеричное представление: %s\nДвоичное представление: %s\n" % (input_pm.bit_length() - 1, input_pm, poly_hex, poly_bin))

    encoder = BCHEncoder(input_pm, input_t) # Инициализоровать код

    print("Характеристики БЧХ-кода:\n(n, k, d) = (%d, %d, %d)\n" % (encoder.n, encoder.k, encoder.d))
    print("Поле Галуа:\n%s\n" % encoder.field) # Вывод поля Галуа
    print("g(x) = %s\n" % (np.binary_repr(encoder.generator))) # Вывод порождающего полинома
    
    maxword = ((1 << encoder.k) - 1)
    word = int_input("Введите слово для кодирования (целое, десятичное, %i макс.):" % maxword, lambda x: x >= 1 and x <= maxword) # Кодируемая комбинация = 370711156493
    print("Кодируемая информация: %s\n" % np.binary_repr(word))

    encoder1 = encoder.encode_non_systematic(word)
    encoder2 = encoder.encode_systematic(word)
    
    print("Результат кодирования (несистематический код): %s\n" % ''.join([str(i) for i in encoder1]),
          "                      (систематический код):   %s\n" % ''.join([str(i) for i in encoder2]), sep="")

    input_errors = int_input("Сколько ошибок вводить (целое, десятичное):", lambda x: x >= 1)

    ers = random.sample(range(0, len(encoder1)), input_errors)
    for i in ers:
        encoder1[i] ^= 1
        encoder2[i] ^= 1

    print("Комбинация с ошибками (несистематический код): %s\n" % ''.join([str(i) for i in encoder1]),
          "                      (систематический код):   %s\n" % ''.join([str(i) for i in encoder2]), sep="")

    fixed1 = encoder.decode(encoder1, False)
    fixed2 = encoder.decode(encoder2, True)

    print("Исправленная комбинация (несистематический код): %s\n" % ''.join([str(i) for i in fixed1[0]]), 
          "                          Позиции ошибочных бит: %s\n" % fixed1[1],
          "                          (систематический код): %s\n" % ''.join([str(i) for i in fixed2[0]]), 
          "                          Позиции ошибочных бит: %s\n" % fixed2[1], sep="")

    word_poly = encoder.field.int_to_poly(word)
    if word_poly != fixed1[0]:
        print(f"Несистематическкий код декодировался с {count_errors(word_poly, fixed1[0])} ошибками!")
    
    if word_poly != fixed2[0]:
        print(f"Систематическкий код декодировался с {count_errors(word_poly, fixed2[0])} ошибками!")

    # Построение графика BER

    L = 1000 # Число испытаний
    num_points = 50 # Число точек
    percent_maxerrors = 100

    test_range = np.arange(0, percent_maxerrors + (percent_maxerrors//num_points), percent_maxerrors//num_points)
    BER1 = []
    BER2 = []
    for probability in test_range:
        print("Построение графика... %d/%d" % (len(BER1), num_points))
        error_probability = (probability/100)*1
        errors1 = 0
        errors2 = 0
        for i in range(0, L):
            encoder1 = encoder.encode_non_systematic(word)
            encoder2 = encoder.encode_systematic(word)

            for bit_i in range(0, encoder.n):
                is_error = np.random.choice(np.arange(0, 2), p=[1 - error_probability, error_probability])
                if is_error:
                    encoder1[bit_i] = not encoder1[bit_i]
                    encoder2[bit_i] = not encoder2[bit_i]
            
            fixed1 = encoder.decode(encoder1, False)
            fixed2 = encoder.decode(encoder2, True)

            if word_poly != fixed1[0]:
                errors1 += count_errors(word_poly, fixed1[0])
            if word_poly != fixed2[0]:
                errors2 += count_errors(word_poly, fixed2[0])

        BER1.append(errors1/(L * encoder.k))
        BER2.append(errors2/(L * encoder.k))

    plt.yscale("log")
       
    plt.plot(test_range, BER1)
    plt.plot(test_range, BER2)
    plt.legend (('Несистематическое кодирование', 'Систематическое кодирование'))
    plt.xlabel('Вероятность ошибки, %')
    plt.ylabel('Коэффициент битовых ошибок BER')
    plt.grid(True, which="both")
    plt.gca().invert_xaxis()
    plt.title("БЧХ-код (%d, %d, %d)" % (encoder.n, encoder.k, encoder.d))
    plt.show()
