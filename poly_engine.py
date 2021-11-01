from cmath import sqrt
import numpy as np


class Polynomial:
    def __init__(self, *args, name='f'):
        self.coefficients = list(args)
        self.name = name
        self.degree = len(self.coefficients) - 1

    def __str__(self):
        printable_list = [f'{self.coefficients[0]}'] if self.coefficients[0] != 0 else []
        for i, c in enumerate(self.coefficients):
            if i == 1 and c != 0:
                if c != 1:
                    printable_list.append(f'{c}x + ')
                else:
                    printable_list.append(f'x + ')
            if i > 1 and c != 0:
                if c != 1:
                    printable_list.append(f'{c}x^{i} + ')
                else:
                    printable_list.append(f'x^{i} + ')

        result = f'{self.name}(x) = ' + ''.join(printable_list[::-1])
        if self.coefficients[0] == 0:
            result = result[:-3]
        return result

    def __len__(self):
        return self.degree + 1

    def __getitem__(self, n):
        if n > self.degree:
            raise IndexError(f'{n} is greater than the degree of the polynomial.')
        else:
            return self.coefficients[n]

    def __eq__(self, other):
        if self.coefficients == other.coefficients:
            return True
        return False

    def __add__(self, other):
        g, l = sorted((self, other), key=lambda x: len(x.coefficients), reverse=True)
        new_coefficients = []
        for i in range(len(g)):
            if i <= l.degree:
                new_coefficients.append(g.coefficients[i] + l.coefficients[i])
            else:
                new_coefficients.append(g.coefficients[i])
        new_poly = Polynomial(*new_coefficients, name='h')
        return new_poly

    def __sub__(self, other):
        other.coefficients = [-coeff for coeff in other.coefficients]
        return self + other

    def __mul__(self, other):
        pass

    def __floordiv__(self, other):
        if self.degree < other.degree:
            raise ValueError('The degree of the numerator should be at least as great '
                             'as the degree of the denominator.')
        result, _ = np.polynomial.polynomial.polydiv(self.coefficients, other.coefficients)
        return Polynomial(*result, name='q')

    @property
    def prime(self):
        new_coefficients = []
        for i, coeff in enumerate(self.coefficients[1:]):
            new_coefficients.append((i+1)*coeff)
        return Polynomial(*new_coefficients, name=self.name + "'")

    def diff(self, n=1):
        if n < 0 or type(n) != int:
            raise ValueError('You need to specify a non-negative whole number.')
        if n == 0:
            return self
        elif n == 1:
            return self.prime
        else:
            result = self
            for _ in range(n):
                result = result.prime
        return Polynomial(*result.coefficients, name=self.name + f'^({n})')

    def __call__(self, number):
        result = self.coefficients[0]
        for i, coeff in enumerate(self.coefficients[1:]):
            result += coeff * number ** (i + 1)
        return result

    @property
    def zeros(self):
        return list(np.roots(self.coefficients[::-1]))

