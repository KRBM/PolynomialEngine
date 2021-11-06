from cmath import sqrt
import numpy as np
from scipy.special import binom
import matplotlib.pyplot as plt


class Polynomial:
    def __init__(self, *args, name='f'):
        self.coefficients = list(args)
        self.name = name
        self.degree = len(self.coefficients) - 1

    def __str__(self):
        printable_list = []
        self.coefficients = list(map(lambda r: round(r, 2), self.coefficients))

        for i, c in enumerate(self.coefficients):

            if c == 0:
                continue

            sign = np.sign(c)
            c = abs(c)
            if sign == -1:
                sign = ' - '
            else:
                sign = ' + '

            if i == 0:
                printable_list.append(f'{sign}{c}')

            elif i == 1:
                if c not in (1, -1):
                    printable_list.append(f'{sign}{c}x')
                else:
                    printable_list.append(f'{sign}x')

            elif i > 1:
                if c not in (1, -1):
                    printable_list.append(f'{sign}{c}x^{i}')
                else:
                    printable_list.append(f'{sign}x^{i}')
        if sign.strip() == '+':
            printable_list[-1] = printable_list[-1][3:]
        else:
            printable_list[-1] = printable_list[-1][1] + printable_list[-1][3:]

        result = f'{self.name}(x) = ' + ''.join(printable_list[::-1])

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
        result = np.polymul(self.coefficients[::-1], other.coefficients[::-1])
        return Polynomial(*result[::-1], name='h')

    def __truediv__(self, number):
        assert type(number) == int or type(number) == float
        coefficients = list(map(lambda x: x/number, self.coefficients))
        return Polynomial(*coefficients)

    def __floordiv__(self, other):
        if self.degree < other.degree:
            raise ValueError('The degree of the numerator should be at least as great '
                             'as the degree of the denominator.')
        result, _ = np.polynomial.polynomial.polydiv(self.coefficients, other.coefficients)
        return Polynomial(*result, name='q')

    def negating_shift(self, alpha):
        """Returns p(alpha - x) as an expanded polynomial object"""
        result = Polynomial(0)
        for n, c in enumerate(self.coefficients):
            if n == 0:
                result += Polynomial(c)
            else:
                coefficients = []
                for k in range(n + 1):
                    coeff = c * (-1)**k * alpha**(n - k) * binom(n, k)
                    coefficients.append(coeff)
                result += Polynomial(*coefficients)
        return result

    def make_reflecsive(self, alpha):
        """For any polynomial f, this method returns (f(x) + f(alpha - x)) / 2"""
        return (self + self.negating_shift(alpha)) / 2

    @property
    def prime(self):
        new_coefficients = []
        for i, coeff in enumerate(self.coefficients[1:]):
            new_coefficients.append((i+1)*coeff)
        return Polynomial(*new_coefficients, name=self.name + "'")

    @property
    def local_extrema(self):
        return self.prime.zeros

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
            if n == 2:
                return Polynomial(*result.coefficients, name=self.name + "''")
        return Polynomial(*result.coefficients, name=self.name + f'^({n})')

    def __call__(self, number):
        result = self.coefficients[0]
        for i, coeff in enumerate(self.coefficients[1:]):
            result += coeff * number ** (i + 1)
        return result

    @property
    def zeros(self):
        return list(np.roots(self.coefficients[::-1]))

    def display_graph(self, size=(-100, 100)):
        x = np.arange(size[0], size[1], 0.1)
        y = np.array([self(t) for t in x])
        plt.grid()
        plt.title(f'{self.__str__()}')
        plt.plot(x, y)
        plt.show()

    def display_zeros(self):
        x_points = [zero.real for zero in self.zeros]
        y_points = [zero.imag for zero in self.zeros]
        plt.scatter(x_points, y_points)
        plt.grid()
        plt.title(f'Zeros of {self.__str__()}')
        plt.show()

