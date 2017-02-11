"""
Code to find x

Gabriel Cardoso 16-03-16
"""

import scipy
from scipy import special
import math
import pickle
import time


def calc_jp(L, delta_eta, n_dict):
    """
    Compute the Bessel Function, it sums the bessel function of the same order, but different moduli
    P (k) = A*k^-3
    :param l: the order of the Bessel function
    :param delta_eta:
    :param n_dict: the dictionary with the moduli and vectors
    :return: the vector x, where each component is the sum of the same order, but different moduli
    """
    partial_time = time.time()
    matrix_jp = []
    for l in range(L-1):
        matrix_jp.append([])
    A = 1.0
    for l in range(2, L+1):
        vector_jp = []
        for k in sorted(n_dict.keys()):
            j = scipy.special.sph_jn(l, 2*math.pi*k*delta_eta)
            p = math.sqrt(A) * math.pow(2*math.pi*k*delta_eta, -1.5)  # Harrizon-Zeldovich powerspec
            vector_jp.append(j[0][l]*p)
        matrix_jp[l-2] = vector_jp
    print('calc_jp: OK in {time} s'.format(time = time.time() - partial_time))
    return matrix_jp


def save_jp(filename, x):
    """
    :param filename:
    :param x:
    :return:
    """
    partial_time = time.time()
    with open(filename, 'wb') as f:
        pickle.dump(x, f, protocol=-1)
    print('save_jp: OK in {time} s'.format(time = time.time() - partial_time))


def load_jp(filename):
    """
    :param filename:
    :return:
    """
    partial_time = time.time()
    with open(filename, 'rb') as f:
        x = pickle.load(f)
    print('load_jp: OK in {time} s'.format(time = time.time() - partial_time))
    return x
