"""
Code to find y

Gabriel Cardoso 16-03-16
"""

import scipy
from scipy import special
import math
import numpy
import pickle
import healpy
import time


def calc_y(L, n_dict):
    """
    REDEFINE
    :param l: the degree of the Spherical Harmonic
    :param n_dict: the dictionary with the moduli and vectors
    :return: REDEFINE
    """
    partial_time = time.time()
    matrix_y = dict()
    for l in range(2, L+1):
        for m in range(l+1):
            index = healpy.sphtfunc.Alm.getidx(L, l, m)
            vector_y_k = dict()
            for k in n_dict.keys():
                angles = numpy.transpose(n_dict[k])
                y = scipy.special.sph_harm(m, l, angles[0], angles[1])
                vector_y_k.update({k: y})
            matrix_y.update({index: vector_y_k})
            # print('calc_y: y{yl},{ym} OK'.format(yl=l, ym=m))
    print('calc_y: OK in {time} s'.format(time = time.time() - partial_time))
    return matrix_y


def calc_y_pi_twisted(L, n_dict):
    """
    REDEFINE
    :param L: the degree of the Spherical Harmonic
    :param n_dict: the dictionary with the moduli and vectors
    :return: REDEFINE
    """
    partial_time = time.time()
    matrix_y = dict()
    for l in range(2, L + 1):
        for m in range(l + 1):
            index = healpy.sphtfunc.Alm.getidx(L, l, m)
            vector_y_k = dict()
            for k in n_dict.keys():
                vector_y = []
                for angles in n_dict[k]:
                    if angles[0] == math.pi/2:
                        if round((k*math.cos(angles[1]))) % 2:
                            y = 2 * scipy.special.sph_harm(m, l, angles[0], angles[1])
                    else:
                        y = ((1 + math.pow(-1, round((k*math.cos(angles[1]))) + m)) / math.sqrt(2)) * scipy.special.sph_harm(m, l, angles[0], angles[1])
                    vector_y.append(y)
                vector_y_k.update({k: vector_y})
            matrix_y.update({index: vector_y_k})
            print('calc_y: y{yl},{ym} OK'.format(yl=l, ym=m))
    print('calc_y: OK in {time} s'.format(time=time.time() - partial_time))
    return matrix_y


def fluctuations(matrix_y, L):
    """
    :param matrix_y:
    :return:
    """
    partial_time = time.time()
    new_matrix_y = dict()
    for l in range(2, L+1):
        for m in range(l+1):
            index = healpy.sphtfunc.Alm.getidx(L, l, m)
            vector_y = []
            for k in sorted(matrix_y[index].keys()):
                y = matrix_y[index][k]
                y_star = numpy.conjugate(y)
                size = len(matrix_y[index][k])
                cn = numpy.vectorize(numpy.complex)
                cn1 = cn(numpy.random.normal(0.0, 1.0, size), numpy.random.normal(0.0, 1.0, size))
                cn1_star = numpy.conjugate(cn1)
                cn2 = cn(numpy.random.normal(0.0, 1.0, size), numpy.random.normal(0.0, 1.0, size))
                cn2_star = numpy.conjugate(cn2)
                cn3 = cn(numpy.random.normal(0.0, 1.0, size), numpy.random.normal(0.0, 1.0, size))
                cn3_star = numpy.conjugate(cn3)
                cn4 = cn(numpy.random.normal(0.0, 1.0, size), numpy.random.normal(0.0, 1.0, size))
                cn4_star = numpy.conjugate(cn4)
                psi1 = cn1 + math.pow(-1, l) * cn1_star + math.pow(-1, m) * (cn2 + math.pow(-1, l) * cn2_star)
                psi2 = cn3 + math.pow(-1, l) * cn3_star + math.pow(-1, m) * (cn4 + math.pow(-1, l) * cn4_star)
                y_sum = numpy.dot(y, psi1) + numpy.dot(y_star, psi2)
                vector_y.append(y_sum)
            new_matrix_y.update({index: vector_y})
    print('fluctuations: OK in {time} s'.format(time = time.time() - partial_time))
    return new_matrix_y


def fluctuations_pi_twisted(matrix_y, L):
    """
    :param matrix_y:
    :param L:
    :return:
    """
    partial_time = time.time()
    new_matrix_y = dict()
    for l in range(2, L+1):
        for m in range(l+1):
            index = healpy.sphtfunc.Alm.getidx(L, l, m)
            vector_y = []
            for k in sorted(matrix_y[index].keys()):
                y = matrix_y[index][k]
                y_star = numpy.conjugate(y)
                size = len(matrix_y[index][k])
                cn = numpy.vectorize(numpy.complex)
                cn1 = cn(numpy.random.normal(0.0, 1.0, size), numpy.random.normal(0.0, 1.0, size))
                cn1_star = numpy.conjugate(cn1)
                cn2 = cn(numpy.random.normal(0.0, 1.0, size), numpy.random.normal(0.0, 1.0, size))
                cn2_star = numpy.conjugate(cn2)
                psi1 = cn1 + math.pow(-1, l)*cn1_star
                psi2 = cn2 + math.pow(-1, l)*cn2_star
                y_sum = numpy.dot(y, psi1) + numpy.dot(y_star, psi2)
                vector_y.append(y_sum)
            new_matrix_y.update({index: vector_y})
    print('fluctuations_pi_twisted: OK in {time} s'.format(time=time.time() - partial_time))
    return new_matrix_y


def same_fluctuation(matrix_y, L):
    """
    :param matrix_y:
    :param L:
    :return:
    """
    partial_time = time.time()
    cn01 = numpy.complex(0.30606296, -0.56927638)
    cn02 = numpy.complex(-1.63311657, 0.04831234)
    cn03 = numpy.complex(-0.26438616, 0.13311217)
    cn04 = numpy.complex(0.59332552, 2.5154771)
    cn05 = numpy.complex(0.85092897, 0.65772095)
    cn06 = numpy.complex(-1.35946998, -1.84545289)
    cn07 = numpy.complex(-1.05840766, 1.99635083)
    cn08 = numpy.complex(-0.66826781, -1.18950684)
    cn09 = numpy.complex(0.32186199, -0.31077011)
    cn10 = numpy.complex(-0.09322461, 0.39734598)
    cn11 = numpy.complex(-0.02240994, -0.7125173)
    cn12 = numpy.complex(-0.32809129, 0.71789816)
    cn13 = numpy.complex(-2.30235563, -1.24859356)
    cn14 = numpy.complex(0.14949876, -0.84828176)
    cn15 = numpy.complex(0.13850881, -0.71363267)
    cn16 = numpy.complex(-1.92133390, 2.17336027)
    cn17 = numpy.complex(0.30844113, 0.50620514)
    cn18 = numpy.complex(-0.73566574, -1.44575845)
    cn19 = numpy.complex(-0.08296753, 1.93513265)
    cn20 = numpy.complex(-0.19886504, -0.52889081)
    cn21 = numpy.complex(0.7004854, -0.76815833)
    cn22 = numpy.complex(1.14066027, -0.56172062)
    cn23 = numpy.complex(0.54962313, -1.28936300)
    cn24 = numpy.complex(1.55558483, -2.21653618)
    cn01_star = numpy.conjugate(cn01)
    cn02_star = numpy.conjugate(cn02)
    cn03_star = numpy.conjugate(cn03)
    cn04_star = numpy.conjugate(cn04)
    cn05_star = numpy.conjugate(cn05)
    cn06_star = numpy.conjugate(cn06)
    cn07_star = numpy.conjugate(cn07)
    cn08_star = numpy.conjugate(cn08)
    cn09_star = numpy.conjugate(cn09)
    cn10_star = numpy.conjugate(cn10)
    cn11_star = numpy.conjugate(cn11)
    cn12_star = numpy.conjugate(cn12)
    cn13_star = numpy.conjugate(cn13)
    cn14_star = numpy.conjugate(cn14)
    cn15_star = numpy.conjugate(cn15)
    cn16_star = numpy.conjugate(cn16)
    cn17_star = numpy.conjugate(cn17)
    cn18_star = numpy.conjugate(cn18)
    cn19_star = numpy.conjugate(cn19)
    cn20_star = numpy.conjugate(cn20)
    cn21_star = numpy.conjugate(cn21)
    cn22_star = numpy.conjugate(cn22)
    cn23_star = numpy.conjugate(cn23)
    cn24_star = numpy.conjugate(cn24)
    new_matrix_y = dict()
    for l in range(2, L+1):
        for m in range(l+1):
            index = healpy.sphtfunc.Alm.getidx(L, l, m)
            vector_y = []
            for k in sorted(matrix_y[index].keys()):
                y_sum = numpy.complex()
                for y in matrix_y[index][k]:
                    y_star = numpy.conjugate(y)
                    if (int(k) % 2) == 0:
                        if index % 3 == 0:
                            psi1 = cn01 + math.pow(-1, l)*cn01_star + math.pow(-1, m)*(cn02 + math.pow(-1, l)*cn02_star)
                            psi2 = cn03 + math.pow(-1, l)*cn03_star + math.pow(-1, m)*(cn04 + math.pow(-1, l)*cn04_star)
                        elif index % 3 == 1:
                            psi1 = cn05 + math.pow(-1, l)*cn05_star + math.pow(-1, m)*(cn06 + math.pow(-1, l)*cn06_star)
                            psi2 = cn07 + math.pow(-1, l)*cn07_star + math.pow(-1, m)*(cn08 + math.pow(-1, l)*cn08_star)
                        else:
                            psi1 = cn09 + math.pow(-1, l)*cn09_star + math.pow(-1, m)*(cn10 + math.pow(-1, l)*cn10_star)
                            psi2 = cn11 + math.pow(-1, l)*cn11_star + math.pow(-1, m)*(cn12 + math.pow(-1, l)*cn12_star)
                    else:
                        if index % 3 == 0:
                            psi1 = cn13 + math.pow(-1, l)*cn13_star + math.pow(-1, m)*(cn14 + math.pow(-1, l)*cn14_star)
                            psi2 = cn15 + math.pow(-1, l)*cn15_star + math.pow(-1, m)*(cn16 + math.pow(-1, l)*cn16_star)
                        elif index % 3 == 1:
                            psi1 = cn17 + math.pow(-1, l)*cn17_star + math.pow(-1, m)*(cn18 + math.pow(-1, l)*cn18_star)
                            psi2 = cn19 + math.pow(-1, l)*cn19_star + math.pow(-1, m)*(cn20 + math.pow(-1, l)*cn20_star)
                        else:
                            psi1 = cn21 + math.pow(-1, l)*cn21_star + math.pow(-1, m)*(cn22 + math.pow(-1, l)*cn22_star)
                            psi2 = cn23 + math.pow(-1, l)*cn23_star + math.pow(-1, m)*(cn24 + math.pow(-1, l)*cn24_star)
                    new_y = psi1*y + psi2*y_star
                    y_sum += new_y
                vector_y.append(y_sum)
            new_matrix_y.update({index: vector_y})
    print('fluctuations: OK in {time} s'.format(time = time.time() - partial_time))
    return new_matrix_y


def save_matrix_y(filename, vector_y):
    """
    :param filename:
    :param vector_y:
    :return:
    """
    partial_time = time.time()
    with open(filename, 'wb') as f:
        pickle.dump(vector_y, f, protocol=-1)
    print('save_vector_y: OK in {time} s'.format(time = time.time() - partial_time))


def load_matrix_y(filename):
    """
    :param filename:
    :return:
    """
    partial_time = time.time()
    with open(filename, 'rb') as f:
        vector_y = pickle.load(f)
    print('load_vector_y: OK in {time} s'.format(time = time.time() - partial_time))
    return vector_y
