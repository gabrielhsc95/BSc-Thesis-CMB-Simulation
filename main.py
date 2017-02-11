"""
The main file to simulate

Gabriel Cardoso 16-03-16
"""

import time
import calcn
import calcjp
import calcy
import numpy
import csv
import healpy
import os
import math
from multiprocessing import Pool
from matplotlib import pyplot


start_time = time.time()


def product_jpy(matrix_jp, matrix_y, L):
    """
    Multiply vector x and vector y
    :param matrix_jp: the vector x
    :param matrix_y:  the vector y
    :return: A matrix with all the a's
    """
    partial_time = time.time()
    jpy = dict()
    for l in range(2, L+1):
        for m in range(l+1):
            index = healpy.sphtfunc.Alm.getidx(L, l, m)
            j_l = matrix_jp[l-2]
            y_f = matrix_y[index]
            alm = numpy.dot(y_f, j_l)
            jpy.update({index: alm})
    print('product_xy: OK in {time} s'.format(time = time.time() - partial_time))
    return jpy


def save_alm_fits(filename, jpy):
    """
    :param filename:
    :param jpy:
    :return:
    """
    size = max(jpy.keys()) + 1
    alms = [complex(0.0, 0.0)] * size
    for index in jpy.keys():
        alms[index] = jpy[index]
    alms = numpy.array(alms)
    healpy.fitsfunc.write_alm(filename, alms)
    return alms


def load_alm_fits(filename):
    """
    :param filename:
    :return:
    """
    alms = healpy.fitsfunc.read_alm(filename)
    return alms


def save_alms(filename, jpy, L):
    """
    :param filename:
    :param jpy:
    :param L:
    :return:
    """
    with open(filename, 'w') as csvfile:
        fieldnames = ['real', 'imaginary']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for l in range(2, L+1):
            for m in range(l+1):
                index = healpy.sphtfunc.Alm.getidx(L, l, m)
                writer.writerow({'real': jpy[index].real, 'imaginary': jpy[index].imag})


def preparation(num_max, l, delta_eta):
    """
    :param num_max:
    :param l:
    :param delta_eta:
    :return:
    """
    n = calcn.calc_n(num_max)
    calcn.save_n('n.pkl', n)
    jp = calcjp.calc_jp(l, delta_eta, n)
    calcjp.save_jp('jp.pkl', jp)
    y = calcy.calc_y(l, n)
    calcy.save_matrix_y('y.pkl', y)


def preparation_pi_twisted(num_max, l, delta_eta):
    """
    :param num_max:
    :param l:
    :param delta_eta:
    :return:
    """
    n = calcn.calc_n(num_max)
    calcn.save_n('n.pkl', n)
    jp = calcjp.calc_jp(l, delta_eta, n)
    calcjp.save_jp('jp.pkl', jp)
    y = calcy.calc_y_pi_twisted(l, n)
    calcy.save_matrix_y('y_pi_twisted.pkl', y)


def simulate(times):
    """
    :param l:
    :param times:
    :return:
    """
    jp = calcjp.load_jp('jp.pkl')
    y = calcy.load_matrix_y('y.pkl')
    l = len(jp) + 1
    partial_time = time.time()
    for i in range(1, times+1):
        y_f = calcy.fluctuations(y, l)
        jpy = product_jpy(jp, y_f, l)
        num_save = time.time()
        alms = save_alm_fits('simulations/alms{run}.fits'.format(run=num_save), jpy)
        map = healpy.sphtfunc.alm2map(alms, 128)
        healpy.mollview(map)
        pyplot.savefig('simulations/map{run}.png'.format(run=num_save))
    print('simulate: OK in {time} s'.format(time = time.time() - partial_time))


def simulate_pi_twisted(times):
    """
    :param times:
    :return:
    """
    jp = calcjp.load_jp('jp.pkl')
    y = calcy.load_matrix_y('y_pi_twisted.pkl')
    l = len(jp) + 1
    partial_time = time.time()
    for i in range(1, times+1):
        y_f = calcy.fluctuations_pi_twisted(y, l)
        jpy = product_jpy(jp, y_f, l)
        num_save = time.time()
        alms = save_alm_fits('simulations/alms{run}.fits'.format(run=num_save), jpy)
        map = healpy.sphtfunc.alm2map(alms, 128)
        healpy.mollview(map)
        pyplot.savefig('simulations/map{run}.png'.format(run=num_save))
    print('simulate: OK in {time} s'.format(time = time.time() - partial_time))


def simulation_parallel(times):
    """
    :param times:
    :return:
    """
    processes = os.cpu_count()
    each = int(times/processes)
    remainder = times % processes
    divided_process = [each] * processes
    divided_process.append(remainder)
    with Pool() as p:
        p.map(simulate, divided_process)


def simulation_pi_twsited_parallel(times):
    """
    :param times:
    :return:
    """
    processes = os.cpu_count()
    each = int(times/processes)
    remainder = times % processes
    divided_process = [each] * processes
    divided_process.append(remainder)
    with Pool() as p:
        p.map(simulate_pi_twisted, divided_process)


def convergence(n_max, l, m):
    """
    :param n_max:
    :param l:
    :param m:
    :return:
    """
    for num_max in range(80, n_max + 1, 5):
        n = calcn.calc_n(num_max)
        jp = calcjp.calc_jp(l, 1.0, n)
        y = calcy.calc_y(l, n)
        y = calcy.same_fluctuation(y, l)
        jpy = product_jpy(jp, y, l)
        index = healpy.sphtfunc.Alm.getidx(l, l, m)
        with open('convergence_3.txt', 'a') as csvfile:
            fieldnames = ['n', 'real', 'imaginary']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
            writer.writerow({'n': num_max, 'real': jpy[index].real, 'imaginary': jpy[index].imag})
        print('convergence: n', num_max, ' OK')

# convergence(100, 15, 4)
# preparation(3, 5, 0.5)
# simulate(1)
# simulation_parallel(16)


# preparation_pi_twisted(10, 10, 1.0)
# simulate_pi_twisted(1)
# simluation_pi_twisted_parallel(8)


print('Everything: OK in {time} s'.format(time = time.time() - start_time))
