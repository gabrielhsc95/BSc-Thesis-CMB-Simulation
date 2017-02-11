"""
Code to find {n} and {n}n

Gabriel Cardoso 16-03-16
"""

import math
import pickle
import time

# no zeros on z

def calc_n(num_max):
    """
    Create a dictionary, where the keys are the modulus of the vectors, and the values are all the vector with the same moduli.
    The vectors components run from -numMax to numMax in x and y, but from 0 to numMax in z(properties of gravity)
    There is not the pair {0: [[0,0,0]]}, because it is trivial and useless
    :param num_max: the maximum number of the components for the vector
    :return: the dictionary
    """
    partial_time = time.time()
    n = create_keys(num_max + 1)
    n = unique_n(num_max + 1, n)
    n = permutations(n)
    n = vectors_to_angles(n)
    print('calc_n: OK in {time} s'.format(time = time.time() - partial_time))
    return n


def create_keys(num):
    """
    :param num:
    :return:
    """
    partial_time = time.time()
    n = {0.0: []}
    for nx in range(num):
        for ny in range(nx, num):
            for nz in range(ny, num):
                modulo = math.sqrt(nx**2 + ny**2 + nz**2)
                if not(modulo in n.keys()):
                    n.update({modulo: []})
    n.pop(0.0)  # remove the pair {0.0 : []}
    print('create_keys: OK in {time} s'.format(time = time.time() - partial_time))
    return n

def unique_n(num, n):
    """
    Create a dictionary, where the keys are the modulus of the vectors, and the values are vectors with the same moduli.
    There is not any negative component
    there is not any symmetry of permutation, e.g.,  there is just [1,0,0], it is missing [0,1,0] and [0,0,1]
    :param num: the maximum number of the components for the vector
    :return: the dictionary
    """
    partial_time = time.time()
    correction = int(num*1.75)  # add vector with the same moduli, but components greater than num
    for nx in range(correction):
        for ny in range(nx, correction):
            for nz in range(ny, correction):
                modulo = math.sqrt(nx**2 + ny**2 + nz**2)
                if modulo in n.keys():
                    if nz != 0:  # remove nz = 0, it also needs correction in the permutation
                        n[modulo].append([nx, ny, nz])
    print('unique_n: OK in {time} s'.format(time = time.time() - partial_time))
    return n


def permutations(n_dict):
    """
    Manipulate a dictionary to create the symmetry of permutations in the vectors, e.g., it adds [0,1,0] and
    [0,0,1] to [1,0,0]
    :param n_dict: the dictionary without symmetry of permutations
    :return: the dictionary with permutations
    """
    partial_time = time.time()
    for same_moduli in n_dict.values():
        permutations_list = []
        for vector in same_moduli:
            if not(vector[0] == vector[1] and vector[1] == vector[2]):  # when nx=ny=nz, none vector need to be added
                # (x, x, z) -> (x, z, x), (z, x, x)
                if vector[0] == vector[1] and vector[1] != vector[2] and vector[0] != 0:
                    # last condition is to avoid nz = 0
                    permutations_list.append([vector[0], vector[2], vector[1]])
                    permutations_list.append([vector[2], vector[0], vector[1]])
                # (x, y, y) - > (y, x, y), (y, y, x)
                if vector[1] == vector[2] and vector[2] != vector[0]:
                    permutations_list.append([vector[1], vector[0], vector[2]])
                    if vector[0] != 0:  # avoid nz = 0
                        permutations_list.append([vector[1], vector[2], vector[0]])
                # (x, y, x) -> (y, x, x), (x, x, y)
                if vector[0] == vector[2] and vector[2] != vector[1]:
                    permutations_list.append([vector[1], vector[0], vector[2]])
                    if vector[1] != 0:  # avoid nz = 0
                        permutations_list.append([vector[0], vector[2], vector[1]])
                # (x, y, z) -> (x, z, y), (z, x, y), (y, x, z), (y, z, x), (z, y, x)
                if vector[0] != vector[1] and vector[0] != vector[2] and vector[1] != vector[2]:
                    if vector[1] != 0:  # avoid nz = 0
                        permutations_list.append([vector[0], vector[2], vector[1]])
                        permutations_list.append([vector[2], vector[0], vector[1]])
                    permutations_list.append([vector[1], vector[0], vector[2]])
                    if vector[0] != 0:  # avoid nz = 0
                        permutations_list.append([vector[1], vector[2], vector[0]])
                        permutations_list.append([vector[2], vector[1], vector[0]])
        modulo = math.sqrt(same_moduli[0][0]**2 + same_moduli[0][1]**2 + same_moduli[0][2]**2)
        n_dict[modulo] += permutations_list
    print('permutations: OK in {time} s'.format(time = time.time() - partial_time))
    return n_dict


def negative_vectors(n_dict):
    """
    Manipulate a dictionary to add the negative components of the vectors
    There is not negative component in z due to gravity properties
    :param n_dict: the dictionary without negative components
    :return:the dictionary with negative components
    """
    # partial_time = time.time()
    for same_moduli in n_dict.values():
        negative_vectors_list = []
        for vector in same_moduli:
            if not(vector[0] == 0 and vector[1] == 0):  # when nx = 0 and ny=0, none vector need to be added
                if vector[0] != 0 and vector[1] != 0:
                    negative_vectors_list.append([vector[0], -1*vector[1], vector[2]])
                    negative_vectors_list.append([-1*vector[0], vector[1], vector[2]])
                negative_vectors_list.append([-1*vector[0], -1*vector[1], vector[2]])
        modulo = math.sqrt(same_moduli[0][0]**2+same_moduli[0][1]**2+same_moduli[0][2]**2)
        n_dict[modulo] += negative_vectors_list
    # print('negative_vectors: OK in {time} s'.format(time = time.time() - partial_time))
    return n_dict

def vectors_to_angles(n_dict):
    """
    :param n_dict:
    :return:
    """
    partial_time = time.time()
    new_n = dict()
    for k in n_dict.keys():
        new_vector = []
        for vector in n_dict[k]:
            # calculate the azimuthal angle
            if vector[0] != 0:
                phi = math.atan(float(vector[1])/float(vector[0]))
            else:
                phi = math.pi/2  # tan^-1(inf) = pi/2
                if vector[1] == 0:  # tan^-1(0/0) = 0.0
                    phi = 0.0
            # calculate the polar angle
            if vector[2] != 0:
                # sqrt returns float, float/int = float
                theta = math.atan((math.sqrt(vector[0]**2 + vector[1]**2))/vector[2])
            else:
                theta = math.pi/2  # tan^-1(inf) = pi/2
            new_vector.append([phi, theta])
        new_n.update({k: new_vector})
    print('vectors_to_angles: OK in {time} s'.format(time = time.time() - partial_time))
    return new_n


def save_n(filename, n_dict):
    """
    Save the dictionary in a pickle format, it is binary
    :param filename: the name of the file, where it will be saved
    :param n_dict: the dictionary that will be saved
    """
    partial_time = time.time()
    with open(filename, 'wb') as f:
        pickle.dump(n_dict, f, protocol=-1)
    print('save_n: OK in {time} s'.format(time = time.time() - partial_time))


def load_n(filename):
    """
    Load a dictionary saved in pickle
    :param filename: the name of the file to be loaded
    :return: the dictionary loaded
    """
    partial_time = time.time()
    with open(filename, 'rb') as f:
        n = pickle.load(f)
    print('load_n: OK in {time} s'.format(time = time.time() - partial_time))
    return n
