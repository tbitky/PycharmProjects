import numpy as np
import sympy
from sympy import re as real
from decimal import Decimal, ROUND_HALF_UP
from . import nitride_semiconductor_data as data
import xrayutilities as xu
import itertools
import pandas

properties = data.properties
wavelength = xu.wavelength('CuKa1')
Q_scan_rlu_value = 0.5


def composition_calculate(material_1, material_2, a_measured, c_measured, means='poisson'):
    x = sympy.Symbol('x')

    def vegard_equation(index_1, index_2, variable_x=x):
        acv = ['a', 'c', 'v']
        yy = [0] * 3
        for n, jj in enumerate(acv):
            yy[n] = variable_x * properties.at[index_1, jj] + (1 - variable_x) * properties.at[index_2, jj]
        return yy

    def elastic_equation(index_1, index_2, variable_x=x):
        acCC = ['a', 'c', 'C13', 'C33']
        yy = [0] * 4
        for n, jj in enumerate(acCC):
            yy[n] = variable_x * properties.at[index_1, jj] + (1 - variable_x) * properties.at[index_2, jj]
        return yy

    if means == 'poisson':
        fa, fc, fv = vegard_equation(material_1, material_2, x)
        equation = sympy.expand(fv * (c_measured - fc) * fa + (a_measured - fa) * fc)
    elif means == 'elastic':
        fa, fc, fC13, fC33 = elastic_equation(material_1, material_2, x)
        equation = sympy.expand(c_measured - fc * (1 - 2 * fc * (fC13 / fC33) * (a_measured - fa) / fa))

    solve = np.array(sympy.solve(equation, x))
    real_solves = list(map(real, solve))

    solutions = [i for i in real_solves if 1 >= i >= 0]

    if solutions:
        solution = solutions[0] * 100
        lattice_constant_a = solution * properties.at[material_1, 'a'] + (1 - solution) * properties.at[material_2, 'a']
        lattice_constant_c = solution * properties.at[material_1, 'c'] + (1 - solution) * properties.at[material_2, 'c']
        delta_a = (a_measured - lattice_constant_a) / lattice_constant_a * 100
        delta_c = (c_measured - lattice_constant_c) / lattice_constant_c * 100
        if delta_a * delta_c > 0:
            return False
        else:
            return solution, lattice_constant_a, lattice_constant_c, delta_a, delta_c
    else:
        return False


def ternary_a_c_r_calculate(qx, qy, miller_h, miller_k, miller_l, a=0.0, c=0.0, xray=wavelength,
                            magnitude=Q_scan_rlu_value):
    qx = qx / magnitude
    qy = qy / magnitude
    if qx != 0 and a == 0:
        a = abs(np.sqrt((miller_h ** 2 + miller_h * miller_k + miller_k ** 2) * 4 / 3) * xray / qx)
    if qy != 0 and c == 0:
        c = abs(miller_l * xray / qy)
    print("測定値　　　      a={0:.3f}Å, c={1:.3f}Å".format(a, c))
    searched_flag = False
    for i, j in itertools.combinations(properties.index.values, 2):
        material = j[:-1] + i
        result = composition_calculate(j, i, a, c)
        if result:
            line = "{0} {1}={2[0]:4.1f}% ,a={2[1]:.3f}Å, c={2[2]:.3f}Å, ⊿a={2[3]:+.2f}%, ⊿c={2[4]:+.2f}%".format(
                material, j, result)
            print(line)
            searched_flag = True
    if not searched_flag:
        print("There isn't such material")

    return a, c


def omega_and_ttheta_calculate(qx, qy):
    omega_ttheta_candidates = np.empty((4,2))
    A = np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)
    B = (qx ** 2 + qx + qy ** 2)

    omega_ttheta_candidates[0][0] = 2 * np.arctan((qy - A) / B)
    omega_ttheta_candidates[1][0] = omega_ttheta_candidates[0][0]
    omega_ttheta_candidates[2][0] = 2 * np.arctan((qy + A) / B)
    omega_ttheta_candidates[3][0] = omega_ttheta_candidates[2][0]

    omega_ttheta_candidates[0][1] = np.arcsin(
        2 * qy - np.sin(2 * np.arctan((qy - A) / B))) + omega_ttheta_candidates[0][0]

    omega_ttheta_candidates[1][1] = - np.arcsin(
        2 * qy - np.sin(2 * np.arctan((qy - A) / B))) + omega_ttheta_candidates[0][1] - np.pi

    omega_ttheta_candidates[2][1] = np.arcsin(
        2 * qy - np.sin(2 * np.arctan((qy + A) / B))) + omega_ttheta_candidates[2][0]

    omega_ttheta_candidates[3][1] = - np.arcsin(
        2 * qy - np.sin(2 * np.arctan((qy + A) / B))) + omega_ttheta_candidates[3][0] - np.pi

    omega_ttheta_candidates *= 180 / np.pi

    sort = range((omega_ttheta_candidates.shape[0]))
    if np.abs(omega_ttheta_candidates[0][0] - 90) >= np.abs(omega_ttheta_candidates[0][2] - 90):
        sort = [2, 3, 0, 1]

    for i in range(int(len(omega_ttheta_candidates) / 2)):
        if np.abs(omega_ttheta_candidates[2 * i][1] - 180) >= np.abs(omega_ttheta_candidates[2 * i + 1][1] - 180):
            temp = sort
            sort[2 * i + 1] = temp[2 + i]
            sort[2 * i] = temp[2 * i + 1]
    omega_ttheta_candidates = omega_ttheta_candidates[sort, sort]

    ttheta_upper_limit = 159.5285
    ttheta_lower_limit = -13.5795
    omega_upper_limit = 117.5285
    omega_lower_limit = -8.2625
    omega = []
    ttheta = []

    for i in range(len(sort)):
        if omega_lower_limit <= omega_ttheta_candidates[i][0] <= omega_upper_limit and \
                ttheta_lower_limit <= omega_ttheta_candidates[i][1] <= ttheta_upper_limit:
            omega.append(omega_ttheta_candidates[i][0])
            ttheta.append(omega_ttheta_candidates[i][1])
    line = ''
    for i in range(len(omega)):
        line += 'omega = {0[i]:.6}[deg] 2theta= {1[i]:.6}[deg]\n'.format(omega, ttheta)
    print(line)
    return omega, ttheta


def composition_and_relaxation_or_strained_lattice_constant_and_hkl_to_qxqy(material_1, material_2, composition,
                                                                             relaxation,
                                                                             lattice_constant_a, miller_h, miller_k,
                                                                             miller_l, xray=wavelength,
                                                                             magnitude=Q_scan_rlu_value):
    alloy_a = composition * properties.at[material_1, 'a'] + (1 - composition) * properties.at[material_2, 'a']
    alloy_c = composition * properties.at[material_1, 'c'] + (1 - composition) * properties.at[material_2, 'c']
    alloy_v = composition * properties.at[material_1, 'v'] + (1 - composition) * properties.at[material_2, 'v']

    try:
        if lattice_constant_a:
            real_a = lattice_constant_a
        elif relaxation:
            real_a = alloy_a * relaxation / 100
    except ValueError:
        print('input relaxation or lattice_constant_a')

    real_c = alloy_c * (1 - (real_a - alloy_a) / alloy_a / alloy_v)
    if miller_h:
        qx = magnitude * miller_h / abs(miller_h) * np.sqrt(
            (miller_h ** 2 + miller_h * miller_k + miller_k ** 2) * 4 / 3) * xray / real_a

    else:
        qx = 0.0
    if miller_l:
        qy = magnitude * miller_l * xray / real_c
    else:
        qy = 0.0
    line = 'qx={0:.5}[nm^-1]\nqy={1:.5}[nm^-1] '.format(qx, qy)
    print(line)
    return qx, qy


def fine_round(x, y=0):
    round_x = Decimal(str(x)).quantize(Decimal(str(y)), rounding=ROUND_HALF_UP)
    int_x = int(round_x)
    return int_x
