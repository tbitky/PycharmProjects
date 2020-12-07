import numpy as np
import sympy
from sympy import re as real
import pandas as pd

"""
GaN |a,c,ν,Psp,e31,e33,C13,C33
AlN |
InN |
"""
index = ['GaN', 'AlN', 'InN']
columns = ['a', 'c', 'v', 'Psp', 'e31', 'e33', 'C13', 'C33']
properties = [[3.1893, 5.1851, 0.203, -0.029, -0.49, 0.73, 100, 392],
              [3.1130, 4.9816, 0.225, -0.081, -0.60, 1.46, 127, 382],
              [3.5380, 5.7020, 0.291, -0.032, -0.57, 0.97, 94, 200]]
properties_pd = pd.DataFrame(properties, index=index, columns=columns)


def equation_calculate(a, b, a_measured, c_measured):
    x = sympy.Symbol('x')
    ax = x * properties[a][0] + (1 - x) * properties[b][0]
    cx = x * properties[a][1] + (1 - x) * properties[b][1]
    vx = x * properties[a][2] + (1 - x) * properties[b][2]
    equation = sympy.expand(vx * (c_measured - cx) * ax + (a_measured - ax) * cx)
    solve = np.array(sympy.solve(equation, x))
    real_solves = list(map(real, solve))

    solutions = [i for i in real_solves if 1 >= i >= 0]
    if solutions:
        solution = solutions[0]
        lattice_constant_a = solution * properties[a][0] + (1 - solution) * properties[b][0]
        lattice_constant_c = solution * properties[a][1] + (1 - solution) * properties[b][1]
        delta_a = float((a_measured - lattice_constant_a) / lattice_constant_a * 100)
        delta_c = float((c_measured - lattice_constant_c) / lattice_constant_c * 100)
        if delta_a * delta_c > 0:
            return False
        else:
            return solution, lattice_constant_a, lattice_constant_c, delta_a, delta_c
    else:
        return False,


def ternary_a_c_r_calculate(qx, qy, miller_h, miller_k, miller_l, xray=1.54 * 10 ** -10, a=0.0, c=0.0):
    if qx != 0 and a == 0:
        a = abs(np.sqrt((miller_h ** 2 + miller_h * miller_k + miller_k ** 2) * 4 / 3) * (xray / 2 * 10 ** 10) / qx)

    if qy != 0 and c == 0:
        c = abs(miller_l * (xray / 2 * 10 ** 10) / qy)

    print("測定値　　　      a={0:.3f}Å, c={1:.3f}Å".format(a, c))
    algan_solution = equation_calculate(1, 0, a, c)
    ingan_solution = equation_calculate(2, 0, a, c)
    inaln_solution = equation_calculate(2, 1, a, c)
    if algan_solution[0]:
        line = "AlGaN Al={1:4.1f}% ,a={0[1]:.3f}Å, c={0[2]:.3f}Å, ⊿a={0[3]:+.2f}%, ⊿c={0[4]:+.2f}%".format(
            algan_solution, algan_solution[0] * 100)
        print(line)
    if ingan_solution[0]:
        line = "InGaN In={1:4.1f}% ,a={0[1]:.3f}Å, c={0[2]:.3f}Å, ⊿a={0[3]:+.2f}%, ⊿c={0[4]:+.2f}%".format(
            ingan_solution, ingan_solution[0] * 100)
        print(line)
    if inaln_solution[0]:
        line = "InAlN In={1:4.1f}% ,a={0[1]:.3f}Å, c={0[2]:.3f}Å, ⊿a={0[3]:+.2f}%, ⊿c={0[4]:+.2f}%".format(
            inaln_solution, inaln_solution[0] * 100)
        print(line)
    if algan_solution[0] and ingan_solution[0] and inaln_solution[0]:
        print("There isn't such matter")
    return a, c


def omega_and_ttheta_calculate(qx, qy, xray=1.54 * 10 ** -10):
    omega_candidate_1 = 360 * np.arctan(
        (qy - np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)) / np.pi
    ttheta_candidate_1 = 180 * np.arcsin(2 * qy - np.sin(2 * np.arctan(
        (qy - np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)))) / np.pi + 360 * np.arctan(
        (qy - np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)) / np.pi
    ttheta_candidate_2 = -180 * np.arcsin(2 * qy - np.sin(2 * np.arctan(
        (qy - np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)))) / np.pi + 360 * np.arctan(
        (qy - np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)) / np.pi - 180

    omega_candidate_2 = 360 * np.arctan(
        (qy + np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)) / np.pi
    ttheta_candidate_3 = 180 * np.arcsin(2 * qy - np.sin(2 * np.arctan(
        (qy + np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)))) / np.pi + 360 * np.arctan(
        (qy + np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)) / np.pi
    ttheta_candidate_4 = -180 * np.arcsin(2 * qy - np.sin(2 * np.arctan(
        (qy + np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)))) / np.pi + 360 * np.arctan(
        (qy + np.sqrt(-qx ** 4 - 2 * qx ** 2 * qy ** 2 + qx ** 2 - qy ** 4 + qy ** 2)) / (
                qx ** 2 + qx + qy ** 2)) / np.pi - 180

    if 0 <= omega_candidate_1 <= 180:
        if 0<=ttheta_candidate_1 <= 180:
            omega = omega_candidate_1
            ttheta = ttheta_candidate_1
        elif 0<=ttheta_candidate_2 <= 180:
            omega = omega_candidate_1
            ttheta = ttheta_candidate_2
    elif 0<=omega_candidate_2 <= 180:
        if 0<=ttheta_candidate_3 <= 180:
            omega = omega_candidate_2
            ttheta = ttheta_candidate_3
        elif 0<=ttheta_candidate_4 <= 180:
            omega = omega_candidate_2
            ttheta = ttheta_candidate_4
    else:
        omega = None
        ttheta = None
    line = 'omega = {0:.6}[deg]\n2theta= {1:.6}[deg]'.format(omega, ttheta)
    print(line)
    return omega, ttheta


def composition_and_relaxationo_or_strained_lattice_constant_and_hkl_to_qxqy(material_1, material_2, composition,
                                                                             relaxation,
                                                                             lattice_constant_a, miller_h, miller_k,
                                                                             miller_l, xray=1.54 * 10 ** -10):
    alloy_a = composition * properties_pd.at[material_1, 'a'] + (1 - composition) * properties_pd.at[material_2, 'a']
    alloy_c = composition * properties_pd.at[material_1, 'c'] + (1 - composition) * properties_pd.at[material_2, 'c']
    alloy_v = composition * properties_pd.at[material_1, 'v'] + (1 - composition) * properties_pd.at[material_2, 'v']

    try:
        if lattice_constant_a:
            real_a = lattice_constant_a
        elif relaxation:
            real_a = alloy_a * relaxation / 100
    except ValueError:
        print('input relaxation or lattice_constant_a')
    real_c = alloy_c * (1 - (real_a - alloy_a) / alloy_a / alloy_v)
    if miller_h:
        qx = miller_h / abs(miller_h) * np.sqrt((miller_h ** 2 + miller_h * miller_k + miller_k ** 2) * 4 / 3) * (
                xray / 2 * 10 ** 10) / real_a
    else:
        qx = 0.0
    if miller_l:
        qy = miller_l * (xray / 2 * 10 ** 10) / real_c
    else:
        qy = 0.0
    line = 'qx={0:.5}[nm^-1]\nqy={1:.5}[nm^-1] '.format(qx, qy)
    print(line)
    return qx, qy


def main():
    """
    input
    """
    try:
        qx, qy = 0, 0
        qx, qy = map(float, input('qx[rlu] qy[rlu]入力:').split())
        # qx=qx*1.54/4/np.pi
        # qy=qy*1.54/4/np.pi
    except:
        pass
    try:
        hh, kk, ll = 0, 0, 0
        hh, kk, ll = map(int, input('h k l入力:').split())
    except:
        pass
    try:
        material_1, material_2, composition = '', '', ''
        material_1, material_2, composition = map(str, input('半導体1 半導体2 組成入力：').split())
        composition = float(composition) / 100
    except:
        pass
    try:
        lattice_constant_a = 0
        lattice_constant_c = 0
        lattice_constant_a, lattice_constant_c = map(float, input('a軸c軸格子定数[Å]入力：').split())
    except:
        pass
    try:
        relaxation = 0
        relaxation = float(input('緩和率[%]入力：'))
    except:
        pass
    try:
        omega, ttheta = 0, 0
        omega, ttheta = map(float, input('omega[deg] ttheta[deg]入力:').split())
    except:
        pass
    print('\n計算結果' + '-' * 20)
    if all((material_1, material_2, composition)) and any((hh, kk, ll)) and \
            any((lattice_constant_a, relaxation)) and not any((qx, qy)):
        qx, qy = composition_and_relaxationo_or_strained_lattice_constant_and_hkl_to_qxqy(material_1, material_2,
                                                                                          composition,
                                                                                          relaxation,
                                                                                          lattice_constant_a, hh, kk,
                                                                                          ll)
    elif all((omega, ttheta)) and not any((qx, qy)):
        qx = (np.cos(omega * np.pi / 180) - np.cos((ttheta - omega) * np.pi / 180)) / 2
        qy = (np.sin(omega * np.pi / 180) + np.sin((ttheta - omega) * np.pi / 180)) / 2
        line = 'qx={0:.5}[nm^-1]\nqy={1:.5}[nm^-1] '.format(qx, qy)
        print(line)
    if any((hh, kk, ll)) and not all(
            (material_1, material_2, composition, lattice_constant_a, lattice_constant_c, relaxation)):
        ternary_a_c_r_calculate(qx, qy, hh, kk, ll, a=lattice_constant_a, c=lattice_constant_c)
    if any((hh, kk, ll)) and not any((omega, ttheta)):
        omega_and_ttheta_calculate(qx, qy)


if __name__ == '__main__':
    main()
