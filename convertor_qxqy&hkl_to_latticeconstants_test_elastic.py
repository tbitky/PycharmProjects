import numpy as np
import sympy
from sympy import re as real
import pandas as pd
import lib.egw.functions as egw

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
        qx, qy = egw.composition_and_relaxation_or_strained_lattice_constant_and_hkl_to_qxqy(material_1, material_2,
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
        egw.ternary_a_c_r_calculate(qx, qy, hh, kk, ll, a=lattice_constant_a, c=lattice_constant_c)
    if any((hh, kk, ll)) and not any((omega, ttheta)):
        egw.omega_and_ttheta_calculate(qx, qy)


if __name__ == '__main__':
    main()
