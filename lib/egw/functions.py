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


def equation_calculate(material_1, material_2, a_measured, c_measured):
    x = sympy.Symbol('x')

    def vegard_equation(index_1, index_2, variable_x=x):
        acv = ['a', 'c', 'v']
        yy = [0] * 3
        for n, jj in enumerate(acv):
            yy[n] = variable_x * properties.at[index_1, jj] + (1 - variable_x) * properties.at[index_2, jj]
        return yy

    fa, fc, fv = vegard_equation(material_1, material_2, x)
    equation = sympy.expand(fv * (c_measured - fc) * fa + (a_measured - fa) * fc)
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


def ternary_a_c_r_calculate(qx, qy, miller_h, miller_k, miller_l, xray=wavelength):
    a = abs(np.sqrt((miller_h ** 2 + miller_h * miller_k + miller_k ** 2) * 4 / 3) * xray / qx)
    c = abs(miller_l / 2 * xray / qy)
    print("測定値　　　      a={0:.3f}Å, c={1:.3f}Å".format(a, c))
    searched_flag = False
    for i, j in itertools.combinations(properties.index.values, 2):
        material = j[:-1] + i
        result = equation_calculate(j, i, a, c)
        if result:
            line = "{0} {1}={2[0]:4.1f}% ,a={2[1]:.3f}Å, c={2[2]:.3f}Å, ⊿a={2[3]:+.2f}%, ⊿c={2[4]:+.2f}%".format(
                material, j, result)
            print(line)
            searched_flag = True
    if not searched_flag:
        print("There isn't such material")

    return a, c


def fine_round(x, y=0):
    round_x = Decimal(str(x)).quantize(Decimal(str(y)), rounding=ROUND_HALF_UP)
    int_x = int(round_x)
    return int_x
