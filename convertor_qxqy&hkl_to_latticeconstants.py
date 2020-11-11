import numpy as np
import sympy
from sympy import re as real

"""
GaN |a,c,ν,Psp,e31,e33,C13,C33
AlN |
InN |
"""
properties = [[3.1893, 5.1851, 0.203, -0.029, -0.49, 0.73, 100, 392],
              [3.1130, 4.9816, 0.225, -0.081, -0.60, 1.46, 127, 382],
              [3.5380, 5.7020, 0.291, -0.032, -0.57, 0.97, 94, 200]]


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
        delta_a = (a_measured - lattice_constant_a) / lattice_constant_a * 100
        delta_c = (c_measured - lattice_constant_c) / lattice_constant_c * 100
        if delta_a * delta_c > 0:
            return False
        else:
            return solution, lattice_constant_a, lattice_constant_c, delta_a, delta_c
    else:
        return False,


def ternary_a_c_r_calculate(qx, qy, miller_h, miller_k, miller_l, xray=1.54 * 10 ** -10):
    a = abs(np.sqrt((miller_h ** 2 + miller_h * miller_k + miller_k ** 2) * 4 / 3) * (xray / 2 * 10 ** 10) / qx)
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
    omega = sympy.Symbol('omega', positive=True)
    ttheta = sympy.Symbol('ttheta', positive=True)
    fx = qx - (sympy.cos(omega * sympy.pi / 180) - sympy.cos((ttheta - omega) * sympy.pi / 180)) / 2
    fy = qy - (sympy.sin(omega * sympy.pi / 180) + sympy.sin((ttheta - omega) * sympy.pi / 180)) / 2
    solves = sympy.solve([fx, fy], [omega, ttheta])

    for i, j in enumerate(solves):
        calculated_omega = None
        calculated_ttheta = None
        for k, l in enumerate(j):
            if l < 0:
                break
        else:
            calculated_omega = solves[i][0]
            calculated_ttheta = solves[i][1]
    line = 'omega = {0:.6}[deg]\n2theta= {1:.6}[deg]'.format(calculated_omega, calculated_ttheta)
    print(line)
    return calculated_omega, calculated_ttheta


def main():
    qx, qy = map(float, input('qx[rlu] qy入力[rlu]:').split())
    hh, kk, ll = map(int, input('h k l入力:').split())
    ternary_a_c_r_calculate(qx, qy, hh, kk, ll)
    omega_and_ttheta_calculate(qx, qy)


if __name__ == '__main__':
    main()
