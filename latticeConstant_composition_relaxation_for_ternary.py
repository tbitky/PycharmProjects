
import tkinter
from tkinter import filedialog
import os
import sympy
from sympy import re as real
import numpy as np
import xrayutilities as xu
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import maximum_filter
from decimal import Decimal, ROUND_HALF_UP

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


def ternary_a_c_r_calculate(qx, qy, miller_h, miller_k, miller_l):
    a = abs(2 * np.sqrt((miller_h ** 2 + miller_h * miller_k + miller_k ** 2) / 3) / qx)
    c = abs(miller_l / qy)
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


def omega_and_2thete_convert_to_qx_and_qy(omega, twothete, wavelength=1.54 * 10 ** -10):
    qx = np.empty(len(omega))
    qy = np.empty(len(omega))
    for i in range(len(omega)):
        qx[i], qy[i] = (np.cos(np.radians(omega[i])) - np.cos(
            np.radians(twothete[i] - omega[i]))) / wavelength * 10 ** -9, (
                               np.sin(np.radians(omega[i])) + np.sin(
                           np.radians(twothete[i] - omega[i]))) / wavelength * 10 ** -9
    return qx, qy


def detect_peaks(image, filter_size=3, order=0.5):
    local_max = maximum_filter(image, footprint=np.ones((filter_size, filter_size)), mode='constant')
    detected_peaks = np.ma.array(image, mask=~(image == local_max))

    temp = np.ma.array(detected_peaks, mask=~(detected_peaks >= detected_peaks.max() * order))
    peaks_index = np.where((temp.mask != True))
    return peaks_index


def fine_round(x, y=0):
    round_x = Decimal(str(float(x))).quantize(Decimal(str(y)), rounding=ROUND_HALF_UP)
    int_x = int(round_x)
    return int_x


def main():
    root = tkinter.Tk()
    root.withdraw()
    fTyp = [("xrdファイル", "*.xrdml")]
    iDir = os.path.abspath(os.path.dirname(r"\\Sirius\carbon-nas\SR4000\003_XRD\RSM"))
    filepath = filedialog.askopenfilename(filetypes=fTyp, initialdir=iDir)

    xrdfileoption = xu.io.panalytical_xml.XRDMLFile(os.path.basename(filepath), path=os.path.dirname(filepath))
    rawxrdfile = xu.io.panalytical_xml.getxrdml_map(os.path.basename(filepath), path=os.path.dirname(filepath))

    detect_range = (3,) + xrdfileoption.scan.ddict['detector'].shape
    twod_datas = np.array(rawxrdfile).reshape((detect_range))
    qx, qy = omega_and_2thete_convert_to_qx_and_qy(rawxrdfile[0], rawxrdfile[1])
    twod_datas[0] = np.array(qx).reshape(xrdfileoption.scan.ddict['detector'].shape)
    twod_datas[1] = np.array(qy).reshape(xrdfileoption.scan.ddict['detector'].shape)
    twod_datas = twod_datas.transpose(0, 2, 1)
    twod_datas = twod_datas[:, ::-1, ::-1]
    nonzero_indices = np.where(twod_datas[2] > 0)
    peak_indices = detect_peaks(twod_datas[2], filter_size=10, order=0.02)
    count_time = xrdfileoption.scan.ddict['countTime'][0]

    hh, kk, ll = map(int, input('hkl入力:').split())

    for i in range((len(peak_indices[0]))):
        print("\nPeak" + str(i))
        xx = twod_datas[0][peak_indices[0][i]][peak_indices[1][i]]
        yy = twod_datas[1][peak_indices[0][i]][peak_indices[1][i]]
        ternary_a_c_r_calculate(float(xx) / 10, float(yy) / 10, hh, kk, ll)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    "plot measured data as axes[0]"
    cm0 = plt.get_cmap("jet")
    cNorm0 = LogNorm(vmin=np.min(twod_datas[2][nonzero_indices]) / xrdfileoption.scan.ddict['countTime'][0],
                     vmax=np.max(twod_datas[2]) / xrdfileoption.scan.ddict['countTime'][0])
    scalarMap0 = matplotlib.cm.ScalarMappable(norm=cNorm0, cmap=cm0)
    ax.set_title('measured')
    ax.scatter(twod_datas[0][nonzero_indices], twod_datas[1][nonzero_indices], s=3, marker='o',
               c=scalarMap0.to_rgba(twod_datas[2][nonzero_indices] / count_time))

    ax.scatter(twod_datas[0][peak_indices], twod_datas[1][peak_indices], s=5, marker='o', facecolor='w',
               edgecolor='k')
    for i in range(len(peak_indices[0])):
        ax.annotate("Peak" + str(i), (
            twod_datas[0][peak_indices[0][i]][peak_indices[1][i]],
            twod_datas[1][peak_indices[0][i]][peak_indices[1][i]]),
                    color='w')

    "color bar"
    fig.colorbar(scalarMap0, ax=ax)
    plt.subplots_adjust(right=0.85)
    plt.subplots_adjust(wspace=0.1)
    ax.set_xlabel('qx(${nm^{-1}}$)')
    ax.set_ylabel('qy(${nm^{-1}}$)')

    plt.show()


if __name__ == '__main__':
    main()
