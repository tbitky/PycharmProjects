from numpy.linalg import solve
import tkinter
from tkinter import filedialog
from tkinter import messagebox
import os
import sympy
from sympy import I
import scipy
import re
from sympy import re as real
from scipy.signal import peak_widths
import pandas as pd
import numpy as np
import xrayutilities as xu
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import maximum_filter
from decimal import Decimal, ROUND_HALF_UP
from lmfit import models
import lmfit
import itertools
from ipywidgets import Output
from matplotlib.patches import Ellipse

"""
GaN |a,c,ν,Psp,e31,e33,C13,C33
AlN |
InN |
"""
properties = [[3.189, 5.185, 0.203, -0.029, -0.49, 0.73, 100, 392],
              [3.112, 4.982, 0.225, -0.081, -0.60, 1.46, 127, 382],
              [3.548, 5.760, 0.291, -0.032, -0.57, 0.97, 94, 200]]


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

    # 小さいピーク値を排除（最大ピーク値のorder倍以下のピークは排除）
    temp = np.ma.array(detected_peaks, mask=~(detected_peaks >= detected_peaks.max() * order))
    peaks_index = np.where((temp.mask != True))
    return peaks_index


def fine_round(x, y=0):
    round_x = Decimal(str(float(x))).quantize(Decimal(str(y)), rounding=ROUND_HALF_UP)
    int_x = int(round_x)
    return int_x


def fitting_gaussian_start_parameters(data, peaks, count_time):
    start = [0] * (len(peaks[0]) * 7)
    for i in range(len(peaks[0])):
        xpeakposi = peak_widths(data[2][peaks[0][i]], [peaks[1][i], ])
        ypeakposi = peak_widths(data[2].T[peaks[1][i]], [peaks[0][i], ])
        x_peak_width = data[0][peaks[0][i]][fine_round(xpeakposi[3])] - data[0][peaks[0][i]][
            fine_round(xpeakposi[2])]
        y_peak_width = data[1][fine_round(ypeakposi[3])][peaks[1][i]] - data[1][
            fine_round(ypeakposi[2])][peaks[1][i]]

        start[i * 7] = data[0][peaks[0][i]][peaks[1][i]]
        start[i * 7 + 1] = data[1][peaks[0][i]][peaks[1][i]]
        start[i * 7 + 2] = x_peak_width / (2 * np.sqrt(2 * np.log(2)))
        start[i * 7 + 3] = y_peak_width / (2 * np.sqrt(2 * np.log(2)))
        start[i * 7 + 4] = data[2][peaks[0][i]][peaks[1][i]] / count_time
        start[i * 7 + 5] = np.min(data[2])
        start[i * 7 + 6] = np.pi / 2
    return start


def fitting_gaussian(x, y, *start):
    func = 0
    for i in range(len(start) // 7):
        func = func + xu.math.functions.Gauss2d(x, y, *start[i * 7:(i + 1) * 7])
    return func


def fitting_voigt_start_parameters(data, peaks, count_time):
    start = [0] * (len(peaks[0]) * 8)
    for i in range(len(peaks[0])):
        xpeakposi = peak_widths(data[2][peaks[0][i]], [peaks[1][i], ])
        ypeakposi = peak_widths(data[2].T[peaks[1][i]], [peaks[0][i], ])
        x_peak_width = data[0][peaks[0][i]][fine_round(xpeakposi[3])] - data[0][peaks[0][i]][
            fine_round(xpeakposi[2])]
        y_peak_width = data[1][fine_round(ypeakposi[2])][peaks[1][i]] - data[1][
            fine_round(ypeakposi[3])][peaks[1][i]]

        start[i * 8] = data[0][peaks[0][i]][peaks[1][i]]
        start[i * 8 + 1] = data[1][peaks[0][i]][peaks[1][i]]
        start[i * 8 + 2] = x_peak_width
        start[i * 8 + 3] = y_peak_width
        start[i * 8 + 4] = data[2][peaks[0][i]][peaks[1][i]] / count_time
        start[i * 8 + 5] = np.min(data[2])
        start[i * 8 + 6] = np.pi / 2
        start[i * 8 + 7] = 0
    return start


def pseudo_voigt_2d(x, y, amplitude, xcenter, ycenter, xsigma, ysigma, xgamma, ygamma, angle, eta):
    gparams = [xcenter, ycenter, xsigma, ysigma, amplitude, 0, angle]
    lparams = [xcenter, ycenter, xgamma, ygamma, amplitude, 0, angle]
    func = eta * xu.math.functions.Lorentz2d(x, y, *lparams) + (1 - eta) * xu.math.functions.Gauss2d(x, y, *gparams)
    return func


def voigt_2d(x, y, amplitude, xcenter, ycenter, xsigma, ysigma, xgamma, ygamma):
    xz = (x - xcenter + 1j * xgamma) / (xsigma * np.sqrt(2.0))
    yz = (y - ycenter + 1j * ygamma) / (ysigma * np.sqrt(2.0))
    xw = scipy.special.wofz(xz)
    yw = scipy.special.wofz(yz)
    model = amplitude * xw.real * yw.real / (xsigma * ysigma * np.sqrt(2.0 * np.pi) ** 2)
    return model


def constant_func(x, y, c):
    return c


def lmfit_voigt_startparams(data, peaks, counttime):
    global values
    values = ['amplitude', 'xcenter', 'ycenter', 'xsigma', 'ysigma', 'xgamma', 'ygamma', 'angle', 'eta']
    startparams = pd.DataFrame(data=[[0.0, ] * len(values)], index=['model0_', ],
                               columns=values)
    costrains = pd.DataFrame(data=[
        [np.inf, 0, True],
        [np.inf, -np.inf, False],
        [np.inf, -np.inf, False],
        [np.inf, 0, True],
        [np.inf, 0, True],
        [np.inf, 0, True],
        [np.inf, 0, True],
        [np.pi, -np.pi, True],
        [1, 0, True]],
        index=values, columns=('max', 'min', 'vary'))
    for i in range(len(peaks[0])):
        amplitude = data[2][peaks[0][i]][peaks[1][i]]

        xcenter = data[0][peaks[0][i]][peaks[1][i]]
        ycenter = data[1][peaks[0][i]][peaks[1][i]]

        xpeakposi = peak_widths(data[2][peaks[0][i]], [peaks[1][i], ])
        xpeak_width = data[0][peaks[0][i]][fine_round(xpeakposi[3])] - data[0][peaks[0][i]][
            fine_round(xpeakposi[2])]
        ypeakposi = peak_widths(data[2].T[peaks[1][i]], [peaks[0][i], ])
        ypeak_width = data[1][fine_round(ypeakposi[2])][peaks[1][i]] - data[1][
            fine_round(ypeakposi[3])][peaks[1][i]]

        xsigma = xpeak_width / (2 * np.sqrt(2 * np.log(2)))
        xgamma = xpeak_width
        ysigma = ypeak_width / (2 * np.sqrt(2 * np.log(2)))
        ygamma = ypeak_width
        angle = np.pi / 2
        eta = 0

        name = 'model' + str(i) + '_'

        startparams.loc[name] = [amplitude, xcenter, ycenter, xsigma, ysigma, xgamma, ygamma, angle, eta]
    return startparams, costrains


def lmfit_voigt_fit(data, startparams, costrains):
    # startparams:XCEN,YCEN,
    model = lmfit.Model(func=constant_func, independent_vars=['x', 'y'], prefix='const_')
    for i in range(len(startparams)):
        mod = lmfit.Model(func=pseudo_voigt_2d, independent_vars=['x', 'y'], prefix='model' + str(i) + '_')
        model += mod

    params = model.make_params()
    params['const_c'].set(value=np.min(data[2]), vary=False)
    for i in range(len(startparams)):
        for j in values:
            ii = 'model' + str(i) + '_'
            params[ii + j].set(
                value=startparams.at[ii, j],
                max=costrains.at[j, 'max'],
                min=costrains.at[j, 'min'],
                vary=costrains.at[j, 'vary']
            )
    output = model.fit(data=data[2].ravel(), params=params, x=data[0].ravel(), y=data[1].ravel())
    return output


def coefficient_of_determination_r2(experiment, expection):
    mean = np.mean(experiment)
    diff = experiment - expection
    numirator = 0
    denominator = 0
    for i in diff.ravel():
        numirator += (i - mean) ** 2
    for j in experiment.ravel():
        denominator += (j - mean) ** 2
    r2 = 1 - numirator / denominator
    return r2


def main():
    root = tkinter.Tk()
    root.withdraw()
    fTyp = [("xrdファイル", "*.xrdml")]
    iDir = os.path.abspath(os.path.dirname(r"\\Altair\SR4000data\003_XRD\ "))
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
    peak_indices = detect_peaks(twod_datas[2], filter_size=10, order=0.01)
    count_time = xrdfileoption.scan.ddict['countTime'][0]
    startparams, costrains = lmfit_voigt_startparams(twod_datas, peak_indices, count_time)
    hh, kk, ll = map(int, input('hkl入力:').split())
    Result = lmfit_voigt_fit(twod_datas, startparams, costrains)
    fitdatas = Result.best_fit.reshape(twod_datas[2].shape)
    r2 = coefficient_of_determination_r2(twod_datas[2], fitdatas)
    print(str(r2))
    fitted_params_dict = Result.best_values

    new_values = ['a', 'c']
    fit_values = new_values + values
    fitted_params = pd.DataFrame(data=[[0.0] * len(fit_values)], index=['Peak0', ], columns=fit_values)

    for i in range(int((len(fitted_params_dict) - 1) / len(values))):
        for j in values:
            fitted_params.at['Peak' + str(i), j] = fitted_params_dict['model' + str(i) + '_' + j]
        lattce_constants = ternary_a_c_r_calculate(float(fitted_params.at['Peak' + str(i), 'xcenter']) / 10,
                                                   float(fitted_params.at['Peak' + str(i), 'ycenter']) / 10,
                                                   hh, kk, ll)
        for k, l in enumerate(new_values):
            fitted_params.at['Peak' + str(i), l] = lattce_constants[k]

    print(Result.fit_report())
    if len(fitted_params) % 2 == 0:
        fig_num = 2 + len(fitted_params) // 2
    elif len(fitted_params) % 2 == 1:
        fig_num = 1 + (len(fitted_params) + 1) // 2
    fig, axes = plt.subplots(2, fig_num, figsize=(10, 4.5))
    axes = axes.ravel()
    "plot measured data as axes[0]"
    cm0 = plt.get_cmap("jet")
    cNorm0 = LogNorm(vmin=np.min(twod_datas[2][nonzero_indices]) / xrdfileoption.scan.ddict['countTime'][0],
                     vmax=np.max(twod_datas[2]) / xrdfileoption.scan.ddict['countTime'][0])
    scalarMap0 = matplotlib.cm.ScalarMappable(norm=cNorm0, cmap=cm0)
    axes[0].set_title('measured')
    axes[0].scatter(twod_datas[0][nonzero_indices], twod_datas[1][nonzero_indices], s=3, marker='o',
                    c=scalarMap0.to_rgba(twod_datas[2][nonzero_indices] / xrdfileoption.scan.ddict['countTime'][0]))

    "plot calculated data as axes[1]"
    nonzero_fitindices = np.where(fitdatas > 0)
    axes[1].set_title('calculated')
    axes[1].scatter(twod_datas[0][nonzero_fitindices], twod_datas[1][nonzero_fitindices], s=3, marker='o', alpha=1,
                    c=scalarMap0.to_rgba(fitdatas[nonzero_fitindices] / xrdfileoption.scan.ddict['countTime'][0]))

    axes[1].scatter(twod_datas[0][peak_indices], twod_datas[1][peak_indices], s=5, marker='o', facecolor='w',
                    edgecolor='k')
    for i in range(len(peak_indices[0])):
        axes[1].annotate("Peak" + str(i), (
            twod_datas[0][peak_indices[0][i]][peak_indices[1][i]],
            twod_datas[1][peak_indices[0][i]][peak_indices[1][i]]),
                         color='w')
    out = Output()
    display(out)
    ellipses = pd.DataFrame(columns=('xy_cordinate', 'horizontal_width', 'vertical_width', 'angle'))
    for i in fitted_params.index.values:
        ellipses.at[i, 'xy_cordinate'] = (fitted_params.at[i, 'amplitude':]

    class DraggableEllipse:
        def __init__(self, ellipse):
            self.ellipse = ellipse

        def connect(self):
            'connect to all the events we need'
            self.cidpress = self.rect.figure.canvas.mpl_connect(
                'button_press_event', self.on_press)
            self.cidrelease = self.rect.figure.canvas.mpl_connect(
                'button_release_event', self.on_release)
            self.cidmotion = self.rect.figure.canvas.mpl_connect(
                'motion_notify_event', self.on_motion)

    "color bar"
    fig.colorbar(scalarMap0, ax=axes[0])
    fig.colorbar(scalarMap0, ax=axes[1])
    plt.subplots_adjust(right=0.85)
    plt.subplots_adjust(wspace=0.1)
    for i in range(len(axes)):
        axes[i].set_xlabel('qx(${nm^{-1}}$)')
        axes[i].set_ylabel('qy(${nm^{-1}}$)')

    "cross section"
    posi = float(input("入力して:"))
    xpoints = np.array([0, ] * len(twod_datas[0]))
    flag = False
    for i in range(len(twod_datas[0])):
        xpoints[i] = np.abs(np.asarray(twod_datas[0][i]) - posi).argmin()
        if flag:
            xpoints = np.delete(xpoints, slice(i, None))
            break
        if xpoints[i] == xpoints[i - 1]:
            flag = True
    axes[2].set_title('cross section')
    axes[2].scatter(twod_datas[1].T[xpoints], twod_datas[2].T[xpoints] / xrdfileoption.scan.ddict['countTime'][0],
                    color='k', s=5, marker='o')
    for i in fitted_params.index.values:
        params = fitted_params.loc[i, 'amplitude':]
        one_fit = pseudo_voigt_2d(twod_datas[0], twod_datas[1], *params)
        axes[2].scatter(twod_datas[1].T[xpoints], one_fit.T[xpoints] / xrdfileoption.scan.ddict['countTime'][0],
                        label=i, s=5, marker='o')
    axes[2].legend()

    for i, j in enumerate(fitted_params.index.values):
        params = fitted_params.loc[j, 'amplitude':]
        one_fit = pseudo_voigt_2d(twod_datas[0], twod_datas[1], *params)
        one_nonzero_fitindices = np.where(one_fit > 0)
        axes[i+3].set_title(j)
        axes[i + 3].scatter(twod_datas[0][one_nonzero_fitindices], twod_datas[1][one_nonzero_fitindices], s=3,
                            marker='o', alpha=1,
                            c=scalarMap0.to_rgba(
                                one_fit[one_nonzero_fitindices] / xrdfileoption.scan.ddict['countTime'][0]))

    plt.show()


if __name__ == '__main__':
    main()
