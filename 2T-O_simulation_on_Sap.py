from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import tkinter
from tkinter import filedialog
import xrayutilities as xu
import numpy as np
import pandas as pd
import os


def func(x, *params):
    # paramsの長さでフィッティングする関数の数を判別。
    num_func = int(len(params) / 3)

    # ガウス関数にそれぞれのパラメータを挿入してy_listに追加。
    y_list = []
    for i in range(num_func):
        y = np.zeros_like(x)
        param_range = list(range(3 * i, 3 * (i + 1), 1))
        amp = params[int(param_range[0])]
        ctr = params[int(param_range[1])]
        wid = params[int(param_range[2])]
        stdrd_dist = np.exp(-((x - ctr) / wid) ** 2)
        notzero_index = np.where(stdrd_dist > 0)
        y[notzero_index] = y[notzero_index] + np.log10(amp) * stdrd_dist[notzero_index] + np.log10(params[-1])
        y_list.append(y)

    # y_listに入っているすべてのガウス関数を重ね合わせる。
    y_sum = np.zeros_like(x)
    for i in y_list:
        y_sum = y_sum + i

    # 最後にバックグラウンドを追加。
    y_sum = y_sum + np.log10(params[-1])

    return y_sum


def fit_plot(x, *params):
    num_func = int(len(params) / 3)
    y_list = []
    for i in range(num_func):
        y = np.zeros_like(x)
        param_range = list(range(3 * i, 3 * (i + 1), 1))
        amp = params[int(param_range[0])]
        ctr = params[int(param_range[1])]
        wid = params[int(param_range[2])]
        stdrd_dist = np.exp(-((x - ctr) / wid) ** 2)
        notzero_index = np.where(stdrd_dist > 0.012)
        y[notzero_index] = y[notzero_index] + np.log10(amp) * stdrd_dist[notzero_index] + np.log10(params[-1])
        y_list.append(y)
    return y_list


def main():
    """
    ファイルダイアログを開いてfilepathに格納
    """
    root = tkinter.Tk()
    root.withdraw()
    fTyp = [("xrdファイル", "*.xrdml")]
    iDir = os.path.abspath(os.path.dirname(r"\\Sirius\carbon-nas\SR4000\003_XRD"))
    filepath = filedialog.askopenfilename(filetypes=fTyp, initialdir=iDir)
    """
    半値幅計算
    """

    xrdfile = xu.io.panalytical_xml.XRDMLFile(os.path.basename(filepath), path=os.path.dirname(filepath))
    scanmot, intensity = (
        xu.io.panalytical_xml.getxrdml_scan(filetemplate=xrdfile.filename, motors=xrdfile.scan.scanmotname,
                                            path=os.path.dirname(xrdfile.full_filename)))
    count_time = xrdfile.scan.ddict['countTime']
    notzero_index = np.where(intensity > 0)
    scanmot = scanmot / 2
    intensity = np.log10(intensity)

    # 初期値のリストを作成
    # [amp,ctr,wid]
    guess = []
    guess.append([837610, 17.27223, 0.01757])
    guess.append([1437, 17.34145, 0.12478])

    # バックグラウンドの初期値
    background = 1

    # 初期値リストの結合
    guess_total = []
    for i in guess:
        guess_total.extend(i)
    guess_total.append(background)

    y_list = fit_plot(scanmot[notzero_index], *guess_total)
    baseline = np.zeros_like(scanmot[notzero_index]) + np.log10(guess_total[-1])
    for n, i in enumerate(y_list):
        plt.fill_between(scanmot[notzero_index], i, baseline, facecolor=cm.rainbow(n / len(y_list)), alpha=0.6)
    plt.show()

    popt, pcov = curve_fit(func, scanmot[notzero_index], intensity[notzero_index], p0=guess_total)

    fit = func(scanmot, *popt)
    plt.scatter(scanmot, intensity, s=5)
    plt.plot(scanmot, fit, ls='-', c='black', lw=1)

    y_list = fit_plot(scanmot, *popt)
    baseline = np.zeros_like(scanmot) + np.log10(popt[-1])
    for n, i in enumerate(y_list):
        plt.fill_between(scanmot, i, baseline, facecolor=cm.rainbow(n / len(y_list)), alpha=0.6)
    qy = (np.sin(popt[4] * np.pi / 180) + np.sin((popt[4]) * np.pi / 180)) / 2
    c = abs(2 * (1.54 / 2 * 10 ** 10) / qy)

    algan=func(scanmot, *popt[3:6], 0)
    notzero_gaus=np.where(algan>0.012)
    omega1=scanmot[notzero_gaus[0]]
    omega2=scanmot[notzero_gaus[-1]]
    thickness = 1.54 / 10 / 2 / (np.sin(omega2/180/np.pi) - np.sin(omega1/180/np.pi))
    print(popt[3:6])
    print(thickness, c)
    plt.show()


if __name__ == '__main__':
    main()
