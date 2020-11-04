import os
import tkinter
import re
import pandas as pd
import xrayutilities as xu
import numpy as np
from tkinter import filedialog
from tkinter import messagebox
from scipy.signal import find_peaks, peak_widths
import openpyxl
from openpyxl.styles.borders import Border, Side
from openpyxl.styles.alignment import Alignment
import subprocess
from time import sleep
from decimal import Decimal, ROUND_HALF_UP
import psutil
from scipy import signal


def distinguish_sample_number_and_direction(x):
    """
    ファイル名からサンプル番号,測定方向を抽出。
    """
    match_template = re.compile(r"フィルタ+ 傾斜G/d_/d.xlsx")
    for line in match_template:
        if re.match(line, x):
            nums = re.findall('\d+', x)
            sample_name = int(nums[0])
            direction = int(nums[1])
            break
    else:
        sample_name, direction = None, None

    return sample_name, direction


def fine_round(x, y=0):
    round_x = Decimal(str(x)).quantize(Decimal(str(y)), rounding=ROUND_HALF_UP)
    int_x = int(round_x)
    return int_x


def bowing_direction(data):
    maxid = signal.argrelmax(data, order=int(len(data) / 10))
    minid = signal.argrelmin(data, order=int(len(data) / 10))
    if len(maxid) + len(minid) == 1:
        if maxid:
            direction = '凸'
        else:
            direction = '凹'
    elif len(maxid) + len(minid) == 2:
        direction = 'N'
    elif len(maxid) + len(minid) == 3:
        if len(maxid) >= len(minid):
            direction = 'M'
        else:
            direction = 'W'
    else:
        direction = False
    return direction


def switching_fwhm_calculate(position, intensity, means='absolute'):
    peaks, properties = find_peaks(intensity, distance=len(intensity) // 2)
    if means == 'absolute':
        try:
            data = (position, intensity)
            fwhm_in_degree = xu.math.misc.fwhm_exp(position, intensity)
        except:
            fwhm_in_degree = None
    elif means == 'relative':
        try:
            widths = peak_widths(intensity, peaks)
            left_edge = position[int(widths[2][0])] + (widths[2][0] - int(widths[2][0])) * (
                    position[int(widths[2][0]) + 1] - position[int(widths[2][0])])
            right_edge = position[int(widths[3][0])] + (widths[3][0] - int(widths[3][0])) * (
                    position[int(widths[3][0]) + 1] - position[int(widths[3][0])])
            fwhm_in_degree = np.abs(left_edge - right_edge)
        except TypeError:
            fwhm_in_degree = None
    if fwhm_in_degree:
        fwhm_in_arcsec = fwhm_in_degree * 3600
    else:
        fwhm_in_arcsec = fwhm_in_degree
    return fwhm_in_arcsec


def main():
    """
    ファイルダイアログを開いてfilepathに格納
    """
    root = tkinter.Tk()
    root.withdraw()
    fTyp = [("excelファイル", "*.xls")]
    iDir = os.path.abspath(os.path.dirname(r"\\Sirius\carbon-nas\SR4000\001_bowing"))
    filepath = filedialog.askopenfilenames(filetypes=fTyp, initialdir=iDir)

    """
    ファイルをデータフレーム化
    """

    Data = pd.DataFrame(columns=('directory', 'file_name', 'sample_name', 'bowing_1', 'shape_1', 'bowing_2', 'shape_2'))
    for i in range(len(filepath)):
        datawb=openpyxl.load_workbook(filepath[i],read_only=True)
        datasheet = wb['測定結果']
        max_column=len(datasheet.iter_cols(min_col=1, min_row=4, max_row=4))
        for col in datasheet.iter_cols(min_col=1, min_row=4, max_row=4):
            col=x
        file_name = os.path.basename(filepath[i])
        sample_name, direction = distinguish_sample_number_and_direction(file_name)
        Data.at[i, 'directory'] = os.path.dirname(filepath[i])
        Data.at[i, 'file_name'] = file_name
        Data.at[i, 'sample_name'] = sample_name

    """
    Excelシートの初期化
    """
    wb = openpyxl.load_workbook(filepath[1])
    sheet = wb['測定結果']
    for col in sheet.iter_cols(min_col=1, min_row=4, max_row=4):
        for j in np.arange(1, 12):
            if j == 1:
                sheet.merge_cells(start_row=1, end_row=2, start_column=j, end_column=j)
            elif j == 2 or j == 6:
                sheet.merge_cells(start_row=1, end_row=1, start_column=j, end_column=j + 3)

    sheet["A1"].value = "サンプル番号"
    sheet["B1"].value = "002"
    sheet["F1"].value = "102"
    sheet["j1"].value = "100"
    for k in np.arange(1, 12):
        if k == 2 or k == 6:
            sheet.cell(row=2, column=k).value = "センター"
        elif k == 3 or k == 7:
            sheet.cell(row=2, column=k).value = "ミドル"
        elif k == 4 or k == 8:
            sheet.cell(row=2, column=k).value = "エッジ"
        elif k == 5 or k == 9:
            sheet.cell(row=2, column=k).value = "平均"
    """
    シートに記入
    """
    complete_data = ~(Data.isnull().any(axis=1))
    maxsample = np.max(Data.loc[complete_data, 'sample_name'])
    minsample = np.min(Data.loc[complete_data, 'sample_name'])
    error = 0
    for i in Data.index.values:
        if complete_data[i]:
            row = int(Data.at[i, 'sample_name']) - minsample + 3
            sheet.cell(row=row, column=1).value = int(Data.at[i, 'sample_name'])
            for j, k in enumerate(['002', '102', '100']):
                if Data.at[i, 'plane'] == k:
                    datacolumn = 2 + j * 4
            unique_positions = np.array(list((map(int, Data.loc[(Data['sample_name'] == Data.at[i, 'sample_name']) & (
                    Data['plane'] == Data.at[i, 'plane']), 'position'].unique()))))
            unique_positions = np.sort(unique_positions)
            for j, k in enumerate(unique_positions):
                if Data.at[i, 'position'] == k:
                    datacolumn += j
        else:
            row = maxsample - minsample + 4 + error
            sheet.cell(row=row, column=1).value = Data.at[i, 'file_name']
            datacolumn = 2
            error += 1
        if Data.at[i, 'FWHM'] is not None:
            sheet.cell(row=row, column=datacolumn).value = float(
                "{:.1f}".format(float(Data.at[i, 'FWHM'])))
            if float(Data.at[i, 'FWHM']) > 4000:  # 4000s以上を誤差注意とする
                fill = openpyxl.styles.PatternFill(patternType="solid", fgColor="FFFF00", bgColor="FFFF00")
                sheet.cell(row=row, column=1).fill = fill

    """
    枠線とか整える
    """
    side = Side(style="thin", color="000000")
    border = Border(top=side, bottom=side, left=side, right=side)
    align = Alignment(horizontal="center", vertical="center", wrap_text=False)
    for col in sheet.iter_cols(min_col=1):
        for cell in col:
            cell.border = border
            cell.alignment = align
            if cell.column > 1:
                cell.number_format = '0.0'
    """
    ファイルを開く
    """
    user_directory = os.path.expanduser('~')
    excel_file = user_directory + '\AppData\Local\Temp\X線の半値幅リスト.xlsx'
    wb.save(excel_file)
    subprocess.Popen([r"C:\Program Files\Microsoft Office\root\Office16\EXCEL.EXE", excel_file])

    """
    Excel閉じたらファイル削除
    """
    sleep(5)
    while True:
        flag = True
        ps = psutil.process_iter(attrs=['name'])
        for i, p in enumerate(ps):
            pi = p.info
            if pi['name'] == 'EXCEL.EXE':
                flag = False
        if flag:
            break
        sleep(1)
    os.remove(excel_file)


if __name__ == '__main__':
    main()