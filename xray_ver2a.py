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


def distinguish_hkl_location_number(x):
    """
    ファイル名からサンプル番号、測定した結晶面、測定位置を抽出。
    """
    match_template = map(re.compile, [r"G\d+\sOmega\s\d{3}_upper_\d\.xrdml",
                                      r"G\d+\s1point_omega\s\d{3}\supper_\d\.xrdml"])
    for line in match_template:
        if re.match(line, x):
            sample_name, plane, position = re.findall('\d+', x)
            sample_name = int(sample_name)
            position = int(position)
            break
    else:
        messagebox.showinfo('エラー', x + 'のファイル名はエラーです')
        sample_name, plane, position = None, None, None

    return sample_name, plane, position


def fine_round(x, y=0):
    round_x = Decimal(str(x)).quantize(Decimal(str(y)), rounding=ROUND_HALF_UP)
    int_x = int(round_x)
    return int_x



def main():
    """
    ファイルダイアログを開いてfilepathに格納
    """
    root = tkinter.Tk()
    root.withdraw()
    fTyp = [("xrdファイル", "*.xrdml")]
    iDir = os.path.abspath(os.path.dirname(r"\\Altair\SR4000data\003_XRD\ "))
    filepath = filedialog.askopenfilenames(filetypes=fTyp, initialdir=iDir)

    """
    ファイルをデータフレーム化
    """
    Data = pd.DataFrame(columns=('directory', 'file_name', 'sample_name', 'plane', 'position', 'FWHM'))
    for i in range(len(filepath)):
        file_name = os.path.basename(filepath[i])
        sample_name, plane, position = distinguish_hkl_location_number(file_name)
        Data.at[i, 'directory'] = os.path.dirname(filepath[i])
        Data.at[i, 'file_name'] = file_name
        Data.at[i, 'sample_name'] = sample_name
        Data.at[i, 'plane'] = plane
        Data.at[i, 'position'] = position
    """
    半値幅計算
    """
    for i in range(len(Data)):
        xrdfile = xu.io.panalytical_xml.XRDMLFile(Data.at[i, 'file_name'], path=Data.at[i, 'directory'])
        scanmot, intensity = (
            xu.io.panalytical_xml.getxrdml_scan(filetemplate=xrdfile.filename, motors=xrdfile.scan.scanmotname,
                                                path=os.path.dirname(xrdfile.full_filename)))
        count_time = xrdfile.scan.ddict['countTime']
        integer_intensity = np.array(intensity / count_time).astype(np.int)
        try:
            peaks, properties = find_peaks(integer_intensity, distance=(len(integer_intensity) / 2))
            widths = peak_widths(integer_intensity, peaks)
            fwhm = np.abs(scanmot[fine_round(widths[3][0])] - scanmot[fine_round(widths[2][0])])
            Data.at[i, 'FWHM'] = fwhm * 3600
        except TypeError:
            print(Data[i][1] + "は半値幅が計算できませんでした。")
            Data.at[i, 'FWHM'] = None

    """
    Excelシートの初期化
    """
    wb = openpyxl.Workbook()
    sheet = wb.active
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
