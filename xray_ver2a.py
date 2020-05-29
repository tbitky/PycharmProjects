import os
import tkinter
import re
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


def distinguish_hkl_location_number(x, y):
    """
    ファイル名がいつものヤツか判別する(ndim=1)とファイル名からサンプル番号、測定した結晶面、測定位置を抽出(ndim=2)。
    """
    kara = ()
    if np.ndim(np.array(x)) == 2:
        if re.match(r"G\d+\sOmega\s\d{3}_upper_\d\.xrdml", x[y][1]):
            kara = re.findall(r"G(\d+)\sOmega\s(\d{3})_upper_(\d)\.xrdml", x[y][1])
        elif re.match(r"G\d+\s1point_omega\s\d{3}\supper_\d\.xrdml", x[y][1]):
            kara = re.findall(r"G(\d+)\s1point_omega\s(\d{3})\supper_(\d)\.xrdml", x[y][1])
        else:
            messagebox.showinfo('エラー', x[y][1] + 'のファイル名はエラーです')
            kara = False
    elif np.ndim(np.array(x)) == 1:
        if (re.match(r"G\d+\sOmega\s\d{3}_upper_\d\.xrdml", os.path.basename(x[y]))
                or re.match(r"G\d+\s1point_omega\s\d{3}\supper_\d\.xrdml", os.path.basename(x[y]))):
            kara = True
        else:
            messagebox.showinfo('エラー', x[y] + 'のファイル名はエラーです')
            kara = False
    return kara


def fine_round(x, y=0):
    round_x = Decimal(str(x)).quantize(Decimal(str(y)), rounding=ROUND_HALF_UP)
    int_x = int(round_x)
    return int_x


def check_123_or_134(xarray, i):
    """
    ファイル名が1,2,3か1,3,4か区別する。
    """
    for n in range(xarray.shape[0]):
        if xarray[n][2] == xarray[i][2]:
            if xarray[n][3] == xarray[i][3]:
                if xarray[n][4] == "4":
                    return True
    return False


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
    ファイル名がいつものでないやつを弾いてリスト化
    """
    fList = np.zeros((len(filepath), 6), dtype="U39")
    errornr = 0
    for i in range(len(filepath)):
        if distinguish_hkl_location_number(filepath, i):
            fList[i][0] = os.path.dirname(filepath[i])
            fList[i][1] = os.path.basename(filepath[i])
    temp = fList
    fList = np.delete(temp, np.where(temp[:, 0] == ""), axis=0)  # fListの無駄な配列をカット

    """
    半値幅計算
    """
    for num in range(len(fList)):
        xrdfile = xu.io.panalytical_xml.XRDMLFile(fList[num][1], path=fList[num][0])
        scanmot, inte = (
            xu.io.panalytical_xml.getxrdml_scan(filetemplate=xrdfile.filename, motors=xrdfile.scan.scanmotname,
                                                path=os.path.dirname(xrdfile.full_filename)))
        integer_inte = np.array(inte * 10).astype(np.int)

        fList[num, 2:5] = np.array(distinguish_hkl_location_number(fList, num))

        try:
            peaks, properties = find_peaks(integer_inte, distance=(len(integer_inte) / 2))
            widths = peak_widths(integer_inte, peaks)
            fwhm = scanmot[fine_round(widths[3][0])] - scanmot[fine_round(widths[2][0])]
            fList[num][5] = fwhm * 3600
        except TypeError:
            print(fList[num][1] + "は半値幅が計算できませんでした。")
            fList[num][5] = ""

    """
    Excelシートの初期化
    """
    wb = openpyxl.Workbook()
    sheet = wb.active
    minsample = np.min(np.array(fList[:, 2], dtype=np.int32), axis=0)
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
    for i in range(len(fList)):
        datacolumn = 0
        sheet.cell(row=int(fList[i][2]) - minsample + 3, column=1).value = int(fList[i][2])
        if fList[i][3] == "002":
            datacolumn = 2
        elif fList[i][3] == "102":
            datacolumn = 6
        elif fList[i][3] == "100":
            datacolumn = 10
        else:
            messagebox.showinfo('エラー', fList[i][1] + 'は(002),(102),(100)ではありませんね？')

        if check_123_or_134(fList, i):
            if fList[i][4] == "3" or fList[i][4] == "4":
                datacolumn += int(fList[i][4]) - 2
            else:
                datacolumn += int(fList[i][4]) - 1
        else:
            datacolumn += int(fList[i][4]) - 1

        if fList[i][5] != "":
            sheet.cell(row=int(fList[i][2]) - minsample + 3, column=datacolumn).value = float(
                "{:.1f}".format(float(fList[i][5])))
            if float(fList[i][5]) > 4000:  # 4000s以上を誤差注意とする
                fill = openpyxl.styles.PatternFill(patternType="solid", fgColor="FFFF00", bgColor="FFFF00")
                sheet.cell(row=int(fList[i][2]) - minsample + 3, column=1).fill = fill
        else:
            sheet.cell(row=int(fList[i][2]) - minsample + 3, column=datacolumn).value = ""
            fill = openpyxl.styles.PatternFill(patternType="solid", fgColor="FFFF00", bgColor="FFFF00")
            sheet.cell(row=int(fList[i][2]) - minsample + 3, column=1).fill = fill
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
