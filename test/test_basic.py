import lib.egw.functions as egw
def main():
    qx, qy = map(float, input('qx[rlu] qy入力[rlu]:').split())
    hh, kk, ll = map(int, input('h k l入力:').split())
    egw.ternary_a_c_r_calculate(qx, qy, hh, kk, ll)


if __name__ == '__main__':
    main()
