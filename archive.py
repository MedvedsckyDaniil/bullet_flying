import numpy as np
import matplotlib.pyplot as plt
import time


def parse(file_name: str):
    file = open(file_name, 'r')
    butcher = np.zeros((16, 16))
    j = 0
    while file:
        line = file.readline()
        if line == "":
            break
        tmp = line.split()
        for i in range(len(tmp)):
            if '/' in tmp[i]:
                ch, zn = tmp[i].split('/')
                butcher[j, i] = float(ch) / float(zn)
            else:
                butcher[j, i] = float(tmp[i])
        j += 1
    butcher[-1] = butcher[j - 1]
    butcher[j - 1] = np.zeros(16)
    return butcher


def cx(mah):
    c_x = 0
    if mah < 0.73:
        c_x = 0.157
    elif mah < 0.82:
        c_x = 0.033 * mah + 0.133
    elif mah < 0.91:
        c_x = 0.161 + 3.9 * (mah - 0.823) ** 2
    elif mah < 1:
        c_x = 1.5 * mah - 1.176
    elif mah < 1.18:
        c_x = 0.386 - 1.6 * (mah - 1.176) ** 2
    elif mah < 1.62:
        c_x = 0.384 * np.sin(1.85 / mah)
    elif mah < 3.06:
        c_x = 0.29 / mah + 0.172
    elif mah < 3.53:
        c_x = 0.316 - 0.016 * mah
    else:
        c_x = 0.259
    # return c_x # (если без нормирования)
    if mah < 0.44:
        return 0.61 * c_x
    elif mah < 0.733:
        return 0.58 * c_x
    elif mah < 0.88:
        return 0.48 * c_x
    elif mah < 1:
        return 0.6 * c_x
    elif mah < 1.173:
        return 0.57 * c_x
    elif mah < 1.466:
        return 0.5 * c_x
    elif mah < 2.053:
        return 0.45 * c_x
    else:
        return 0.48 * c_x


def Force(speed):
    if F_para:
        S = 6.5 * 10 ** (-5)  # S – площадь миделева сечения тела
        i = 0.47  # i – коэффициент формы
        # cx – функция закона сопротивления воздуха;
        a = 335  # a – скорость звука; ну 335 м/сек
        p = 1.2754  # ρ – плотность атмосферы. ну 1,2754 кг/м³
        return i * S * p * speed ** 2 * cx(speed / a) / 2
    else:
        return 0


def air(startX, startY, mass, speed, teta, file):
    teta = np.deg2rad(teta)
    butcher = parse(file)
    print(butcher)
    g = 9.80665
    x = startX
    y = startY
    t = 0
    adge = np.tan(teta)
    Xspeed = speed * np.cos(teta)

    # du_dx = -Force(speed) / (mass * speed)
    # dadge_dx = -g / (Xspeed ** 2)
    # dy_dx = adge
    # dt_dx = 1 / Xspeed
    # speed = Xspeed * np.sqrt(1 + adge ** 2)
    h = 0.002
    # h = 1
    bounds = 10000
    ys = np.zeros(int(bounds / h))
    ys[0] = y
    ts = np.arange(0, bounds, h)

    # Система дифференциальных уравнений:
    # du/dx = -F(v)/(m*v);
    # dγ/dx = -g/(u^2);
    # dy/dx = γ;
    # dt/dx = 1/u;
    # v = u*sqrt(1+γ^2);
    itter = 0
    while y > 0:  # y
        # k1
        dXspeed_dx1 = h * (-Force(speed) / (mass * speed))
        dadge_dx1 = h * (-g / (Xspeed ** 2))
        dy_dx1 = h * (adge)
        dt_dx1 = h / Xspeed
        speed1 = dXspeed_dx1 * np.sqrt(1 + dadge_dx1 ** 2)
        # k2
        sp = speed + speed1 * butcher[1, 0]
        spx = Xspeed + dXspeed_dx1 * butcher[1, 0]
        dXspeed_dx2 = h * (-Force(sp) / (mass * sp))
        dadge_dx2 = h * (-g / ((spx) ** 2))
        dy_dx2 = h * (adge + dadge_dx1 * butcher[1, 0])
        dt_dx2 = h / (spx)
        speed2 = (dXspeed_dx2) * np.sqrt(1 + (dadge_dx2) ** 2)
        # k3
        sp = speed + speed1 * butcher[2, 0] + speed2 * butcher[2, 1]
        spx = Xspeed + dXspeed_dx1 * butcher[2, 0] + dXspeed_dx2 * butcher[2, 1]
        dXspeed_dx3 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx3 = h * (-g / ((spx) ** 2))
        dy_dx3 = h * (adge + dadge_dx1 * butcher[2, 0] + dadge_dx2 * butcher[2, 1])
        dt_dx3 = h / (spx)
        speed3 = (dXspeed_dx3) * np.sqrt(1 + (dadge_dx3) ** 2)
        # k4
        sp = speed + speed1 * butcher[3, 0] + speed2 * butcher[3, 1] + speed3 * butcher[3, 2]
        spx = Xspeed + dXspeed_dx1 * butcher[3, 0] + dXspeed_dx2 * butcher[3, 1] + dXspeed_dx3 * butcher[3, 2]
        dXspeed_dx4 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx4 = h * (-g / ((spx) ** 2))
        dy_dx4 = h * (adge + dadge_dx1 * butcher[3, 0] + dadge_dx2 * butcher[3, 1] + dadge_dx3 * butcher[3, 2])
        dt_dx4 = h / (spx)
        speed4 = (dXspeed_dx4) * np.sqrt(1 + (dadge_dx4) ** 2)
        # k5
        sp = speed + speed1 * butcher[4, 0] + speed2 * butcher[4, 1] + speed3 * butcher[4, 2] + speed4 * butcher[
                4, 3]
        spx = Xspeed + dXspeed_dx1 * butcher[4, 0] + dXspeed_dx2 * butcher[4, 1] + dXspeed_dx3 * butcher[
                4, 2] + dXspeed_dx4 * butcher[4, 3]
        dXspeed_dx5 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx5 = h * (-g / ((spx) ** 2))
        dy_dx5 = h * (adge + dadge_dx1 * butcher[4, 0] + dadge_dx2 * butcher[4, 1] + dadge_dx3 * butcher[
                4, 2] + dadge_dx4 * butcher[4, 3])
        dt_dx5 = h / (spx)
        speed5 = (dXspeed_dx5) * np.sqrt(1 + (dadge_dx5) ** 2)
        # k6
        sp = speed + speed1 * butcher[5, 0] + speed2 * butcher[5, 1] + speed3 * butcher[5, 2] + speed4 * butcher[
                5, 3] + speed5 * butcher[5, 4]
        spx = (Xspeed + dXspeed_dx1 * butcher[5, 0] + dXspeed_dx2 * butcher[5, 1] + dXspeed_dx3 *
                    butcher[5, 2] + dXspeed_dx4 * butcher[5, 3] + dXspeed_dx5 * butcher[5, 4])
        dXspeed_dx6 = h * (-Force(sp) / (mass * sp))
        dadge_dx6 = h * (-g / (spx ** 2))
        dy_dx6 = h * (adge + dadge_dx1 * butcher[5, 0] + dadge_dx2 * butcher[5, 1] + dadge_dx3 * butcher[
                5, 2] + dadge_dx4 * butcher[5, 3] + dadge_dx5 * butcher[5, 4])
        dt_dx6 = h / (spx)
        speed6 = (dXspeed_dx6) * np.sqrt(1 + (dadge_dx6) ** 2)

        # k7
        sp = speed + speed1 * butcher[6, 0] + speed2 * butcher[6, 1] + speed3 * butcher[6, 2] + speed4 * butcher[
                6, 3] + speed5 * butcher[6, 4] + speed6 * butcher[6, 5]
        spx = (Xspeed + dXspeed_dx1 * butcher[6, 0] + dXspeed_dx2 * butcher[6, 1] + dXspeed_dx3 *
                    butcher[6, 2] + dXspeed_dx4 * butcher[6, 3] + dXspeed_dx5 * butcher[6, 4] + dXspeed_dx6 * butcher[
                        6, 5])
        dXspeed_dx7 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx7 = h * (-g / (spx ** 2))
        dy_dx7 = h * (adge + dadge_dx1 * butcher[6, 0] + dadge_dx2 * butcher[6, 1] + dadge_dx3 * butcher[
                6, 2] + dadge_dx4 * butcher[6, 3] + dadge_dx5 * butcher[6, 4] + dadge_dx6 * butcher[6, 5])
        dt_dx7 = h / (spx)
        speed7 = (dXspeed_dx7) * np.sqrt(1 + (dadge_dx7) ** 2)

        # k8
        sp = speed + speed1 * butcher[7, 0] + speed2 * butcher[7, 1] + speed3 * butcher[7, 2] + speed4 * butcher[
                7, 3] + speed5 * butcher[7, 4] + speed6 * butcher[7, 5] + speed7 * butcher[7, 6]
        spx = (Xspeed + dXspeed_dx1 * butcher[7, 0] + dXspeed_dx2 * butcher[7, 1] + dXspeed_dx3 *
                    butcher[7, 2] + dXspeed_dx4 * butcher[7, 3] + dXspeed_dx5 * butcher[7, 4] + dXspeed_dx6 * butcher[
                        7, 5] + dXspeed_dx7 * butcher[7, 6])
        dXspeed_dx8 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx8 = h * (-g / (spx ** 2))
        dy_dx8 = h * (adge + dadge_dx1 * butcher[7, 0] + dadge_dx2 * butcher[7, 1] + dadge_dx3 * butcher[
                7, 2] + dadge_dx4 * butcher[7, 3] + dadge_dx5 * butcher[7, 4] + dadge_dx6 * butcher[7, 5] + dadge_dx7 *
                          butcher[7, 6])
        dt_dx8 = h / (spx)
        speed8 = (dXspeed_dx8) * np.sqrt(1 + (dadge_dx8) ** 2)
        # k9
        sp = speed + speed1 * butcher[8, 0] + speed2 * butcher[8, 1] + speed3 * butcher[8, 2] + speed4 * butcher[
                8, 3] + speed5 * butcher[8, 4] + speed6 * butcher[8, 5] + speed7 * butcher[8, 6] + speed8 * butcher[
                      8, 7]
        spx = (Xspeed + dXspeed_dx1 * butcher[8, 0] + dXspeed_dx2 * butcher[8, 1] + dXspeed_dx3 *
                    butcher[8, 2] + dXspeed_dx4 * butcher[8, 3] + dXspeed_dx5 * butcher[8, 4] + dXspeed_dx6 * butcher[
                        8, 5] + dXspeed_dx7 * butcher[8, 6] + dXspeed_dx8 * butcher[8, 7])
        dXspeed_dx9 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx9 = h * (-g / (spx ** 2))
        dy_dx9 = h * (adge + dadge_dx1 * butcher[8, 0] + dadge_dx2 * butcher[8, 1] + dadge_dx3 * butcher[
                8, 2] + dadge_dx4 * butcher[8, 3] + dadge_dx5 * butcher[8, 4] + dadge_dx6 * butcher[8, 5] + dadge_dx7 *
                          butcher[8, 6] + dadge_dx8 * butcher[8, 7])
        dt_dx9 = h / (spx)
        speed9 = (dXspeed_dx9) * np.sqrt(1 + (dadge_dx9) ** 2)
        # k10
        sp = speed + speed1 * butcher[9, 0] + speed2 * butcher[9, 1] + speed3 * butcher[9, 2] + speed4 * butcher[
                9, 3] + speed5 * butcher[9, 4] + speed6 * butcher[9, 5] + speed7 * butcher[9, 6] + speed8 * butcher[
                       9, 7] + speed9 * butcher[9, 8]
        spx = (Xspeed + dXspeed_dx1 * butcher[9, 0] + dXspeed_dx2 * butcher[9, 1] + dXspeed_dx3 *
                     butcher[9, 2] + dXspeed_dx4 * butcher[9, 3] + dXspeed_dx5 * butcher[9, 4] + dXspeed_dx6 * butcher[
                         9, 5] + dXspeed_dx7 * butcher[9, 6] + dXspeed_dx8 * butcher[9, 7] + dXspeed_dx9 * butcher[
                         9, 8])
        dXspeed_dx10 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx10 = h * (-g / (spx ** 2))
        dy_dx10 = h * (adge + dadge_dx1 * butcher[9, 0] + dadge_dx2 * butcher[9, 1] + dadge_dx3 * butcher[
                9, 2] + dadge_dx4 * butcher[9, 3] + dadge_dx5 * butcher[9, 4] + dadge_dx6 * butcher[9, 5] + dadge_dx7 *
                           butcher[9, 6] + dadge_dx8 * butcher[9, 7] + dadge_dx9 * butcher[9, 8])
        dt_dx10 = h / (spx)
        speed10 = (dXspeed_dx10) * np.sqrt(1 + (dadge_dx10) ** 2)
        # k11
        sp = speed + speed1 * butcher[10, 0] + speed2 * butcher[10, 1] + speed3 * butcher[10, 2] + speed4 * \
                   butcher[
                       10, 3] + speed5 * butcher[10, 4] + speed6 * butcher[10, 5] + speed7 * butcher[10, 6] + speed8 * \
                   butcher[
                       10, 7] + speed9 * butcher[10, 8] + speed10 * butcher[10, 9]
        spx = (Xspeed + dXspeed_dx1 * butcher[10, 0] + dXspeed_dx2 * butcher[10, 1] + dXspeed_dx3 *
                     butcher[10, 2] + dXspeed_dx4 * butcher[10, 3] + dXspeed_dx5 * butcher[10, 4] + dXspeed_dx6 *
                     butcher[
                         10, 5] + dXspeed_dx7 * butcher[10, 6] + dXspeed_dx8 * butcher[10, 7] + dXspeed_dx9 * butcher[
                         10, 8] + dXspeed_dx10 * butcher[10, 9])
        dXspeed_dx11 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx11 = h * (-g / (spx ** 2))
        dy_dx11 = h * (adge + dadge_dx1 * butcher[10, 0] + dadge_dx2 * butcher[10, 1] + dadge_dx3 * butcher[
                10, 2] + dadge_dx4 * butcher[10, 3] + dadge_dx5 * butcher[10, 4] + dadge_dx6 * butcher[
                               10, 5] + dadge_dx7 *
                           butcher[10, 6] + dadge_dx8 * butcher[10, 7] + dadge_dx9 * butcher[10, 8] + dadge_dx10 *
                           butcher[
                               10, 9])
        dt_dx11 = h / (spx)
        speed11 = (dXspeed_dx11) * np.sqrt(1 + (dadge_dx11) ** 2)
        # k12
        sp = speed + speed1 * butcher[11, 0] + speed2 * butcher[11, 1] + speed3 * butcher[11, 2] + speed4 * \
                   butcher[
                       11, 3] + speed5 * butcher[11, 4] + speed6 * butcher[11, 5] + speed7 * butcher[11, 6] + speed8 * \
                   butcher[
                       11, 7] + speed9 * butcher[11, 8] + speed10 * butcher[11, 9] + speed11 * butcher[11, 10]
        spx = (Xspeed + dXspeed_dx1 * butcher[11, 0] + dXspeed_dx2 * butcher[11, 1] + dXspeed_dx3 *
                     butcher[11, 2] + dXspeed_dx4 * butcher[11, 3] + dXspeed_dx5 * butcher[11, 4] + dXspeed_dx6 *
                     butcher[
                         11, 5] + dXspeed_dx7 * butcher[11, 6] + dXspeed_dx8 * butcher[11, 7] + dXspeed_dx9 * butcher[
                         11, 8] + dXspeed_dx10 * butcher[11, 9] + dXspeed_dx11 * butcher[11, 10])
        dXspeed_dx12 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx12 = h * (-g / (spx ** 2))
        dy_dx12 = h * (adge + dadge_dx1 * butcher[11, 0] + dadge_dx2 * butcher[11, 1] + dadge_dx3 * butcher[
                11, 2] + dadge_dx4 * butcher[11, 3] + dadge_dx5 * butcher[11, 4] + dadge_dx6 * butcher[
                               11, 5] + dadge_dx7 *
                           butcher[11, 6] + dadge_dx8 * butcher[11, 7] + dadge_dx9 * butcher[11, 8] + dadge_dx10 *
                           butcher[
                               11, 9] + dadge_dx11 * butcher[11, 10])
        dt_dx12 = h / (spx)
        speed12 = (dXspeed_dx12) * np.sqrt(1 + (dadge_dx12) ** 2)
        # k13
        sp = speed + speed1 * butcher[12, 0] + speed2 * butcher[12, 1] + speed3 * butcher[12, 2] + speed4 * \
                   butcher[
                       12, 3] + speed5 * butcher[12, 4] + speed6 * butcher[12, 5] + speed7 * butcher[12, 6] + speed8 * \
                   butcher[
                       12, 7] + speed9 * butcher[12, 8] + speed10 * butcher[12, 9] + speed11 * butcher[
                       12, 10] + speed12 * butcher[12, 11]
        spx = (Xspeed + dXspeed_dx1 * butcher[12, 0] + dXspeed_dx2 * butcher[12, 1] + dXspeed_dx3 *
                     butcher[12, 2] + dXspeed_dx4 * butcher[12, 3] + dXspeed_dx5 * butcher[12, 4] + dXspeed_dx6 *
                     butcher[
                         12, 5] + dXspeed_dx7 * butcher[12, 6] + dXspeed_dx8 * butcher[12, 7] + dXspeed_dx9 * butcher[
                         12, 8] + dXspeed_dx10 * butcher[12, 9] + dXspeed_dx11 * butcher[12, 10] + dXspeed_dx12 *
                     butcher[12, 11])
        dXspeed_dx13 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx13 = h * (-g / (spx ** 2))
        dy_dx13 = h * (adge + dadge_dx1 * butcher[12, 0] + dadge_dx2 * butcher[12, 1] + dadge_dx3 * butcher[
                12, 2] + dadge_dx4 * butcher[12, 3] + dadge_dx5 * butcher[12, 4] + dadge_dx6 * butcher[
                               12, 5] + dadge_dx7 *
                           butcher[12, 6] + dadge_dx8 * butcher[12, 7] + dadge_dx9 * butcher[12, 8] + dadge_dx10 *
                           butcher[
                               12, 9] + dadge_dx11 * butcher[12, 10] + dadge_dx12 * butcher[12, 11])
        dt_dx13 = h / (spx)
        speed13 = (dXspeed_dx13) * np.sqrt(1 + (dadge_dx13) ** 2)
        # k14
        sp = speed + speed1 * butcher[13, 0] + speed2 * butcher[13, 1] + speed3 * butcher[13, 2] + speed4 * \
                   butcher[
                       13, 3] + speed5 * butcher[13, 4] + speed6 * butcher[13, 5] + speed7 * butcher[13, 6] + speed8 * \
                   butcher[
                       13, 7] + speed9 * butcher[13, 8] + speed10 * butcher[13, 9] + speed11 * butcher[
                       13, 10] + speed12 * butcher[13, 11] + speed13 * butcher[13, 12]
        spx = (Xspeed + dXspeed_dx1 * butcher[13, 0] + dXspeed_dx2 * butcher[13, 1] + dXspeed_dx3 *
                     butcher[13, 2] + dXspeed_dx4 * butcher[13, 3] + dXspeed_dx5 * butcher[13, 4] + dXspeed_dx6 *
                     butcher[
                         13, 5] + dXspeed_dx7 * butcher[13, 6] + dXspeed_dx8 * butcher[13, 7] + dXspeed_dx9 * butcher[
                         13, 8] + dXspeed_dx10 * butcher[13, 9] + dXspeed_dx11 * butcher[13, 10] + dXspeed_dx12 *
                     butcher[13, 11] + dXspeed_dx12 * butcher[13, 12])
        dXspeed_dx14 = h * (-Force(sp) / (mass * (sp)))
        dadge_dx14 = h * (-g / (spx ** 2))
        dy_dx14 = h * (adge + dadge_dx1 * butcher[13, 0] + dadge_dx2 * butcher[13, 1] + dadge_dx3 * butcher[
                13, 2] + dadge_dx4 * butcher[13, 3] + dadge_dx5 * butcher[13, 4] + dadge_dx6 * butcher[
                               13, 5] + dadge_dx7 *
                           butcher[13, 6] + dadge_dx8 * butcher[13, 7] + dadge_dx9 * butcher[13, 8] + dadge_dx10 *
                           butcher[
                               13, 9] + dadge_dx11 * butcher[13, 10] + dadge_dx12 * butcher[13, 11] + dadge_dx13 * butcher[13, 12])
        dt_dx14 = h / (spx)
        speed14 = (dXspeed_dx14) * np.sqrt(1 + (dadge_dx14) ** 2)
        # x
        Xspeed += dXspeed_dx1 * butcher[-1, 0] + dXspeed_dx2 * butcher[-1, 1] + dXspeed_dx3 * butcher[
            -1, 2] + dXspeed_dx4 * butcher[-1, 3] + dXspeed_dx5 * butcher[-1, 4] + dXspeed_dx6 * butcher[
                      -1, 5] + dXspeed_dx7 * butcher[-1, 6] + dXspeed_dx8 * butcher[-1, 7] + dXspeed_dx9 * butcher[
                      -1, 8] + dXspeed_dx10 * butcher[-1, 9] + dXspeed_dx11 * butcher[-1, 10] + dXspeed_dx12 * butcher[
                      -1, 11] + dXspeed_dx13 * butcher[-1, 12] + dXspeed_dx14 * butcher[-1, 13]
        adge += dadge_dx1 * butcher[-1, 0] + dadge_dx2 * butcher[-1, 1] + dadge_dx3 * butcher[-1, 2] + dadge_dx4 * \
                butcher[-1, 3] + dadge_dx5 * butcher[-1, 4] + dadge_dx6 * butcher[-1, 5] + dadge_dx7 * butcher[
                    -1, 6] + dadge_dx8 * butcher[-1, 7] + dadge_dx9 * butcher[-1, 8] + dadge_dx10 * butcher[
                    -1, 9] + dadge_dx11 * butcher[-1, 10] + dadge_dx12 * butcher[-1, 11] + dadge_dx13 * butcher[-1, 12] + dadge_dx14 * butcher[-1, 13]
        y += dy_dx1 * butcher[-1, 0] + dy_dx2 * butcher[-1, 1] + dy_dx3 * butcher[-1, 2] + dy_dx4 * butcher[
            -1, 3] + dy_dx5 * butcher[-1, 4] + dy_dx6 * butcher[-1, 5] + dy_dx7 * butcher[-1, 6] + dy_dx8 * butcher[
                 -1, 7] + dy_dx9 * butcher[-1, 8] + dy_dx10 * butcher[-1, 9] + dy_dx11 * butcher[-1, 10] + dy_dx12 * \
             butcher[-1, 11] + dy_dx13 * butcher[-1, 12] + dy_dx14 * butcher[-1, 13]
        t = dt_dx1 * butcher[-1, 0] + dt_dx2 * butcher[-1, 1] + dt_dx3 * butcher[-1, 2] + dt_dx4 * butcher[
            -1, 3] + dt_dx5 * butcher[-1, 4] + dt_dx6 * butcher[-1, 5] + dt_dx7 * butcher[-1, 6] + dt_dx8 * butcher[
                -1, 7] + dt_dx9 * butcher[-1, 8] + dt_dx10 * butcher[-1, 9] + dt_dx11 * butcher[-1, 10] + dt_dx12 * \
            butcher[-1, 11] + dt_dx13 * butcher[-1, 12] + dt_dx14 * butcher[-1, 13]
        speed += speed1 * butcher[-1, 0] + speed2 * butcher[-1, 1] + speed3 * butcher[-1, 2] + speed4 * butcher[
            -1, 3] + speed5 * butcher[-1, 4] + speed6 * butcher[-1, 5] + speed7 * butcher[-1, 6] + speed8 * butcher[
                     -1, 7] + speed9 * butcher[-1, 8] + speed10 * butcher[-1, 9] + speed11 * butcher[-1, 10] + speed12 * \
                 butcher[-1, 11] + speed13 * butcher[-1, 12] + speed14 * butcher[-1, 13]

        itter += 1
        if y >= 0 and itter < int(bounds / h):
            ys[itter] = y
        else:
            break
    num = min(len(ts), itter)
    return ts[:num], ys[:num]


# int main()
# {
startX = 0
startY = 1.6
mass = 0.0061
speedi = 315
teta = 0
F_para = True
# F_para=False
file = 'butcher/DP8'
first_stamp = time.time()
a = air(startX, startY, mass, speedi, teta, file)
end_stamp = time.time()
print((end_stamp - first_stamp) * 1000)
F_para = False
file = 'butcher/DP8'
first_stamp = time.time()
b = air(startX, startY, mass, speedi, teta, file)
end_stamp = time.time()
print((end_stamp - first_stamp) * 1000)


def dd(i, argX, argY):
    if i == argX:
        return argY
    else:
        return 0


ay_max = max(a[1])
ax_max = a[0][a[1].argmax()]

cy = [dd(i, ax_max, ay_max) for i in a[0]]
by_max = max(b[1])
bx_max = b[0][b[1].argmax()]
dy = [dd(i, bx_max, by_max) for i in b[0]]
plt.plot(a[0], cy, 'g')
plt.plot(b[0], dy, 'y')

plt.plot(a[0], a[1], label='normal')
plt.plot(b[0], b[1], 'r', label='F(air)=0')
plt.legend()
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(a[0], a[1], color='tab:blue', label='angle normal')
ax.plot(b[0], b[1], color='tab:red', label='angle F(air)=0')
ax.set_xlim(0, 2)
ax.set_ylim(startY - 1, startY + 1)
plt.legend()
plt.show()
# return 0
# }
