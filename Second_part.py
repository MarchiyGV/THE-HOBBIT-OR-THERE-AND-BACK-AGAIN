import matplotlib.pyplot as plt
import numpy as np
from math import sin, cos, sqrt, copysign, atan, pi
 
data = np.loadtxt("1_to_2_.txt", delimiter='\t', dtype = np.float)

m = 143 * 10 ** 3  # Масса ракеты в начале перелёта
M_E = 5.97 * 10 ** 24  # Масса Земли
M_M = 7.35 * 10 ** 22  # Масса Луны
R_E = 6375 * 10 ** 3  # Радиус Земли
R_M = 1738 * 10 ** 3  # Радиус Луны
R_orb = 384405 * 10 ** 3  # Радиус орбиты Луны
V_M = 1.023 * 10 ** 3  # Скорость
G = 6.67 * 10 ** (-11)  # Гравитационная постоянная
F = [1016 * 10 ** 3, 97.75 * 10 ** 3]  # Массив со значениями сил тяги разных двигателей 
U = [4130, 3050]  # Массив со скоростями истечения продуктов сгорания
fuel = [F[0] / U[0], F[1] / U[1]]  # Массив со значениями расхода топлива
dt = 0.5  # Промежуток времени для численного интегрирования
Moon_pos = [R_orb, 0]  # Начальное положение Луны
x0 = data[0]
y0 = R_E + data[1]
V_x0 = - sqrt(G * M_E / sqrt(x0 ** 2 + y0 ** 2))  # Проекция скорости ракеты в начале перелёта
V_y0 = 0  # Проекция скорости ракеты в начале перелёта
v_final = sqrt(G * M_M / (R_M + 50 * 10 ** 3))  # Желаемая скорость на окололунной орбите
V = [[V_x0, V_y0]]  # Массив, содержащий массивы проекций скоростей ракеты на протяжении полёта
coord = [[x0, y0]]  # Массив, содержащий массивы координат ракеты на протяжении полёта
MOON_pos = [Moon_pos]  # Аналогичный массив для положения Луны
MOON_V = []
c_in_m_fr = []

def alpha(pos):  # Функция для определения угла между радиус-вектором положения корабля и осью Ox
    if pos[0] == 0:
        return pi / 2 if pos[1] > 0 else 3 * pi / 2
    elif pos[0] > 0:
        return atan(pos[1] / pos[0]) if pos[1] > 0 else 2 * pi + atan(pos[1] / pos[0])
    else:
        return pi + atan(pos[1] / pos[0])

def moon_coord(pos, t):  # Определение координат Луны через время t.
    w = V_M / R_orb
    return [R_orb * cos(alpha(pos) - w * t), R_orb * sin(alpha(pos) - w * t)]

def coord_v_in_dif_fr(r_pos, v, m_pos, frame='Moon'):  # Перевода координат и проекций скорости из одной СО в другую.
    if frame == 'Moon':  # Перевод в СО "Луна"
        new_pos = [r_pos[i] - m_pos[i] for i in range(2)]
        new_v = [v[i] - moon_vel(m_pos)[i] for i in range(2)]
    else:# Перевод в СО "Земля"
        new_pos = [r_pos[i] + m_pos[i] for i in range(2)]
        new_v = [v[i] + moon_vel(m_pos)[i] for i in range(2)]
    return new_pos, new_v

def moon_vel(pos):  # Функция, возвращающая скорость Луны в указанном положении.
    return [V_M * sin(alpha(pos)), - V_M * cos(alpha(pos))]

def arrival(r_pos, v, m_pos, mass, t):  # Изменение траектории ракеты при подлёте к Луне.
    global MOON_pos, coord, Moon_pos
    r_pos, v = coord_v_in_dif_fr(r_pos, v, m_pos, 'Moon')
    while r_pos[0] ** 2 + r_pos[1] ** 2 >= (R_M + 200 * 10 ** 3) ** 2:
        t += dt
        main_angle = alpha(v) + pi
        possible_var = []
        for ang in np.linspace(main_angle - 50 * pi / 180, main_angle + 50 * pi / 180, 500):
            fx, fy = F[1] * cos(ang), F[1] * sin(ang)
            n_pos, n_v = pos_and_velocity(r_pos, v, m_pos, fx, fy, mass, dt, 'Moon')
            sqr_n_v = n_v[0] ** 2 + n_v[1] ** 2
            n_dist = sqrt(n_pos[0] ** 2 + n_pos[1] ** 2)
            n_cos_vect = abs(n_v[0] * n_pos[0] + n_v[1] * n_pos[1]) / sqrt(sqr_n_v * (n_pos[0] ** 2 + n_pos[1] ** 2))
            possible_var.append([ang, n_cos_vect, sqrt(sqr_n_v), n_dist])
        found = 0
        while found == 0 and len(possible_var) != 0:
            optim = possible_var[0]
            i_optim = 0
            for i in range(len(possible_var)):
                if possible_var[i][1] < optim[1] and possible_var[i][2] >= v_final:
                    optim = possible_var[i]
                    i_optim = i
            if optim[2] <= v_final:
                del possible_var[i_optim]
            else:
                found = 1
        if len(possible_var) != 0:
            fx, fy = F[1] * cos(optim[0]), F[1] * sin(optim[0])
            mass -= fuel[1] * dt
            r_pos, v = pos_and_velocity(r_pos, v, m_pos, fx, fy, mass, dt, 'Moon')
            coord.append(coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')[0])
            V.append(coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')[1])
        else:
            r_pos, v = pos_and_velocity(r_pos, v, m_pos, 0, 0, mass, dt, 'Moon')
            m_pos = moon_coord(m_pos, dt)
            coord.append(coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')[0])
            V.append(coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')[1])
    return coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')[0], coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')[1], mass, t

def pos_and_velocity(r_pos, v, m_pos, fx, fy, m, dt, frame='Earth'):  # Суммирование
    # Вычисление положения ракеты и проекций её скоростей через время dt в различных СО.
    x, y = r_pos[0], r_pos[1]
    if frame == 'Earth':
        x_m, y_m = m_pos[0], m_pos[1]
        a_x = fx / m - G * M_E * x / (x ** 2 + y ** 2) ** (3 / 2) - G * M_M * (x - x_m) / ((x - x_m) ** 2 +
                                                                                       (y - y_m) ** 2) ** (3 / 2)
        a_y = fy / m - G * M_E * y / (x ** 2 + y ** 2) ** (3 / 2) - G * M_M * (y - y_m) / ((x - x_m) ** 2 +
                                                                                       (y - y_m) ** 2) ** (3 / 2)
    else:
        a_x = fx / m - G * M_M * x / (x ** 2 + y ** 2) ** (3 / 2)
        a_y = fy / m - G * M_M * y / (x ** 2 + y ** 2) ** (3 / 2)
    v_new = [v[0] + a_x * dt, v[1] + a_y * dt]
    pos_new = []
    for i in range(2):
        pos_new.append(r_pos[i] + (v[i] + v_new[i]) * 0.5 * dt)
    return pos_new, v_new

def flying(r_pos, v, m_pos, mass, fx, fy, detachment): # Функция, осуществляющая полет
    if mass >= (20 + 5.5 + 22.5 + 15) * 10 ** 3 + fuel[0] * dt:
        r_pos, v = pos_and_velocity(r_pos, v, m_pos, fx, fy, mass, dt)
        mass -= fuel[0] * dt
    else:
        if detachment == 0:
            mass -= 20 * 10 ** 3
            detachment = 1
        if (fx ** 2 + fy ** 2 - F[1] ** 2) < 1:
            r_pos, v = pos_and_velocity(r_pos, v, m_pos, fx, fy, mass, dt)
            mass -= fuel[1] * dt
        else:
            r_pos, v = pos_and_velocity(r_pos, v, m_pos, 0, 0, mass, dt)
    m_pos = moon_coord(m_pos, dt)
    if abs(fx ** 2 + fy ** 2 - F[1] ** 2) < 1:
        r_pos, v = pos_and_velocity(r_pos, v, m_pos, fx, fy, mass, dt)
        mass -= fuel[0] * dt
    return r_pos, v, m_pos, mass, detachment

def approx_time_of_flight_and_rot(start_pos, start_v):  # Поиск примерного времени полёта до Луны
    # и угла поворота радиус-вектора ракеты.
    global MOON_pos, coord, Moon_pos
    r_pos = [0, sqrt(start_pos[0] ** 2 + start_pos[1] ** 2)]
    start_angle = alpha(r_pos)
    v = [sqrt(start_v[0] ** 2 + start_v[1] ** 2), 0]
    t = 0
    local_m = m
    detachment = 0
    m_pos = [0, R_orb]
    while local_m >= (20 + 5.5 + 22.5 + 15) * 10 ** 3 + fuel[0] * dt:
        fx = F[0] * v[0] / sqrt(v[0] ** 2 + v[1] ** 2)
        fy = F[0] * v[1] / sqrt(v[0] ** 2 + v[1] ** 2)
        r_pos, v, m_pos, local_m, detachment = flying(r_pos, v, m_pos, local_m, fx, fy, detachment)
        m_pos = [0, R_orb]
        t += dt
    while r_pos[0] ** 2 + r_pos[1] ** 2 < R_orb ** 2:
        r_pos, v, m_pos, local_m, detachment = flying(r_pos, v, m_pos, local_m, 0, 0, detachment)
        m_pos = [0, R_orb]
        t += dt
    finish_angle = alpha(r_pos)
    if finish_angle - start_angle > 0:
        return t, finish_angle - start_angle
    else:
        return t, 2 * pi + finish_angle - start_angle

def waiting(r_pos, v, m_pos, time_and_rot):  # Поиск времени ожидания на орбите Земли для оптимального полёта.
    global MOON_pos, coord, Moon_pos
    waiting_time = 0
    f_time, rot = time_and_rot[0], time_and_rot[1]
    st_angle = alpha(r_pos)
    while (cos(st_angle - rot) * R_orb - moon_coord(m_pos, f_time )[0]) ** 2 + \
            (sin(st_angle - rot) * R_orb - moon_coord(m_pos, f_time)[1]) ** 2 > 10 ** 10:
        m_pos = moon_coord(m_pos, dt)
        r_pos, v = pos_and_velocity(r_pos, v, m_pos, 0, 0, m, dt)
        st_angle = alpha(r_pos)
        waiting_time += dt
        MOON_pos.append(m_pos)
        coord.append(r_pos)
        V.append(v)
    return waiting_time, r_pos, v

def flight(r_pos, v, m_pos):
    global m, MOON_pos, coord, Moon_pos, MOON_V, c_in_m_fr
    print('СТАРТ ЭТАПА №2. Начинаем подготовку к перелету на Луну.')
    t, r_pos, v = waiting(r_pos, v, m_pos, approx_time_of_flight_and_rot(r_pos, v))
    m_pos = moon_coord(m_pos, t)
    detachment = 0
    print('1) Достигнута наиболее выгодная точка на орбите. Корабль покидает орбиту Земли.')
    print("С момента начала ожидания прошло:", "%.0f" % (t), "секунд")
    while sqrt((r_pos[0] - m_pos[0]) ** 2 + (r_pos[1] - m_pos[1]) ** 2) > 2 * R_M:
        fx = F[0] * v[0] / sqrt(v[0] ** 2 + v[1] ** 2)
        fy = F[0] * v[1] / sqrt(v[0] ** 2 + v[1] ** 2)
        r_pos, v, m_pos, m, detachment = flying(r_pos, v, m_pos, m, fx, fy, detachment)
        t += dt
        MOON_pos.append(m_pos)
        V.append(v)
        coord.append(r_pos)
        MOON_V.append(moon_vel(m_pos))
    r_pos, v = coord_v_in_dif_fr(r_pos, v, m_pos, 'Moon')
    print('2) Корабль приближается к Луне. Начинается заход на лунную орбиту.')
    print("С момента начала маневра прошло:", "%.0f" % (t), "секунд;", "Осталось топлива:", "%.0f" % (m), "килограмм")
    while v[0] ** 2 + v[1] ** 2 >= v_final ** 2 * 2:
        fx = - F[1] * v[0] / sqrt(v[0] ** 2 + v[1] ** 2)
        fy = - F[1] * v[1] / sqrt(v[0] ** 2 + v[1] ** 2)
        m -= fuel[1] * dt
        r_pos, v = pos_and_velocity(r_pos, v, m_pos, fx, fy, m, dt, 'Moon')
        m_pos = moon_coord(m_pos, dt)
        coord.append(coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')[0])
        V.append(coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')[1])
        t += dt
    print('3) Начинается активная стадия захода на лунную орбиту.')
    print("С момента начала маневра прошло:", "%.0f" % (t), "секунд;", "Осталось топлива:", "%.0f" % (m), "килограмм")
    r_pos, v = coord_v_in_dif_fr(r_pos, v, m_pos, 'Earth')
    r_pos, v, m, t = arrival(r_pos, v, m_pos, m, t)
    print('4) Корабль зашел на лунную орбиту. Миссия выполнена.')
    print("С момента начала маневра прошло:", "%.0f" % (t), "секунд;", "Осталось топлива:", "%.0f" % (m), "килограмм")
    while t < 180000:
        r_pos, v, m_pos, m, detachment = flying(r_pos, v, m_pos, m, 0, 0, detachment)
        t += dt
        MOON_pos.append(m_pos)
        V.append(v)
        coord.append(r_pos)
        MOON_V.append(moon_vel(m_pos))
        c_in_m_fr.append(coord_v_in_dif_fr(r_pos, v, m_pos, 'Moon')[0])
    return

flight(coord[0], V[0], Moon_pos)
fig = plt.figure(num = None, figsize = (12, 9), dpi = 100)
ax1 = fig.add_subplot(221)
ax1.set_title(r'Движение системы в СО "Земля"')
ax1.set_xlabel(r'$x$, м')
ax1.set_ylabel(r'$y$, м')
plt.grid()
ax1.plot([coord[i][0] for i in range(len(coord))], [coord[i][1] for i in range(len(coord))])
ax1.plot([MOON_pos[i][0] for i in range(len(MOON_pos))], [MOON_pos[i][1] for i in range(len(MOON_pos))])
ax2 = fig.add_subplot(222)
ax2.set_title(r'Скорость ракеты в СО "Земля"')
ax2.set_xlabel(r'$v_x$, м/c')
ax2.set_ylabel(r'$v_y$, м/c')
plt.grid()
ax2.plot([V[i][0] for i in range(len(V))], [V[i][1] for i in range(len(V))])
plt.show()
    
x_fM, y_fM = MOON_pos[-1]
x_f, y_f = coord[-1]
x_fM = str(x_fM / 1000)
y_fM = str(y_fM / 1000)
x_f = str(x_f / 1000)
y_f = str(y_f / 1000)
f = open('2_to_3.txt', 'w')
f.write(x_fM + '\t' + y_fM + '\t' + x_f + '\t' + y_f)
f.close()