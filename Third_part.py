from math import sqrt, pi, acos, cos, sin

G = 6.67e-11  # Гравитационная постоянная
M = 7.36e22  # Масса Луны
Radius = 1738000  # Радиус Луны
mass_space_ship = 6835 # Сухая масса 
mass_engine = 8135  # Масса топлива
connection_of_lm_and_lk = True
F_crit = 45040 # Максимальная тяга двигателя
u_fuel = 3050 # Скорость истечения продуктов
int_period = 0.01 # Малый период времени

string = open('2_to_3.txt').readlines()
m = list([[float(i) for i in string[k].split()] for k in range((len(string)))])

x_moon = m[0][0]
y_moon = m[0][1]
x_ship = m[0][2]
y_ship = m[0][3]

result = open('output landing moon.txt', 'w')
result.close()

def changing_system_of_coordinates(x_moon, y_moon, x_ship,
                                   y_ship):  # Функция перехода к системе координат с точкой отсчета в центре луны, OY сонапрвлена с радиус-вектором, соединяющим землю и луну
    cos = y_moon / sqrt(x_moon ** 2 + y_moon ** 2)
    sin = x_moon / sqrt(x_moon ** 2 + y_moon ** 2)
    x_move = x_ship - x_moon
    y_move = y_ship - y_moon
    x_turn = x_move * cos - y_move * sin
    y_turn = x_move * sin + y_move * cos
    return x_turn, y_turn


x_turn, y_turn = changing_system_of_coordinates(x_moon, y_moon, x_ship, y_ship)  # Переход в СО "Луна"
r_orbit = sqrt(x_turn ** 2 + y_turn ** 2) * 1000  # Радиус лунной орбиты 
orbital_speed = sqrt(G * M / r_orbit) # Подсчет скорости вращения по орбите


def a_x_and_a_y(x, y, grav_const=G,
                mass_moon=M):  # Определяет значения гравитационного ускорения
    a_x = - (grav_const * mass_moon * x) / ((x ** 2 + y ** 2) ** 1.5)
    a_y = - (grav_const * mass_moon * y) / ((x ** 2 + y ** 2) ** 1.5)
    return a_x, a_y


def speed_and_coordinate(x, y, v_x, v_y, a_x, a_y,
                         t = int_period):  # Определяет изменение скорости и координаты после одной итерации
    x = x + v_x * t + a_x * (t ** 2) / 2
    y = y + v_y * t + a_y * (t ** 2) / 2
    v_x = v_x + a_x * t
    v_y = v_y + a_y * t
    return x, y, v_x, v_y


def free_flight(x, y, v_x, v_y, coordinate = 0):  # Осуществляет свободный полет
    while x <= coordinate:
        a_x, a_y = a_x_and_a_y(x, y)
        x, y, v_x, v_y = speed_and_coordinate(x, y, v_x, v_y, a_x, a_y)
    return x, y, v_x, v_y


def speed_after_impulse(x = 0, y = r_orbit, v_x = - orbital_speed, v_y = 0,
                        rad=Radius):  # Определяет конечную скорость для перехода на эллиптическую орбиту
    height_finish = y
    min_speed = v_x # Минимально возможная скорость
    x_probe = x
    y_probe = y
    v_y_probe = v_y
    while height_finish > rad + ((y_probe - rad) / 2):
        min_speed += 0.25
        x = x_probe
        y = y_probe
        v_x = min_speed
        v_y = v_y_probe
        x, y, v_x, v_y = free_flight(x, y, v_x, v_y)
        height_finish = abs(y)
    return min_speed

def defining_coordinate_of_giving_impulse(orb_speed=orbital_speed):  # Определяет координату х импульса
    speed = speed_after_impulse()
    x_impulse = (-speed ** 2 + orb_speed ** 2) / 6  # В формуле учитывается максимальная тяга двигателя
    return x_impulse

def position_on_the_orbit(x, y, radius):  # Переопределяет позицию ЛМ на орбите через длину дуги
    if x <= 0:
        alpha = acos(y / radius)
    else:
        alpha = acos(- y / radius) + pi
    position = radius * alpha  # Длина дуги от верхней точки
    return position


def speed_and_coordinate_after_impulse(x_ship_orbit, y_ship_orbit, orb_speed = orbital_speed, rad = r_orbit, t = int_period,
    m_ship = mass_space_ship, m_engine = mass_engine, u = u_fuel):  # Определение скорости и координаты после подачи импульса
    x = defining_coordinate_of_giving_impulse()
    y = sqrt(rad ** 2 - x ** 2)
    position_of_impulse = position_on_the_orbit(x, y, rad)
    position_of_ship = position_on_the_orbit(x_ship_orbit, y_ship_orbit, rad)
    if position_of_impulse >= position_of_impulse:
        t_init = (position_of_impulse - position_of_ship) / orb_speed
    else:
        t_init = (position_of_impulse + 2 * pi * rad - position_of_ship) / orb_speed
    t_new = 0
    v_x = - (orb_speed * y) / rad
    v_y = (orb_speed * x) / rad
    a_x = 3
    result = open('output landing moon.txt', 'a')
    while x > 0:
        a_grav_x, a_y = a_x_and_a_y(x, y)
        a_engine = a_x - a_grav_x
        fuel_consumption = (m_ship + m_engine) * a_engine / u
        x, y, v_x, v_y = speed_and_coordinate(x, y, v_x, v_y, a_x, a_y)
        m_engine -= fuel_consumption * t
        t_new += t
        t_new = round(t_new * 100) / 100
        result.write(str(x) + '\t' + str(y) + '\t' + str(sqrt(x ** 2 + y ** 2)) + '\t ' + str(
            sqrt(v_x ** 2 + v_y ** 2)) + '\t' + str(sqrt(a_x ** 2 + a_y ** 2)) + '\t' + str(t_new) + '\n')
    result.close()
    return m_engine, v_x, v_y, x, y, t_init, t_new, position_of_impulse
m_engine_after_impuls, v_x, v_y, x, y, t_init, t_new, position_of_impulse = speed_and_coordinate_after_impulse(x_turn, y_turn)

print('1) Сообщение импульса.')
print('Осталось топлива: ', m_engine_after_impuls, 'кг;')
print('Абсцисса: ', x, 'км;')
print('Ордината: ', y, 'км;')
print('Скорость на Ох: ', v_x, 'км/с;')
print('Скорость на Оу: ', v_y, 'км/с;')


def writing_coordinates(x, y, v_x, v_y, mass_engine, t_new, f_max=F_crit,
                        t=int_period):  # Свободное падение по эллиптической орбите
    n = 1
    result = open('output landing moon.txt', 'a')
    while y > 0:  # Половина траектории
        a_x, a_y = a_x_and_a_y(x, y)
        x, y, v_x, v_y = speed_and_coordinate(x, y, v_x, v_y, a_x, a_y)
        t_new += t
        t_new = round(t_new * 100) / 100
        n += 1
        if n % 10 == 0:
            result.write(str(x) + '\t' + str(y) + '\t' + str(sqrt(x ** 2 + y ** 2)) + '\t' + str(
                sqrt(v_x ** 2 + v_y ** 2)) + '\t' + str(sqrt(a_x ** 2 + a_y ** 2)) + '\t' + str(t_new) + '\n')
            # Для записи в файл
    a_x, a_y = a_x_and_a_y(x, y)
    while (v_x ** 2) < abs(2 * x * (- (f_max / (
            mass_space_ship + mass_engine)) + a_x)):  # Рассчитывается скорость для которой погасится горизонтальная проекция
        a_x, a_y = a_x_and_a_y(x, y)
        x, y, v_x, v_y = speed_and_coordinate(x, y, v_x, v_y, a_x, a_y)
        t_new += t
        t_new = round(t_new * 100) / 100
        n += 1
        if n % 10 == 0:
            result.write(str(x) + '\t' + str(y) + '\t' + str(sqrt(x ** 2 + y ** 2)) + '\t' + str(
                sqrt(v_x ** 2 + v_y ** 2)) + '\t' + str(sqrt(a_x ** 2 + a_y ** 2)) + '\t' + str(t_new) + '\n')
    result.close()
    return x, y, v_x, v_y, t_new


x_st, y_st, v_x_st, v_y_st, t_new = writing_coordinates(x, y, v_x, v_y, m_engine_after_impuls, t_new)
print('2) В конце свободного падения')
print('Абсцисса: ', x_st, 'км;')
print('Ордината: ', y_st, 'км;')
print('Скорость на Ох: ', v_x_st, 'км/с;')
print('Скорость на Оу: ', v_y_st, 'км/с;')


def making_stop(x, y, v_x, v_y, m_engine, t_new, f_max=F_crit, m_ship=mass_space_ship, u=u_fuel,
                t=int_period):  # осуществление торможения
    a_x_grav, a_y_grav = a_x_and_a_y(x, y)
    a_x = - (F_crit / (m_engine + m_ship)) + a_x_grav
    result = open('output landing moon.txt', 'a')
    while v_x > 0:
        a_x_grav, a_y_grav = a_x_and_a_y(x, y)
        a_x_engine = - a_x + a_x_grav
        a_y_engine = - sqrt(((f_max / (m_ship + m_engine)) ** 2) - (a_x_engine ** 2))
        a_y = a_y_grav + a_y_engine
        if m_engine < 7800:
            if v_y > 100:
                m_engine -= (f_max / u) * t
            else:
                a_y = a_y_grav
                m_engine -= (a_x_engine * (m_ship + m_engine) / u) * t
        else:
            m_engine -= (f_max / u) * t
        t_new += t
        t_new = round(t_new * 100) / 100
        x, y, v_x, v_y = speed_and_coordinate(x, y, v_x, v_y, a_x, a_y)
        result.write(str(x) + '\t' + str(y) + '\t' + str(sqrt(x ** 2 + y ** 2)) + '\t' + str(
            sqrt(v_x ** 2 + v_y ** 2)) + '\t' + str(sqrt(a_x ** 2 + a_y ** 2)) + '\t' + str(t_new) + '\n')
        # Запись в файл координат и скоростей
    result.close()
    return x, y, v_x, v_y, m_engine, t_new


print('3) Перед началом посадки')
x_fin, y_fin, v_x_fin, v_y_fin, m_engine_fin, t_new = making_stop(x_st, y_st, v_x_st, v_y_st, m_engine_after_impuls,

                                                                  t_new)
print('Осталось топлива: ', m_engine_fin, 'кг;')
print('Абсцисса: ', x_fin, 'км;')
print('Ордината: ', y_fin, 'км;')
print('Скорость на Ох: ', v_x_fin, 'км/с;')
print('Скорость на Оу: ', v_y_fin, 'км/с;')


def total_stop(x, y, v_x, v_y, m_engine, t_new, u=u_fuel, m_ship=mass_space_ship, rad=Radius,
               t=int_period):  # Посадка
    result = open('output landing moon.txt', 'a')
    while v_y > 0:
        a_x_grav, a_y_grav = a_x_and_a_y(x, y)
        a_y = - (v_y ** 2 / (2 * abs(abs(y) - rad)))
        a_engine = a_y_grav - a_y
        fuel_consumption = a_engine * (m_ship + m_engine) / u
        x, y, v_x, v_y = speed_and_coordinate(x, y, v_x, v_y, a_x_grav, a_y)
        m_engine -= fuel_consumption * t
        t_new += t
        t_new = round(t_new * 100) / 100
        result.write(str(x) + '\t' + str(y) + '\t' + str(sqrt(x ** 2 + y ** 2)) + '\t' + str(
            sqrt(v_x ** 2 + v_y ** 2)) + '\t' + str(sqrt(a_x_grav ** 2 + a_y ** 2)) + '\t' + str(t_new) + '\n')
        # Запись в файл координат и скоростей
    result.close()
    return x, y, v_x, v_y, m_engine, t_new


print('4) Заканчиваем посадку')
x_landing, y_landing, v_x_landing, v_y_landing, mass_engine_left, t_new = total_stop(x_fin, y_fin, v_x_fin, v_y_fin,
                                                                                     m_engine_fin, t_new)
print('Координаты прилунения: x = ' + str(x_landing) + ' y = ' + str(y_landing) + ' ; посадочные скорости: v_x = ' + str(
    v_x_landing) + ' v_y = ' + str(v_y_landing) + ' ; остаток топлива = ' + str(-mass_engine_left))

t_overall = t_init + t_new
position_of_lk = (position_of_impulse + orbital_speed * t_new) % (2 * pi * r_orbit)
y_lk = r_orbit * cos(position_of_lk / r_orbit)
x_lk = - r_orbit * sin(position_of_lk / r_orbit)  # определение координат лунного корабля после посадки

print('Координаты ЛК на орбите: x = ' + str(x_lk) + ' y = ' + str(y_lk) + ' ; Время посадки = ' + str(
    t_overall))
result = open('output landing moon2.txt', 'a')
result.write('Coordinates landing:' + '\n' + 'x = ' + str(x_landing) + '\n' + 'y = ' + str(
    y_landing) + '\n' + 'speed landing:' + '\n' + 'v_x = ' + str(v_x_landing) + '\n' + ' v_y = ' + str(
    v_y_landing) + '\n' + 'fuel left = ' + str(mass_engine_left) + '\n')
result.write('Coordinates of moon ship on orbit:' + '\n' 'x = ' + str(x_lk) + '\n' + 'y = ' + str(
    y_lk) + '\n' + 'time_of_the_stage = ' + str(t_overall))
result.close()

import matplotlib.pyplot as plt  # Визуализация
import pylab
from numpy import *

string = open('output_landing_moon.txt').readlines()
m = array([[float(i) for i in string[k].split()] for k in range((len(string)))])
# порядок записи в файл -- координата x, y, sqrt(x^2+y^2), sqrt(Vx^2+Vy^2), sqrt(ax^2+ay^2), t, координата x, y ЛК на орбите
from matplotlib.pyplot import *

plt.title(' Положение ЛМ в пространстве относительно луны ', size=11)
plot(list(m[:, 0] / 1000), list(m[:, 1] / 1000), "blue", markersize=0.1)
plt.xlabel('Координата x, км')
plt.ylabel('Координата y, км')
pylab.xlim(-2000, 2000)
pylab.ylim(-2000, 2000)
plt.grid()
show()

plt.title(' r(t) ', size=11)
plot(list(m[:, 5]), list(m[:, 2] / 1000), "blue", markersize=0.1)
plt.ylabel('Расстояние до центра луны, км ')
plt.xlabel('Время, c')
plt.grid()
show()

plt.title(' V(t) ', size=11)
plot(list(m[:, 5]), list(m[:, 3]), "blue", markersize=0.1)
plt.ylabel('Скорость, м/с ')
plt.xlabel('Время, с')
plt.grid()
show()

x_lk = str(x_lk) 
y_lk = str(y_lk) 
t_overall = str(t_overall) 
f = open('3_to_4.txt', 'w') 
f.write(x_lk + '\t' + y_lk + '\t' + t_overall) 
f.close()
