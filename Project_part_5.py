import numpy as np
import matplotlib.pyplot as plt
import math as m

from scipy.integrate import ode

'''
Физические константы
'''
G=6.67e-11
Re=6371000
Me=5.97e24
R=384405000
Rm = 1738000
Mm = 7.35e22
Vm = 1022
Wm = Vm/R
'''
Параметры ЛК
'''
fuelMass= 17700
emptyMass= 22500 - 17700
Force= 95750
Vfuel= 3050
TotalMass = 5500 + 22500
'''
Пользовательские параметры и параметры системы
'''
h1 = 70000
phi0 = 0
moon_crush1 = False #индикаторы врезаний
moon_crush2 = False
moon_crush3 = False
earth_crush3 = False
moon_crush4 = False
earth_crush4 = False

'''
Входные данные (СО Луны)
'''
f = open('4_to_5.txt', 'r')
input_data = f.read()
f.close()
input_data = input_data.replace("[", "")
input_data = input_data.replace("]", "")
input_data = input_data.replace(" ", "")
x, Vx, y, Vy, totalTime = input_data.split(",")
x = float(x)
Vx = float(Vx)
y = float(y)
Vy = float(Vy)
totalTime = float(totalTime)
'''
Перевод в СО Земли и расчёты
'''
def moon(t):    #координаты луны
        x = R*m.cos(Wm*t + phi0)
        y = R*m.sin(Wm*t + phi0)
        vx = -Vm*m.sin(Wm*t + phi0)
        vy = Vm*m.cos(Wm*t + phi0)
        
        return [x, vx, y, vy]
    
xm, vxm, ym, vym = moon(totalTime)
x_start,  Vx_start, y_start, Vy_start = x + xm, Vx + vxm, y + ym, Vy + vym

h2 = (x**2 + y**2)**0.5 - Rm
R1 = Re + h1
R2 = Rm + h2
VI_m = (G*Mm/R2)**0.5
'''
Обработчики шагов интегрирования
'''
def fout1(t, y):# ожидание
        global moon_crush1
        e = 0.1
        ts.append(t)
        ys.append(list(y.copy()))
        y1, y2, y3, y4 = y
        Xm, VXm, Ym, VYm = moon(totalTime + t)
        moon_crush1 = False
        if ((y1-Xm)**2 + (y3-Ym)**2)**0.5 - Rm <= 0:
          moon_crush1 = True
          return -1
        if abs(m.atan2(-y4, -y2) - m.atan2(y3, y1)) < e:
          print(y4/y2, y3/y1)
          return -1
      
def fout2(t, y):# разгон
        global moon_crush2
        t1s.append(t)
        y1s.append(list(y.copy()))
        y1, y2, y3, y4 = y
        Xm, VXm, Ym, VYm = moon(totalTime + t)
        moon_crush2 = False
        if ((y1-Xm)**2 + (y3-Ym)**2)**0.5 - Rm <= 0:
          moon_crush2 = True
          return -1
        if (fuelMass-t*Force/Vfuel<=0):
          return -1
        
        
def fout3(t, y):# перелёт
        global moon_crush3
        global earth_crush3
        t2s.append(t)
        y2s.append(list(y.copy()))
        y1, y2, y3, y4 = y
        R2s.append((y1**2 + y3**2)**0.5)
        Xm, VXm, Ym, VYm = moon(totalTime + t)
        moon_crush3 = False
        earth_crush3 = False
        
        if ((y1-Xm)**2 + (y3-Ym)**2)**0.5 - Rm <= 0:
          moon_crush3 = True
          return -1
        if (y1**2 + y3**2)**0.5 - Re <= 0:
          earth_crush3 = True
          return -1
        if (y1**2 + y3**2)**0.5 > 10*min(R2s):
          return -1
      
def fout4(t, y):# визуализация финишной орбиты
        global moon_crush4
        global earth_crush4
        ts.append(t)
        ys.append(list(y.copy()))
        y1, y2, y3, y4 = y
        Xm, VXm, Ym, VYm = moon(totalTime + t)
        moon_crush4 = False
        earth_crush4 = False
        if ((y1-Xm)**2 + (y3-Ym)**2)**0.5 - Rm <= 0:
          moon_crush4 = True
          return -1
        if (y1**2 + y3**2)**0.5 - Re <= 0:
          earth_crush4 = True
          return -1
        
'''
Функции
'''        
def launcher_angle_to_Earth(x,y): #направление на землю
    
    return m.atan2(-y, -x)

        
# функции правых частей системы ОДУ
    
def f1(t, y):  #полёт по инерции
         y1, y2, y3, y4 = y
         xm, vxm, ym, vym = moon(t+totalTime)
         
         ax=-y1*G*Me/((y1*y1+y3*y3)**1.5) + (xm-y1)*G*Mm/(((xm-y1)*(xm-y1)+(ym-y3)*(ym-y3))**1.5)
         ay=-y3*G*Me/((y1*y1+y3*y3)**1.5) + (ym-y3)*G*Mm/(((xm-y1)*(xm-y1)+(ym-y3)*(ym-y3))**1.5)
         
         return [y2,ax, y4,ay] 

def f2(t, y):   #разгон с орбиты луны
         y1, y2, y3, y4 = y
         xm, vxm, ym, vym = moon(t+totalTime)      
         
         dM_fuel=Force/Vfuel*t
         jet_a=Force/(TotalMass-dM_fuel)
         angle=launcher_angle_to_Earth(y1,y3)
         
         ax=-y1*G*Me/((y1*y1+y3*y3)**1.5) + (xm-y1)*G*Mm/(((xm-y1)*(xm-y1)+(ym-y3)*(ym-y3))**1.5) + jet_a*m.cos(angle)
         ay=-y3*G*Me/((y1*y1+y3*y3)**1.5) + (ym-y3)*G*Mm/(((xm-y1)*(xm-y1)+(ym-y3)*(ym-y3))**1.5) + jet_a*m.sin(angle)
         return [y2,ax, y4,ay]       
       
'''
Создание поверхности Земли
'''
xe,ye=[],[]
for i in range(0, 630):
    xe.append(Re*m.cos(i/100))
    ye.append(Re*m.sin(i/100))        

'''
Ожидание момента, когда скорость ЛК будет направлена на землю
'''
print("t =", totalTime)
print("Mass =", TotalMass)

y0,t0=[x_start,  Vx_start, y_start, Vy_start], 0
ODE=ode(f1)
ODE.set_integrator('dopri5')
ODE.set_solout(fout1)
ODE.set_initial_value(y0, t0)
ts, ys = [ ],[ ]
tau_up = 2*m.pi/Wm
ODE.integrate(tau_up)
Y=np.array(ys)

totalTime+=ts[-1]

Xm, VXm, Ym, VYm = moon(totalTime)
xm,ym=[],[]

for i in range(0, 630):
    xm.append(Xm + Rm*m.cos(i/100))
    ym.append(Ym + Rm*m.sin(i/100))

plt.plot(Y[:,0],Y[:,2],linewidth=3)
plt.axis('equal')
plt.plot(xe,ye,linewidth=2)
plt.plot(xm,ym,linewidth=2)
plt.title("Waiting \n ")

plt.grid(True)

plt.show()
print("x =", Y[-1:,0], "\ny =", Y[-1:,2], "\nVx =", Y[-1:,1], "\nVy =", Y[-1:,3])
print("dt =", ts[-1])
print("t =", totalTime)
h = (((Xm-Y[-1:,0])*(Xm-Y[-1:,0])+(Ym-Y[-1:,2])*(Ym-Y[-1:,2]))**0.5)-Rm
print("h_moon =", h)
print("V =", (Y[-1:,1]*Y[-1:,1] + Y[-1:,3]*Y[-1:,3])**0.5)

'''
Разгон
'''

k_bottom = 0
k_top = 1
R_finish = R1
e = R_finish

while e > 1000:
    k = (k_bottom + k_top)/2
    
    y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3]], 0 
    ODE1=ode(f2)
    ODE1.set_integrator('dopri5')
    ODE1.set_solout(fout2)
    ODE1.set_initial_value(y0, t0) 
    t1s, y1s = [ ],[ ]
    ODE1.integrate((fuelMass/(Force/Vfuel))*k)
    Y1=np.array(y1s)
    
    '''
    Перелёт
    '''
    y0,t0=[Y1[-1:,0],Y1[-1:,1],Y1[-1:,2],Y1[-1:,3]], t1s[-1]
    ODE2=ode(f1)
    ODE2.set_integrator('dopri5')
    ODE2.set_solout(fout3)
    ODE2.set_initial_value(y0, t0)
    t2s, y2s, R2s = [ ],[ ],[ ]
    ODE2.integrate(1000000)
    R_min = min(R2s)
    e = abs(R_finish - R_min)
    if k == 1:
        break
    if R_min < R_finish:
        k_top = k
    elif R_min > R_finish:
        k_bottom = k
    else:
        break

TotalMass-=fuelMass*k
fuelMass*=(1-k)

Y2=np.array(y2s)
h_min = R_min - Re
index_of_R_min = R2s.index(R_min)
tau = t2s[index_of_R_min]


    
totalTime += t1s[-1]

xm,ym=[],[]
Xm, VXm, Ym, VYm = moon(totalTime)

for i in range(0, 630):
    xm.append(Xm + Rm*m.cos(i/100))
    ym.append(Ym + Rm*m.sin(i/100))

plt.plot(Y1[:,0],Y1[:,2],linewidth=3)
plt.axis('equal')
plt.plot(xm,ym,linewidth=1)
plt.title("Launch \n ")

plt.grid(True)

  
plt.show()
print("x =", Y1[-1:,0], "\ny =", Y1[-1:,2], "\nVx =", Y1[-1:,1], "\nVy =", Y1[-1:,3])
print("dt =", t1s[-1])
print("t =", totalTime)
h = (((Xm-Y1[-1:,0])*(Xm-Y1[-1:,0])+(Ym-Y1[-1:,2])*(Ym-Y1[-1:,2]))**0.5)-Rm
print("h_moon =", h)
print("V =", (Y1[-1:,1]*Y1[-1:,1] + Y1[-1:,3]*Y1[-1:,3])**0.5)
print("Mass =", TotalMass)
print("fuel: dM/M =", k)

'''
Вывод полёта
'''
totalTime += tau

xm,ym=[],[]
Xm, VXm, Ym, VYm = moon(totalTime)

for i in range(0, 630):
    xm.append(Xm + Rm*m.cos(i/100))
    ym.append(Ym + Rm*m.sin(i/100))


plt.plot(Y2[0:index_of_R_min,0],Y2[0:index_of_R_min,2],linewidth=3)
plt.axis('equal')
plt.plot(xe,ye,linewidth=1)
plt.plot(xm,ym,linewidth=1)
plt.title("Flight \n ")

plt.grid(True)

  
plt.show()
print("x =", Y2[index_of_R_min,0],"\ny =", Y2[index_of_R_min,2],"\nVx =", Y2[index_of_R_min,1],"\nVy =",Y2[index_of_R_min,3]  )
print("dt =", tau);
print("t =", totalTime);
h = ((Y2[index_of_R_min,0]**2+Y2[index_of_R_min,2]**2)**0.5)-Re
print("h_earth =", h)
print("V =", (Y2[index_of_R_min,1]**2 + Y2[index_of_R_min,3]**2)**0.5)


'''
Вывод выходных данных
'''
f = open('5_to_6.txt', 'w')
f.write(str(Y2[index_of_R_min]) + ", " + str(totalTime) + ", " + str(TotalMass) + ", " + str(fuelMass))
f.close()

'''
Визуализация получившейся орбиты
'''

y0,t0=[Y2[index_of_R_min,0],Y2[index_of_R_min,1],Y2[index_of_R_min,2],Y2[index_of_R_min,3]], 0
ODE=ode(f1)
ODE.set_integrator('dopri5')
ODE.set_solout(fout4)
ODE.set_initial_value(y0, t0)
ts, ys = [ ],[ ]
ODE.integrate(1000000)
Y=np.array(ys)
totalTime+=ts[-1]
Xm, VXm, Ym, VYm = moon(totalTime)
xm,ym=[],[]

for i in range(0, 630):
    xm.append(Xm + Rm*m.cos(i/100))
    ym.append(Ym + Rm*m.sin(i/100))

plt.plot(Y[:,0],Y[:,2],linewidth=3)
plt.axis('equal')
plt.plot(xe,ye,linewidth=2)
plt.title("Orbit \n ")

plt.grid(True)

  
plt.show()

if moon_crush1 or moon_crush2 or moon_crush3 or moon_crush4 :
    print("Moon crush!!!")
    
if earth_crush3 or earth_crush4:
    print("Earth crush!!!")