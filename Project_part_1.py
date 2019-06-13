import numpy as np
import matplotlib.pyplot as plt
import math as m

from scipy.integrate import ode
h=0
Re=6371000
G=6.67e-11
Me=5.97e24
coeff_Cx=0.1
totalTime=0

S1_fuelMass=2010000
S1_emptyMass=135000
S1_Force=34350000
S1_Vfuel=2580
S1_Diametr=10.1

S2_fuelMass=458700-37600
S2_emptyMass=37600
S2_Force=5115000
S2_Vfuel=4130
S2_Diametr=10.1

S3_fuelMass=100000
S3_emptyMass=20000
S3_Force=1016000
S3_Vfuel=4130
S3_Diametr=6.6

SpaceshipMass=5500+22500+15000
TotalMass=SpaceshipMass+S3_fuelMass+S2_fuelMass+S1_fuelMass+S3_emptyMass+S2_emptyMass+S1_emptyMass


def fout1(t, y):# обработчик шага 
    ts.append(t)
    ys.append(list(y.copy()))
    y1, y2, y3, y4 = y
    if (S1_fuelMass-S1_Force/S1_Vfuel*t<=0):
        return -1
 
def fout2(t, y):# обработчик шага 
    ts.append(t)
    ys.append(list(y.copy()))
    y1, y2, y3, y4 = y
    if (S2_fuelMass-S2_Force/S2_Vfuel*t<=0):
        return -1    
 
def fout3(t, y):# обработчик шага 
    ts.append(t)
    ys.append(list(y.copy()))
    y1, y2, y3, y4 = y
    if (y4>=m.sqrt(G*Me/(m.sqrt(y1*y1+y3*y3)))):
        return -1 

def fout4(t,y):
    ts.append(t)
    ys.append(list(y.copy()))
    y1, y2, y3, y4 = y
    if (m.sqrt(y1*y1+y3*y3)-Re<185000):
        print("Кажется, у нас проблемы")
        return -1
    if y1<0:
        return -1
    
        
def rho(x,y):
    '''с помощью линейной апроксимации определяет плотность воздуха на необходимой нам высоте'''
    Space = [0, 1.85 * 0.00001, 1.5*0.0001, 3 * 0.0001, 1.03 * 0.001, 4 * 0.001, 7.26 * 0.001, 
             0.0136, 0.0251, 0.0469, 0.0889, 0.1216, 0.1665, 0.2279, 0.3119, 0.3648, 0.4135, 0.4671, 
             0.5258, 0.59, 0.6601, 0.7365, 0.8194, 0.9093, 1,1.1]
            # плотность для разных высот
    Space_lst = [100000, 80000, 70000, 60000, 50000, 40000, 36000, 32000, 28000, 24000, 20000, 
                 18000, 16000, 14000, 12000, 11000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000,2000,1000]    
   
    h = m.sqrt(x*x +y*y) - Re
    
    i = 25
    while h > Space_lst[i]:
        i-=1
        if i < 0:
            return 0
           
    delta = h - Space_lst[i]
# разница между высотой и ближайшим значением
    delta_h = Space_lst[i-1] - Space_lst[i]
# разница между ближайшими соседями
    otn = delta/delta_h
# относительное отклонение
    p = Space[i] + ((Space[i-1] - Space[i]) * otn)
    return p

def launcher_angle(x,y):
    Space_lst = [180000,150000,100000,80000,  60000,  40000,  30000, 15000, 12000,10000,9000,6000,2000,0]
    angle_lst = [70, 68, 66,64,63,61,59,57,55,45,40,25,20,0]
    h = m.sqrt(x*x +y*y) - Re
    i = 0
    while h < Space_lst[i]:
        i+=1
        if i >=13:
           break
    return angle_lst[i]*m.pi/180

def launcher_angle2(x,y,t):
    if x>0:
        if (70 + t)*m.pi/180 < m.atan(y/x) + m.pi/2:
            return ( 70+ t)*m.pi/180
        else:
            return m.atan(y/x) + m.pi/2 
    else:
        if (70 + t)*m.pi/180 < m.pi - m.atan(y/-x) + m.pi/2:
            return (70 + t)*m.pi/180
        else:
            return m.pi - m.atan(y/-x) + m.pi/2

# функция правых частей системы ОДУ
def f1(t, y):
    global TotalMass,S1_fuelMass    
    y1, y2, y3, y4 = y                  
    dM_fuel=S1_Force/S1_Vfuel*t
    jet_a=S1_Force/(TotalMass-dM_fuel)
    vv=y2*y2+y4*y4
    angle=launcher_angle(y1,y3)         
    S=m.pi*S1_Diametr**2/4
    resistant_a=coeff_Cx*rho(y1,y3)*vv*S/(TotalMass- dM_fuel)/2
    ax=-y1*G*Me/((y1*y1+y3*y3)**1.5)+(jet_a-resistant_a)*m.cos(angle)
    ay=-y3*G*Me/((y1*y1+y3*y3)**1.5)+(jet_a-resistant_a)*m.sin(angle)
    return [y2,ax, y4,ay] 

def f2(t, y):
    global TotalMass,S2_fuelMass    
    y1, y2, y3, y4 = y                  
    dM_fuel=S2_Force/S2_Vfuel*t
    jet_a=S2_Force/(TotalMass-dM_fuel)         
    angle=launcher_angle(y1,y3)          
    vv=y2*y2+y4*y4
    S=m.pi*S2_Diametr**2/4
    resistant_a=coeff_Cx*rho(y1,y3)*vv*S/(TotalMass- dM_fuel)/2
    ax=-y1*G*Me/((y1*y1+y3*y3)**1.5)+(jet_a-resistant_a)*m.cos(angle)
    ay=-y3*G*Me/((y1*y1+y3*y3)**1.5)+(jet_a-resistant_a)*m.sin(angle)
    return [y2,ax, y4,ay]

def f3(t, y):
    global TotalMass,S3_fuelMass    
    y1, y2, y3, y4 = y                  
    dM_fuel=S3_Force/S3_Vfuel*t
    jet_a=S3_Force/(TotalMass-dM_fuel)                             
    angle=launcher_angle2(y1,y3,t)                 
    vv=y2*y2+y4*y4
    S=m.pi*S3_Diametr**2/4
    resistant_a=coeff_Cx*rho(y1,y3)*vv*S/(TotalMass- dM_fuel)/2
    ax=-y1*G*Me/((y1*y1+y3*y3)**1.5)+(jet_a-resistant_a)*m.cos(angle)
    ay=-y3*G*Me/((y1*y1+y3*y3)**1.5)+(jet_a-resistant_a)*m.sin(angle)                          
    return [y2,ax, y4,ay]   

def f4(t,y):
    global TotalMass    
    y1, y2, y3, y4 = y                 
    dM_fuel=0                    
    angle = launcher_angle2(y1,y3,t+t3)         
    vv=y2*y2+y4*y4
    S=m.pi*S3_Diametr**2/4
    resistant_a=coeff_Cx*rho(y1,y3)*vv*S/(TotalMass- dM_fuel)/2
    ax=-y1*G*Me/((y1*y1+y3*y3)**1.5)-resistant_a*m.cos(angle)
    ay=-y3*G*Me/((y1*y1+y3*y3)**1.5)-resistant_a*m.sin(angle)                 
    return [y2,ax, y4,ay]   
       

x_start=Re
Vx_start=0
y_start=0
Vy_start=2*m.pi/(24*23600) * Re 

xc,yc=[],[]
for i in range(0, 630):
    xc.append(Re*m.cos(i/100))
    yc.append(Re*m.sin(i/100))

print("Первая ступень")
y0,t0=[x_start,  Vx_start, y_start, Vy_start], 0 # начальные условия 
ODE=ode(f1)
ODE.set_integrator('dopri5')#, max_step=0.01)
ODE.set_solout(fout1)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(S1_fuelMass/(S1_Force/S1_Vfuel))      # решение ОДУ
Y=np.array(ys)
plt.plot(Y[:,0],Y[:,2],linewidth=3)#,label='k=%.1f'% k)
print(Y[-1:,0],Y[-1:,2],Y[-1:,1],Y[-1:,3]  )
print(m.sqrt(Y[-1:,0]**2 +Y[-1:,2]**2) - Re)
print(ts[-1]);
print(TotalMass)
totalTime+=ts[-1]
TotalMass-=(S1_fuelMass+S1_emptyMass)
print(TotalMass)
print()

print("Вторая ступень")
y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3]], 0 # начальные условия 
ODE=ode(f2)
ODE.set_integrator('dopri5')
ODE.set_solout(fout2)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(S2_fuelMass/(S2_Force/S2_Vfuel))      # решение ОДУ
Y=np.array(ys)
plt.plot(Y[:,0],Y[:,2],linewidth=3)#,label='k=%.1f'% k)
print(Y[-1:,0],Y[-1:,2],Y[-1:,1],Y[-1:,3]  )
print(m.sqrt(Y[-1:,0]**2 +Y[-1:,2]**2) - Re)
print(ts[-1]);
print(TotalMass)
totalTime+=ts[-1]
TotalMass-=(S2_fuelMass+S2_emptyMass)
print(TotalMass)
print()

print("Третья ступень")
t3=255.2
y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3]], 0 # начальные условия
ODE=ode(f3)
ODE.set_integrator('dopri5')#, max_step=0.01)
ODE.set_solout(fout3)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(t3)      # решение ОДУ
Y=np.array(ys)
plt.plot(Y[:,0],Y[:,2],linewidth=3)#,label='k=%.1f'% k)
print(Y[-1:,0],Y[-1:,2],Y[-1:,1],Y[-1:,3]  )
print(m.sqrt(Y[-1:,0]*Y[-1:,0] + Y[-1:,2]*Y[-1:,2]) - Re )
print(ts[-1]);
print(TotalMass)
totalTime+=ts[-1]
TotalMass-=((S3_Force/S3_Vfuel)*ts[-1])
print(TotalMass)
print()

print("Свободный полет")
y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3]], 0 # начальные условия
ODE=ode(f4)
ODE.set_integrator('dopri5',max_step=1.883913)
ODE.set_solout(fout4)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(10000)      # решение ОДУ
Y=np.array(ys)
plt.plot(Y[:,0],Y[:,2],linewidth=3)#,label='k=%.1f'% k)
plt.axis('equal')
plt.plot(xc,yc,linewidth=2)
plt.title("Траектория до перехода ко второму этапу")
plt.grid(True)  
plt.show()
print(Y[-1:,0],Y[-1:,2],Y[-1:,1],Y[-1:,3]  )
print(m.sqrt(Y[-1:,0]*Y[-1:,0] + Y[-1:,2]*Y[-1:,2]) - Re )
print(ts[-1]);
print(TotalMass)
totalTime+=ts[-1]
TotalMass-=(0)
print(TotalMass)
print(totalTime)
print()

#Конечные данные и вывод в файл:
print("Вот это мои конечные данные")
print("x =", Y[-1:,0][0])
print("y =", Y[-1:,2][0])
print("h =",m.sqrt(Y[-1:,0]*Y[-1:,0]+Y[-1:,2]*Y[-1:,2])-Re)
print("v_x =",Y[-1:,1][0])
print("v_y =",Y[-1:,3][0])
print("Mass =",TotalMass)
print("Time =", totalTime)

l = [Y[-1:,0][0],Y[-1:,2][0],Y[-1:,1][0],Y[-1:,3][0],TotalMass,totalTime]


f = open('1_to_2.txt', 'w') 
f.write(str(l)) 
f.close()