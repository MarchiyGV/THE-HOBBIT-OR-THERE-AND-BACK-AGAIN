import matplotlib.pyplot as plt
import math
import numpy as np 
g=1.62
r1=65899
r2=1738000
R=r1+r2
V1=math.sqrt(g*R)
dm=0.014
dm1=0.008
dm2=0.005
M=4670
dt=0.001
m=15000
v=3050
data = np.loadtxt("3_to_4_.txt", delimiter='\t', dtype = np.float)
def go_orbit():
    t0=21600
    t1=6567
    x=data[0]
    y=data[1]
    t=data[2]+t0+t1
    alpha=math.atan(y/x)
    V=0
    M=4670
    mass01=[]
    mass02=[]
    while V<50:
        M=M-dm2
        V=V-v*math.log(1-dm/m)
        Vx=V*math.cos(alpha)-g*math.cos(alpha)*dt
        Vy=V*math.sin(alpha)-g*math.sin(alpha)*dt
        x=x+Vx*dt
        y=y+Vy*dt
        t=t+dt
        mass01.append(x)
        mass02.append(y)
    if V>=50:
        while V<900:
            M=M-dm2
            V=V-v*math.log(1-dm2/M)
            Vx=V*math.cos(alpha)-g*math.cos(alpha)*dt
            Vy=V*math.sin(alpha)-g*math.sin(alpha)*dt
            x=x+Vx*dt
            y=y+Vy*dt
            t=t+dt
            mass01.append(x)
            mass02.append(y)
            alpha=alpha+0.0065*dt
    if V>=900:
        while V<V1:
            M=M-dm2
            V=V-v*math.log(1-dm2/M)
            Vx=V*math.cos(alpha)-g*math.cos(alpha)*dt
            Vy=V*math.sin(alpha)-g*math.sin(alpha)*dt
            x=x+Vx*dt
            y=y+Vy*dt
            mass01.append(x)
            mass02.append(y)
            alpha=alpha+0.002*dt
            t=t+dt
    mass0=[]
    mass0.append(V)
    mass0.append(x)
    mass0.append(y)
    mass0.append(t)
    mass0.append(M)
    mass0.append(-math.pi/2+alpha)
    mass0.append(math.sqrt(x**2+y**2)-r2)
    return [mass01, mass02, mass0]


a=[]
o=[]
b=[]
c=[]
#picture for landing
for i in range(-100000, 800000, 30):
    a.append(i)
for i in range(-R, R, 30):
    o.append(i)
for i in range(len(a)):
    b.append(-math.sqrt(r2**2-(a[i]**2)))
for i in range(len(a)):
    c.append(-math.sqrt(R**2-(a[i]**2)))
fig = plt.figure()
plt.plot(a, b, color='r')
plt.plot(a, c, color='r')


go_orbit()
#picture for going to orbit
y=go_orbit()[0]
u=go_orbit()[1]
plt.plot(y, u, color='k')
TR=go_orbit()[2]

print("WHERE I AM ON ORBIT")
print("V=", '%.2f' % TR[0], "x=", '%.2f' % TR[1], "y=", '%.2f' % TR[2], "alpha=", '%.2f' % TR[5], "M=", '%.2f' % TR[4], "t=", '%.2f' % TR[3], "r=", '%.2f' % TR[6])
f = open('4_to_5.txt', 'w') 
f.write(str([TR[1], -TR[0]*math.sin(TR[5]), TR[2], TR[0]*math.cos(TR[5]), TR[3]]))
f.close()



