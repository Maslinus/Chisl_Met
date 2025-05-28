from cProfile import label
from pdb import run
from re import A
from turtle import color, title
import numpy as np
from matplotlib import pyplot as plt

h=0.001 # Шаг
xLeft=0
xRight=5
y1 = list()
y2 = list()
x = 0
z  = list()
z2 = list()
n=(int)((xRight-xLeft)/h)
t = np.arange(xLeft,xRight,h)

x0 = 0.6
Teta0 = 0.25
Gamma = 10000
Beta = 0.8
Alfa = 0.5
Teta_C = 0
B = 10
Da = np.arange(0.05, 0.12, 0.01)

def rungeT(h,n):
    y1.clear()
    y2.clear()
    def F(y1,y2):
        return -y2+y1*(y1**2+y2**2-1), y1+y2*(y1**2+y2**2-1)

    r=0
    y1.append(np.cos(xLeft)/(np.sqrt(1+np.exp(2*xLeft))))
    y2.append(np.sin(xLeft)/(np.sqrt(1+np.exp(2*xLeft))))
    for i in range(1,n):
        k11,k21 = F(y1[r],y2[r])
        k12, k22 = F(y1[r] + h*k11/2,y2[r] + h*k21/2)
        k13,k23 = F(y1[r] + h*k12/2,y2[r] + h*k22/2)
        k14, k24 = F(y1[r] + h*k13,y2[r] + h*k23)
    
        y1.append ( y1[r] + h*(k11 + 2*k12 + 2*k13 + k14)/6)
        y2.append ( y2[r] + h*(k21 + 2*k22 + 2*k23 + k24)/6)
        r+=1

def runge(h, n, Da):
    y1.clear()
    y2.clear()
    def F(x,Teta, Da):
            return -Alfa*x+Da*(1-x)*np.exp(Teta/(1+Teta/Gamma)),-Alfa*Teta+Da*B*(1-x)*np.exp(Teta/(1+Teta/Gamma)) - Beta*(Teta - Teta_C)

    r=0
    y1.append(x0)
    y2.append(Teta0)
    for i in range(1,n):
        k11,k21 = F(y1[r],y2[r], Da)
        k12, k22 = F(y1[r] + h*k11/2,y2[r] + h*k21/2, Da)
        k13,k23 = F(y1[r] + h*k12/2,y2[r] + h*k22/2, Da)
        k14, k24 = F(y1[r] + h*k13,y2[r] + h*k23, Da)
    
        y1.append ( y1[r] + h*(k11 + 2*k12 + 2*k13 + k14)/6)
        y2.append ( y2[r] + h*(k21 + 2*k22 + 2*k23 + k24)/6)
        r+=1

plt.figure(figsize=(6,6))
def CalculZ(n,x):
    z.clear()
    z2.clear()
    for i in range(0,n): # Вывод функции к которой стремимся в тестовой части
        z.append(np.cos(x)/(np.sqrt(1+np.exp(2*x))))
        z2.append(np.sin(x)/(np.sqrt(1+np.exp(2*x))))
        x+=h

CalculZ(n,x)
plt.subplot(211)
plt.plot(t,z, linewidth=4, label="exact solution y1")
rungeT(h,n)
plt.plot(t,list(y1), linewidth=2, label="found solution y1")
plt.legend()
plt.subplot(212)
plt.plot(t, z2, linewidth=4, label="exact solution y2")
plt.plot(t,list(y2), linewidth=2, label="found solution y2")
plt.legend()
plt.show()

ht = 0.1
es1 = list()
es2 =list()
hts = np.arange(h,ht,h)
hts = hts[::-1]
print(hts)
while ht >= h:
    nt=(int)((xRight-xLeft)/ht)
    CalculZ(nt,x)
    rungeT(ht,nt)
    e1 = np.max(np.abs(np.subtract(z, y1)))
    e2 = np.max(np.abs(np.subtract(y2, z2)))
    # print(e1, "  -  ", e2, "  :  ", ht)
    es1.append(e1)
    es2.append(e2)

    ht -= h
    
plt.subplot(211)
plt.xlim(0.1,0.0005)
plt.plot(hts,es1, marker = '.', label="error 1")
plt.legend()
plt.subplot(212)
plt.xlim(0.1,0.0005)
plt.plot(hts,es2, marker = '.', label="error 2")
plt.legend()
plt.show()
    

for i in range(1, len(Da)+1):
    runge(h, n , Da[i-1])
    plt.subplot(len(Da)//2+1, 2, i)
    plt.plot(t,list(y1), linewidth=2, label = 'Y1', color = 'red')

    # plt.subplot(len(Da)//2, 2, i)
    plt.plot(t,list(y2), linewidth=2, label= 'Y2', color = 'blue')
    plt.legend()
    plt.grid()


plt.show()
