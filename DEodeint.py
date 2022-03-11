#无阻尼是否会稳定呢？2
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def diff_equation(y_list,x):
    y,w=y_list
    return np.array([(w-18000)*np.pi/180,((18000/8.18)*(1-1.35*np.sin(y)))])
x=np.linspace(0,15,num=100000)
#y0=[60.1090*np.pi/180,18280.35]
y0=[34.53*np.pi/180,18000]
result=odeint(diff_equation,y0,x)
plt.plot(x,result[:,0]*180/np.pi,label='y')
#plt.plot(x,result[:,1],label='w')
plt.legend()
plt.grid()
plt.show()