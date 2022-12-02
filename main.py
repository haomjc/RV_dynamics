import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_bvp, solve_ivp
from rv_dyna import rv_dyna

x_fig= []
y_fig = []

#x0=[0, 0, 0, 0, 0, 0,-0.01, -0.1, -0.01, -0.1, -0.01, -0.1];
x0=[0, 0, 0, 0, 0, 0, -0, -0, -0, -0, -0, -0];
T=2*np.pi
omega_hx = 2

T_all = 60*T
tspan=np.arange(0,T_all,T/10000)

y =odeint(rv_dyna,x0,tspan, (omega_hx, T_all));

#y=solve_ivp(rv_dyna,tspan,x0,  args = (omega_hx[i],),method='RK23');

#plt.plot(np.full((300, 1), omega_hx[i]),y[300000:y.shape[0]-5:1000])

#sc.set_data(np.c_[x_fig,y_fig])

#b_d = 10**(-6)*10**(4)
b_d = 10**(-4)
#omega_d = 107029.37
#omega_d = 4.7671*10**2
m = 1.5*10**(-3)
z_1 = 15
z_2  = 40

r1=m*z_1/2;#分度圆半径
r2=m*z_2/2;

z_g = 39
z_b = z_g+1

i_16_5 = z_b*z_2/z_1+1

theta_0 = omega_hx*tspan

x_01 = y[:, 0]*(b_d)
x_12 = y[:, 1]*(b_d)
x_23 = y[:, 3]*(b_d)
x_25 = y[:, 5]*(b_d)



theta_1 = theta_0 - x_01

theta_2 = (r1*theta_1 - x_12)/r2
theta_3 = (theta_2 - x_23)/(i_16_5/(z_2/z_1))
theta_5 =  (theta_2 - x_25)/(i_16_5/(z_2/z_1))   #真实的输出角度
theta_R = omega_hx*tspan/i_16_5   #理论输出角度

TR =  theta_R - theta_5

#plt.plot(tspan, theta_0/theta_5)
plt.plot(tspan, TR)
#plt.plot(tspan, theta_2/theta_3)
#plt.plot(tspan, y[:, 0]*(b_d))
#plt.plot(tspan, y)
plt.xlabel('t(s)');plt.ylabel('p');
plt.show()
