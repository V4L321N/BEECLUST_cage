import numpy as np
from scipy import signal
from scipy import integrate
from scipy.integrate import odeint
import matplotlib.pyplot as plt

tmin = 0
tmax = 105
t = np.linspace(tmin, tmax, tmax, endpoint=False)
print(t)

def PULSE(start, end, height, t):
    return np.heaviside(t - start, 0) * height - np.heaviside(t - end, 1) * height

arena_area = 30 ** 2 * np.pi                      #in cm^2
bee_detection_area = 6                      #in cm^2
#x1=bee_detection_area /(arena_area * 0.5)   # probability of a free bee to meet another free bee in a second
x1 = 0.1
x2 = 0.85                                     # bees in cage stop a moving bee more or less likely (due to cage size?)
Cr = 6                                        # number of bees in the left cage
Cl = 0                                        # number of bees in the right cage
L0 = 0                                        # Initial number of bees aggregated on left side
F0 = 24                                       # Initial number of bees running freely
R0 = 0                                        # Initial number of bees aggregated on right side
# tmin = 0                                      # start with this value of the independent variable (t)
# tmax = 100                                    # run for this number of time steps

#two functions for the temperature profile
def Tl(t):
    return 32 + PULSE(30, tmax + 1, 4, t) #step function for left side. "height" argument is at 4 which raises temperatre to 36

def Tr(t):
    return 32 + PULSE(30, tmax + 1, 0, t) #step function for right side. "height" argument is at 0.

# plt.title('temperature_BEECLUST_cage')
# plt.plot(t, Tl(t), label = 'Temp left', linestyle = '--', color = 'blue')
# plt.plot(t, Tr(t), label = 'Temp right', linestyle = '-.', color = 'red')
# plt.xlabel('time [min]')
# plt.ylabel('Temperature [°C]')
# plt.show()

#two functions that convert temperatures into waiting times
a = -2.28
b = 0.08
c = 5.15
d = 0.32
e = 28.5
f = 7

def Wl(t):
    return (((a + b * Tl(t)) ** c / ((a + b * Tl(t)) ** c + d ** c)) * e + f) / 15

def Wr(t):
    return (((a + b * Tr(t)) ** c / ((a + b * Tr(t)) ** c + d ** c)) * e + f) / 15

# plt.plot(t, Wl(t))
# plt.plot(t, Wr(t))
# plt.show()

Temp = np.linspace(28, 40, 100, endpoint=False)

def Wartezeit(Temp):
    return (((a + b * Temp) ** c / ((a + b * Temp) ** c + d ** c)) * e + f) / 15

# plt.plot(Temp, Wartezeit(Temp))
# plt.show()
initial =  [L0, F0, R0]

def Diff_EQ(t, y):
    L, F, R = y
    deL = x1 * (F ** 2) + x1 * F * (L + x2 * Cl) - L / (Wl(t))
    deF = L / (Wl(t)) + R / (Wr(t)) - x1 * (F ** 2) - x1 * F * (L + x2 * Cl) - x1 * (F ** 2) - x1 * F * (R + x2 * Cr)
    deR = x1 * (F ** 2) + x1 * F * (R + x2 * Cr) - R / (Wr(t))
    # deL = np.diff(L,t) == x1 * (F ** 2) + x1 * F * (L + x2 * Cl) - L / (Wl(t))
    # deF = np.diff(F,t) == L / (Wl(t)) + R / (Wr(t)) - x1 * (F ** 2) - x1 * F * (L + x2 * Cl) - x1 * (F ** 2) - x1 * F * (R + x2 * Cr)
    # deR = np.diff(R,t) == x1 * (F ** 2) + x1 * F * (R + x2 * Cr) - R / (Wr(t))
    return [deL, deF, deR]

# def deL(t, L):
#     return np.diff(np.array([x1 * (F ** 2) + x1 * F * (L + x2 * Cl) - L / (Wl(t))]))
# #
# def deF(t, F):
#     return np.diff(np.array([L / (Wl(t)) + R / (Wr(t)) - x1 * (F ** 2) - x1 * F * (L + x2 * Cl) - x1 * (F ** 2) - x1 * F * (R + x2 * Cr)]))
# #
# def deR(t, R):
#     return np.diff(np.array([x1 * (F ** 2) + x1 * F * (R + x2 * Cr) - R / (Wr(t))]))


P = integrate.solve_ivp(Diff_EQ, [tmin, tmax], initial, method='RK45')
# P = integrate.RK45(Diff_EQ, tmin, initial, tmax, 0.01)
print(P.t, P.y[0], P.y[1], P.y[2])
print(P)


# Q1=[[t,L] for t,L,F,R in P]
# LP1=list_plot(Q1, color='red', plotjoined=True, legend_label='Left aggregated bees', linestyle='-.' )
#
# Q2=[[t,F] for t,L,F,R in P]
# LP2=list_plot(Q2, color='green', plotjoined=True, legend_label='Free running bees', linestyle='--')
#
# Q3=[[t,R] for t,L,F,R in P]
# LP3=list_plot(Q3, color='blue', plotjoined=True, legend_label='Right aggregated bees', linestyle=':')
#
# LP = LP1+LP2+LP3
# LP.axes_labels(['time [min]','bees'])
# #LP.axes_width(2)
# show(LP)
plt.plot(P.t, P.y[0])
plt.plot(P.t, P.y[1])
plt.plot(P.t, P.y[2])
plt.show()
