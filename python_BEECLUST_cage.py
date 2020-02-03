import numpy as np
from scipy import signal
from scipy import integrate
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.stats import gompertz

tmin = 0
tmax = 105
t = np.linspace(tmin, tmax, tmax, endpoint=False)
print(t)

def PULSE(start, end, height, t):
    return np.heaviside(t - start, 0) * height - np.heaviside(t - end, 1) * height

# def tanh_DECREASE(start, end, height, t):
#     delta = end - start
#     return 36 - 6 * (np.sin(t/(11) + 47.39) ** 2) * np.heaviside(t - start, 0) * PULSE(start, end, 1, t) - 6 * PULSE(49, tmax + 1, 1, t)

# The following function is actually a pretty dirty hack, but seems to be working.
# The first term(36 in this case) describes the starting temperature, the second term(6 * (np.sin(t/(11.00) + 47.39) ** 2) * PULSE(start, end, 1, t))
# defines the descent via a sinus function that has been stretched to fit into the needed 20 minutes, starting at roughly t=30 until t=50.
# The third term comes into play from there on, holding the temperature at the point where the sinus term left it.
def tanh_DECREASE(start, end, height, t):
    return 36 - 6 * (np.sin(t/(11.00) + 47.39) ** 2) * PULSE(start, end, 1, t) - 6 * PULSE(49.99, tmax + 1, 1, t)


arena_area = 30 ** 2 * np.pi                      #in cm^2
bee_detection_area = 6                      #in cm^2
#x1=bee_detection_area /(arena_area * 0.5)   # probability of a free bee to meet another free bee in a second
x1 = 0.1
x2 = 0.85                                     # bees in cage stop a moving bee more or less likely (due to cage size?)
Cr = 0                                        # number of bees in the left cage
Cl = 0                                        # number of bees in the right cage
L0 = 0                                        # Initial number of bees aggregated on left side
F0 = 64                                       # Initial number of bees running freely
R0 = 0                                        # Initial number of bees aggregated on right side
# tmin = 0                                      # start with this value of the independent variable (t)
# tmax = 105                                    # run for this number of time steps # different at caged experiments (30?) !!!

#two functions for the temperature profile

# def Tl(t):
#     return 36 - PULSE(30, tmax + 1, 6, t) #step function for left side. "height" argument is at 6 which raises temperatre to 36

def Tl(t):
    return tanh_DECREASE(30.00, 49.98, 6, t)

def Tr(t):
    return 32 + PULSE(30, 50, 0, t) #step function for right side. starting time is 30 ; "height" argument is at 0.

plt.title('temperature_BEECLUST_off')
plt.plot(t, Tl(t), label = 'temp. left', linestyle = '--', color = 'blue')
plt.plot(t, Tr(t), label = 'temp. right', linestyle = '-.', color = 'red')
plt.xlabel('time [min]')
plt.ylabel('temperature [°C]')
plt.show()

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

# def Wartezeit(Temp):
#     return (((a + b * Temp) ** c / ((a + b * Temp) ** c + d ** c)) * e + f) / 15

# plt.plot(Temp, Wartezeit(Temp))
# plt.xlabel('temperature / °Celsius', fontsize=18)
# plt.ylabel('waiting time / minutes', fontsize=16)
# plt.show()

initial =  [L0, F0, R0]

def Diff_EQ(t, y):
    L, F, R = y
    deL = x1 * (F ** 2) + x1 * F * (L + x2 * Cl) - L / (Wl(t))
    deF = L / (Wl(t)) + R / (Wr(t)) - x1 * (F ** 2) - x1 * F * (L + x2 * Cl) - x1 * (F ** 2) - x1 * F * (R + x2 * Cr)
    deR = x1 * (F ** 2) + x1 * F * (R + x2 * Cr) - R / (Wr(t))
    return [deL, deF, deR]

P = integrate.solve_ivp(Diff_EQ, [tmin, tmax], initial, method='RK45')
# print(P.t, P.y[0], P.y[1], P.y[2])

plt.plot(P.t, (P.y[0]/F0) * 100, label = "left")
plt.plot(P.t, (P.y[1]/F0) * 100, label = "free")
plt.plot(P.t, (P.y[2]/F0) * 100, label = "right")
plt.legend(loc="upper right")
plt.xlabel('time / minutes', fontsize=18)
plt.ylabel('number / bees', fontsize=16)
# plt.plot(P.t, P.y[0], P.y[1], P.y[2])
plt.show()
