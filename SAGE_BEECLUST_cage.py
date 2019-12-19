forget()

#import the runge-kutta solver from library
from sage.calculus.desolvers import desolve_system_rk4

#declare the used variables
L,F,R,t = var('L F R t')
Tl, Tr, Wl, Wr = var('Tl Tr Wl Wr')
Wartezeit = var ('Wartezeit')

#deL, deR, deF = var('deL deR deF')
#x1,x2,Cl,Cr = var('x1 x2 Cl Cr')
#L0,F0,R0 = var('L0 F0 R0')
#a,b,c,d,e,f = var('a b c d e f')
#tmin,tmax = var('tmin tmax')
#arena_area, bee_detection_area = var('arena_area bee_detection_area')

PULSE(start, end, height, t) =   unit_step(t - start) * height - unit_step(t - end) * height

#inititalize the values
arena_area = 30^2*pi()                      #in cm^2
bee_detection_area = 6                      #in cm^2
#x1=bee_detection_area /(arena_area * 0.5)   # probability of a free bee to meet another free bee in a second 
x1=0.1
x2=0.85                                      # bees in cage stop a moving bee more or less likely (due to cage size?)    
Cr=6                                        # number of bees in the left cage
Cl=0                                        # number of bees in the right cage
L0=0                                        # Initial number of bees aggregated on left side
F0=24                                       # Initial number of bees running freely
R0=0                                        # Initial number of bees aggregated on right side
tmin=0                                      # start with this value of the independent variable (t)
tmax=105                                    # run for this number of time steps

#two functions for the temperature profile      
Tl(t) = 32 + PULSE(30, tmax+1, 4, t)   
Tr(t) = 32

#plot temperature profile
TP1=plot(Tl(t), (0,tmax), color='red', legend_label='Temperature left', linestyle='--')
TP2=plot(Tr(t), (0,tmax), color='blue', legend_label='Temperature right', linestyle='-.')
TP=TP1+TP2
TP.axes_labels(['time [min]','temperature [Celsius]'])
show(TP)

#two functions that convert temperatures into waiting times
a = -2.28
b = 0.08
c = 5.15
d = 0.32
e = 28.5
f = 7
Wl(t)= (((a+b*Tl(t))^c/((a+b*Tl(t))^c+d^c))*e+f) / 15
Wr(t)= (((a+b*Tr(t))^c/((a+b*Tr(t))^c+d^c))*e+f) /15

Wartezeit(Temp) = (((a+b*Temp)^c/((a+b*Temp)^c+d^c))*e+f) / 15
WPx=plot(Wartezeit(Temp), (28,40), color='red', legend_label='waiting time', linestyle='--')
WPx.axes_labels(['temperature (C)','waiting time of bees [min]'])
show(WPx)


#plot waiting-time profile
WP1=plot(Wl(t), (0,tmax), color='red', legend_label='waiting time left', linestyle='--')
WP2=plot(Wr(t), (0,tmax), color='blue', legend_label='waiting time right', linestyle='-.')
WP=WP1+WP2
WP.axes_labels(['time [min]','waiting time of bees [min]'])
show(WP)

#here are the main model equations describing the bees
deL = diff(L,t) == x1*(F^2) + x1*F*(L+x2*Cl) - L/(Wl(t))
deF = diff(F,t) == L/(Wl(t)) + R/(Wr(t)) - x1*(F^2) - x1*F*(L+x2*Cl) - x1*(F^2) - x1*F*(R+x2*Cr)
deR = diff(R,t) == x1*(F^2) + x1*F*(R+x2*Cr) - R/(Wr(t))

#solve them with runge kutta
P=desolve_system_rk4([
 deL.rhs(), 
 deF.rhs(), 
 deR.rhs()
 ],[L,F,R],ics=[tmin,L0,F0,R0],ivar=t,end_points=tmax )

#plot them 
Q1=[[t,L] for t,L,F,R in P]
LP1=list_plot(Q1, color='red', plotjoined=True, legend_label='Left aggregated bees', linestyle='-.' )

Q2=[[t,F] for t,L,F,R in P]
LP2=list_plot(Q2, color='green', plotjoined=True, legend_label='Free running bees', linestyle='--')

Q3=[[t,R] for t,L,F,R in P]
LP3=list_plot(Q3, color='blue', plotjoined=True, legend_label='Right aggregated bees', linestyle=':')

LP = LP1+LP2+LP3
LP.axes_labels(['time [min]','bees'])
#LP.axes_width(2)
show(LP)


#print x1

#solve([32==F-((x1*F*F)/(x1*F-1/34.5))-((x1*F*F)/(x1*F-1/16.5))],F)
