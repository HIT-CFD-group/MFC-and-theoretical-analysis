# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:32:36 2023

@author: Administrator
"""

import matplotlib.pyplot as plt
import math
import numpy as np
#parameter calculation
Ma=2.4
p0a=101325
rho_g=1.204
rho_l=657#657#1000
gam_g = 1.4
rho = 1
gam_l = 5.394#5.394#6.12
pi_l = 1.45E8#1.45E8#3.43E8
c_a = math.sqrt( gam_g * ( p0a  ) / rho_g )
psOp0a = ( Ma ** 2 -1 ) * 2 * gam_g / ( gam_g + 1 ) + 1
rhosOrho0a = ( 1 + ( gam_g + 1 ) / ( gam_g - 1) * psOp0a ) / ( ( gam_g + 1 ) / ( gam_g - 1) + psOp0a )
Ms = math.sqrt( ( gam_g + 1. ) / ( 2. * gam_g ) * ( psOp0a - 1. ) * ( p0a / ( p0a ) ) + 1. )
vel = c_a/gam_g * (psOp0a - 1.) * p0a / ( p0a  ) / Ms
MPD=100 #mesh per diameter
Ny = 12*MPD/2
Nx = 14*MPD
D=0.022
dx = D/MPD 
ps=p0a*psOp0a
rho_post_g=rhosOrho0a*rho_g
cfl = 0.2
c_dt = math.sqrt( 1.4*ps/rho )
dt = cfl * dx/c_dt
print('dt',dt)
print('Ma',Ma)
#c_l =1000.196537
c_l=12.064*Ma+1045
######################################
#c_l=1060.96
#1129.626
c_air=343
r=0.011
pi=np.pi
'''
t_start=1#0.228 beginning of Mach stem
t=t_start*2*r/c_l
'''
timestep=1173
t=timestep*dt
n=c_l/(Ma*c_air)
#critical angle when shock deattached the droplet
alpha_critical=np.arcsin(1/n)
print('critical',alpha_critical*180/pi)
print(n)
#location
if n <=1 :
    critical = 180
else:
    critical= alpha_critical*180/pi
alpha=np.linspace(0, critical,100)
#alpha=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28])
#refraction angle
theta=np.arcsin(n*np.sin(alpha*pi/180))
theta=np.degrees(theta)

gamma_1=2*theta-alpha #different reflection point
gamma_2=180+4*theta-alpha
x_0=r*np.cos(pi-alpha*pi/180)
y_0=r*np.sin(pi-alpha*pi/180)
x_1=r*np.cos(gamma_1*pi/180)
y_1=r*np.sin(gamma_1*pi/180)
x_2=r*np.cos(gamma_2*pi/180)
y_2=r*np.sin(gamma_2*pi/180)
x_inside_01=np.linspace(x_0,x_1,100)
y_inside_01=np.linspace(y_0,y_1,100)
x_inside_12=np.linspace(x_1,x_2,200)#first dimension indicates the point index between two reflection point
                                 #second dimension indicates different incident angle
y_inside_12=np.linspace(y_1,y_2,200)
Mach_stem=np.array([[250,-0.00164882, 0.011421],[300,0.000549607, 0.011421],
                   [350,0.00252819, 0.0112013],[400,0.00428694, 0.0105424],
                   [450,0.005606, 0.00988353],[500,0.00670521, 0.00900499],
                   [550,0.00780443, 0.00812646],[600,0.0086838, 0.00724792],
                   [650,0.00956317, 0.00636938],[700,0.0100029, 0.00527121],
                   [750,0.0104425, 0.00439268],[800,0.0108822, 0.00329451],
                   [850,0.0111021, 0.00241597],[900,0.0113219, 0.00153744],
                   [950,0.0113219, 0.000658902],[980,0.011,0]])

for i in range(len(Mach_stem)):
    Mach_stem[i][0]=Mach_stem[i][0]*dt
    
def yi_yuan_er_ci_equation(a,b,c):
    d=b**2-4*a*c
    if (d<0):
        print("无解")
    else:
        e = math.sqrt(d)
        x1=((-b+e)/(2*a))#调用math模块中sqrt开平方函数
        x2=((-b-e)/(2*a))
    return x1,x2


def get_beginning_time_by_Mach_stem_velocity(Mach_stem,x,y,alpha):#time used for starting disturbation from t=0
    #improved model by mach stem motion
    '''
    velocity=np.zeros(len(Mach_stem)-1)
    for i in range(0,len(Mach_stem)-1):
        velocity[i]=(np.arccos(-Mach_stem[i+1][1]/np.sqrt(Mach_stem[i+1][1]**2+Mach_stem[i+1][2]**2))-
                    np.arccos(-Mach_stem[i][1]/np.sqrt(Mach_stem[i][1]**2+Mach_stem[i][2]**2)))/(Mach_stem[i+1][0]-Mach_stem[i][0])
        #velocity[i]=np.sqrt((Mach_stem[i+1][2]-Mach_stem[i][2])**2+     \
                         #(Mach_stem[i+1][1]-Mach_stem[i][1])**2)/(Mach_stem[i+1][0]-Mach_stem[i][0])
    '''
    #improved model with fitted line
    #time=r*(1-np.cos(alpha*pi/180))/(Ma*c_air)#location without mach stem effect
    if alpha*pi/180 <=np.arccos(-Mach_stem[0][1]/np.sqrt(Mach_stem[0][1]**2+Mach_stem[0][2]**2)):
        time=r*(1-np.cos(alpha*pi/180))/(Ma*c_air)#location without mach stem effect
    

    else:
        dt=(1.8098e-7/math.sqrt(2.8*Ma*Ma-0.4))
        a=(-8e-5)/dt/dt
        b=0.2334/dt
        c=30.531
        # a=(-9e-5)/dt/dt
        # b=0.2442/dt
        # c=27.132
        time=yi_yuan_er_ci_equation(a,b,c-alpha)[0]
    

    '''
    else:
        for i in range(0,len(Mach_stem)-1):#location under mach stem effect
            if alpha*pi/180 >= np.arccos(-Mach_stem[i][1]/np.sqrt(Mach_stem[i][1]**2+Mach_stem[i][2]**2)) \
                and alpha*pi/180 <= np.arccos(-Mach_stem[i+1][1]/np.sqrt(Mach_stem[i+1][1]**2+Mach_stem[i+1][2]**2)):
                    time=Mach_stem[i][0]+(alpha*pi/180-np.arccos(-Mach_stem[i][1]/np.sqrt(Mach_stem[i][1]**2+Mach_stem[i][2]**2)))/velocity[i]
                    #time=Mach_stem[i][0]+np.sqrt((x- Mach_stem[i][1])**2+(y- Mach_stem[i][2])**2)/velocity[i]
    #'''
    #original model
    #time=r*(1-np.cos(alpha*pi/180))/(Ma*c_air)
    
        
    return time
               
def intersection_point(x_icw,y_icw,x1,y1,x2,y2):
                                   #x1,y1 the disturbation of reflect wave
                                   #x2,y2 the droplet shape
    distance_min=1
    if x_icw <=0:
        for i in range(len(x1)):
            for j in range(len(x2)):
                if x1[i]>=x_icw and y1[i]>=y_icw:#First quadrant is only considered before 
                                                     #incident shock wave move across the half of the droplet
                    distance=(y2[j]-y1[i])**2+(x2[j]-x1[i])**2
                    if distance<=distance_min:
                          distance_min=distance
                          index_i=i
                          index_j=j
    else:
        for i in range(len(x1)):
            for j in range(len(x2)):
                if x1[i]>=x_icw and y1[i]<=y_icw:#Forth quadrant is only considered after 
                                                     #incident shock wave move across the half of the droplet
                    distance=(y2[j]-y1[i])**2+(x2[j]-x1[i])**2
                    if distance<=distance_min:
                          distance_min=distance
                          index_i=i
                          index_j=j
        '''
            elif distance > 2*r:
                index_point==2
                distance_min=1
            elif distance<=distance_min and index_point==2:
                print('AAA')
                distance_min=distance
                index_i_2=i
                index_j_2=j
    if x1[index_i_1]<x1[index_i_2]:
        index_i=index_i_2
        index_j=index_j_2
    else:
        index_i=index_i_1
        index_j=index_j_1
        '''
    return x1[index_i],y1[index_i]

def plot_circle(x_o,y_o,r):
    theta=np.arange(0,2*pi,0.01)
    x=x_o+r*np.cos(theta)
    y=y_o+r*np.sin(theta)
    plt.plot(x,y,color='black')
    return x,y


def detect_if_disturb_in_circle(origin,x_end,y_end,r_circle):
    included_alpha_inside=[]
    if origin == 100:#origin is the end of reflective disturbance
        origin = x_end[0]
    for i in range(len(x_end)-1):
        if (x_end[i+1]-origin)**2+(y_end[i+1]-y_end[0])**2<r_circle**2:
           included_alpha_inside.append(alpha[i])
    return(len(included_alpha_inside)/len(alpha))

def main(timestep,mode):
    global x_trasm,y_trasm,x_icw,y_icw,x_end,y_end
    t=timestep*dt
    x_end=[]
    y_end=[]
    for i in range(len(x_1)): #i represents different line
        beginning_time=get_beginning_time_by_Mach_stem_velocity(Mach_stem, x_0[i], y_0[i], alpha[i])
        l_pm_0=c_l*(t-beginning_time)
        l_pm_1=c_l*(t-beginning_time)-2*r*np.cos(theta[i]*pi/180)
        x_trasm=[]
        y_trasm=[]
        x_trasm.append(x_0[i])
        y_trasm.append(y_0[i])
        #if c_l*t>=n*r*(1-np.cos(alpha[i]*pi/180)): #rays begin to emit
        if t-beginning_time>=0:
            x_icw=x_0[i]#######
            y_icw=y_0[i]
            #x_rw,y_rw=reflect_wave(t_rw, x_icw, y_icw)
            #y_rw=reflect_wave(t_rw, x_icw, y_icw)[1]
            #Mach stem effect
            #x_interintersection_point=intersection_point(x_icw,y_icw,x_rw, y_rw, x_droplet, y_droplet)[0]
            #y_interintersection_point=intersection_point(x_icw,y_icw,x_rw, y_rw, x_droplet, y_droplet)[1]
            #if c_l*t-n*r*(1-np.cos(alpha[i]*pi/180)) >= 2*r*np.cos(theta[i]*pi/180): #reflect more than once
            if c_l*(t-beginning_time)-2*r*np.cos(theta[i]*pi/180) >=0 :
                x_trasm.append(x_1[i])
                y_trasm.append(y_1[i])
                for j in range(len(x_inside_12)):
                    if (x_inside_12[j][i]-x_1[i])**2+(y_inside_12[j][i]-y_1[i])**2 < l_pm_1**2:
                        x_trasm.append(x_inside_12[j][i])
                        y_trasm.append(y_inside_12[j][i])
                
                x_end.append(x_trasm[-1])
                y_end.append(y_trasm[-1])
            #'''
            else:   #no reflection
                for j in range(len(x_inside_01)):
                    if(x_inside_01[j][i]-x_0[i])**2+(y_inside_01[j][i]-y_0[i])**2 < l_pm_0**2:
                        x_trasm.append(x_inside_01[j][i])
                        y_trasm.append(y_inside_01[j][i])
            #'''
        #x_reflect.append(x_2[i])
        #y_reflect.append(y_2[i])
        #x_end.append(x_trasm[-1]) #ending line including one or zero time reflection
        #y_end.append(y_trasm[-1])
        y_trasm_2=np.array(y_trasm)*(-1)
        #plt.plot(x_trasm,y_trasm_2)
        if mode==0:     
            plot(x_trasm,y_trasm,x_end,y_end)
        
        #end line distrubution
    #print(x_end)
    x_icw_line=[-0.011+Ma*c_air*t,-0.011+Ma*c_air*t]
    y_icw_line=[y_icw-10,y_icw+10]
    #incident shock wave
    #plt.plot(x_icw_line,y_icw_line)
    #for i in range(10,35,5):
        #print(i,detect_if_disturb_in_circle(D/8,x_end, y_end, D/i))

def plot(x_trasm,y_trasm,x_end,y_end):
    x_droplet=plot_circle(0,0,0.011)[0]
    y_droplet=plot_circle(0,0,0.011)[1]
    plt.plot(x_trasm,y_trasm,linewidth=1)
    #incident shock wave's position

    #symmetry
    y_end_2=np.array(y_end)
    y_end_2=y_end_2*(-1)
    
    #Prediction points
    loc=np.array([0.333,0.296225918,0.2448,0.191])
    x_point=np.array([0.00976439])
    y_point=np.zeros(len(x_point))  
    
    plt.plot(x_end, y_end,label=t,color='blue')   
    plt.plot(x_end, y_end_2,color='blue')
    #plt.scatter(x_point, y_point,label='point',s=20,color='red',zorder=2)  
    plt.axis('equal')
    plt.xlim(-0.022,0.022)
    plt.ylim(-0.022,0.022)
    #plot_circle(0.00976439, 0, D/20)
    #plot_circle(0.00976439, 0, D/40)
    
def detect():
    timestep=np.linspace(D/c_l/dt,2*D/c_l/dt,100)
    #sth=[D/4,0.00681513,0.00976439]
    #Ms2.4
    #sth=[-D/4,-D/8,0,D/8,D/4,0.00571592]
    #Ms4.2
    #sth=[-D/4,-D/8,0,D/8,D/4,0.0072548,0.0075885]
    #Ms4.8
    #sth=[-D/4,-D/8,0,D/8,D/4,0.0072548,0.0085242]
    #Ms5.4
    #sth=[-D/4,-D/8,0,D/8,D/4,0.00681513,0.00976439]
    #sth=[0.00976439]
    sth=[D/4]
    k=len(sth)
    j=len(timestep)
    included_alpha=np.arange(k*j).reshape((k, j))
    for k in range(len(sth)):
        #for i in range(20,45,5):
            for j in range(len(timestep)):#time iteration
                main(timestep[j],1)
                print(detect_if_disturb_in_circle(sth[k],x_end, y_end, D/40))
                included_alpha[k][j]=detect_if_disturb_in_circle(sth[k],x_end, y_end, D/40)
                #included_alpha[k][j]=(detect_if_disturb_in_circle(sth[k],x_end, y_end, D/40))
            '''
            main((0.033-k)/dt/c_l,1)
            included_alpha.append(detect_if_disturb_in_circle(k,x_end, y_end, D/i))
            '''
            #included_alpha=np.array(included_alpha)
            #print('timestep of',i,timestep)
            #print('inc_alpha of',i,included_alpha)
    print(included_alpha)
    print('max',max(included_alpha))
            #print('timestep',(0.033-k)/dt/c_l)
    index = np.unravel_index(included_alpha.argmax(),included_alpha.shape)
    print('timestep when max',timestep[index[1]])
    print('location of focus',timestep[index[0]])
    #print('timestep when max',i,timestep[np.argmax(included_alpha)])
    #print('location of focus',0.033-timestep[np.argmax(included_alpha)]*dt*c_l)
            #print('min',min(included_alpha))
            #print('timestep when min',i,timestep[np.argmin(included_alpha)])
            #print(0.033-timestep[np.argmax(included_alpha)]*dt*c_l)
    #plt.plot(timestep,included_alpha)
        
    # x_line=[2.875E-05/dt,2.875E-05/dt]
    # y_line=[0,1.8*max(included_alpha)]
    # plt.plot(x_line,y_line,linewidth=1,c='black')
    # x_line=[3.012E-05/dt,3.012E-05/dt]
    # y_line=[0,1.8*max(included_alpha)]
    # plt.plot(x_line,y_line,linewidth=1,c='black')
    #plt.xlabel('timestep')
    #plt.ylabel('included_alpha_percent') 
    #plt.legend()

detect()
#main((0.033-0.007588491)/dt/c_l,0)
#main(0.022/dt/c_l,0)