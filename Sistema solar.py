# -*- coding: cp1252 -*-
# 20 de Noviembre 2015
#Movimiento tres cuerpos:generalización a n cuerpos mediante la implementacion de vectores en python
#En esta version se simulara el comportamiento en el sistema solar.

import math
from math import *
from visual import *

##############################################################################################################################################################
##############################################################################################################################################################

#SECCION 1:
#Definicion de funciones:

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#   1ra funcion: Defino segunda derivada de la (dv/dt), ecuación de movimiento

def fv1(t,x,v,m,n):             #Funcion que encuentra la velocidad final de cada cuerpo, luego de transcurrido un tiempo, en donde x es un vector,
                                #que tiene en sus components todas las posiciones de los diferentes objetos, v las velocidades iniciales respectivas a cada objeto
                                #m las masas y n la cantidad de objetos
    vec=vector(0,0,0)
    dxi=[vec]*n;
    G=-6.673884e-11
    d=5e7                       #radio minimo en el cual la velocidad no aumentará mas para evitar divergencia por division por cero o numeros muy pequeños

    for i in range(0,n):
        for j in range(0,n):
            if j!=i:            #evita restar el mismo vector
                rij=sqrt(dot(x[i]-x[j],x[i]-x[j]))
                if rij>d:       #Si la distancia entre los dos cuerpos es menor a d, el incremento de la velocidad sera cero.
                    dxi[i]=dxi[i]+G*m[j]*(x[i]-x[j])/pow(rij,3)
                else :
                    dxi[i]=dxi[i]
        
    return dxi
#----------------------------------------------------------------------------------------------------------------------------------------------------------

#   2da funcion: Defino primera derivada  de la (dx/dt), ecuación de movimiento

def fx(v):                    #Función que da la Posición del sistema
    return v
#-----------------------------------------------------------------------------------------------------------------------------------------------------------               

#   3ra función: Operacion necesaria para hallar k2,k3 y k4 del metodo Ruge-Kutta

def operacion(kix,kiv,xio,vio,n,var,h):               #Función que me realiza las operaciones del metodo ruge-kutta, ya que las entradas son listas y no puedo manipularlas diectamente
    vec=vector(0,0,0)
    xip=[vec]*n;
    vip=[vec]*n;
    kipx=[vec]*n;
    kipv=[vec]*n;
    for i in range (0,n):
        kipx[i]=h*kix[i]
        kipv[i]=h*kiv[i]
        if var==1:                                  #Operación realizada para hallar k2 y k3
            xip[i]=xio[i]+kipx[i]/2
            vip[i]=vio[i]+kipv[i]/2
        else:                                       #Operacioón utilizada para hallar k4
            xip[i]=xio[i]+kipx[i]
            vip[i]=vio[i]+kipv[i]

    return xip,vip,kipx,kipv
#------------------------------------------------------------------------------------------------------------------------------------------------------------

#   4ta función: Operacion necesaria para hallar los nuevos valores de velocidad y posición (vf y xf) en el metodo Ruge-kutta
def operacion2(k1,k2,k3,k4,xio2,n):
    vec=vector(0,0,0)
    xip2=[vec]*n;
        
    for i in range(0,n):
        xip2[i]=xio2[i]+(k1[i]+2*(k2[i]+k3[i])+k4[i])/6.
    
    return xip2
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

#############################################################################################################################################################

# METODO RUGE-KUTTA

def rk4(fv,fx,xi,vi,tf,ti,m,n):         #xi, vi vectores con todas las posiciones y velocidades iniciales de los diferentes cuerpos
        
        var=1;                          #Variable utilizada en la función operación, para definir la actualización para los casos de k2 y k3

        h=tf-ti
           
        k1x=fx(vi)
        k1v=fv(ti,xi,vi,m,n)
         
        xiv,viv,k1x,k1v=operacion(k1x,k1v,xi,vi,n,var,h)                        #Valores que debo entrar al metodo ruge kutta 
        k2x=fx(viv)
        k2v=fv(ti+h/2.,xiv,viv,m,n)
        
         
        xiv,viv,k2x,k2v=operacion(k2x,k2v,xi,vi,n,var,h)
        k3x=fx(viv)
        k3v=fv(ti+h/2.,xiv,viv,m,n)
         
        var=0
        xiv,viv,k3x,k3v=operacion(k3x,k3v,xi,vi,n,var,h)
        k4x=fx(viv)
        k4v=fv(ti+h,xiv,viv,m,n)

        xiv,viv,k4x,k4v=operacion(k4x,k4v,xi,vi,n,var,h)
        
        
        xiv=operacion2(k1x,k2x,k3x,k4x,xi,n)
        xf=xiv

        viv=operacion2(k1v,k2v,k3v,k4v,vi,n)
        vf=viv
        
        return xf,vf
        
##########################################################################################################################################################
##########################################################################################################################################################
#SECCION 2:
#CUERPO DEL PROGRAMA


#SECCION 2-1:
#Nombramiento de variables   
    
nc=11                                           #Numero de cuerpos

ni=1.                                           #Numero de iteraciones en Ttot
Ttot=20000.
dt=Ttot/ni                                      #Cantidad de "tiempo" en que se actializará el sistema
ti=0.
es=30.                                          #Valor de escala para poder apreciar translacion de la luna
ca=pi/180.                                      #Conversion angulo a radianes

#Inicializacion de listas necesarias para el sistema, el numero de componentes de la lista la define nc

mc=[0]*nc               #Lista de las masas
xic=[0]*nc              #Lista de la posición inicial
xfc=[0]*nc              #Lista de la posición final
vic=[0]*nc              #Lista de la velocidad inicial
vfc=[0]*nc              #Lista de la velocidad final

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#SECCION 2-2:
#Condiciones iniciales

#Inicializacion masas

mc[0]=1.989e30                  #Masa Sol
mc[1]=3.302e23                  #Masa Mercurio
mc[2]=4.869e24                  #Masa Venus    
mc[3]=5.9736e24                 #Masa Tierra
mc[4]=6.4185e24                 #Masa Marte
mc[5]=7.349e22                  #Masa Luna
mc[6]=1.899e27                  #Masa Jupiter
mc[7]=5.688e26                  #Masa Saturno
mc[8]=8.686e25                  #Masa Urano
mc[9]=1.024e26                  #Masa Neptuno
mc[10]=1.25e22                  #Masa Pluton
#------------------------------------------------------------------------------------

#Inicialización posiciones iniciales y finales

#A cada componente de la lista xic y vic se inicializa con el vector posición inicial y velocidad inicial respectivamente para evitar problemas la posicion
#final se iguala a la posicion inicial antes de comenzar el programa          

#Cuerpo 1  Sol
vic[0]=vector(0.,0.,0.)
xic[0]=vector(0.,0.,0.)
xfc[0]=xic[0]

#Cuerpo 2   Mercurio

vp=59.25e3                 #velocidad perihelio
pp=46.00e9                 #posicion perihelio
a=29.12478*ca               #angulo del perihelio
i=7.00487*ca
vic[1]=vector(-vp*sin(a),vp*cos(a),0)
xic[1]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[1]=xic[1]

#cuerpo 3   Venus
vp=35.26e3                  #velocidad perihelio
pp=107.48e9                 #posicion perihelio
a=54.8522*ca
i=3.39*ca
vic[2]=vector(-vp*sin(a),vp*cos(a),0)
xic[2]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[2]=xic[2]

#cuerpo4    Tierra

vp=30.288e3                  #velocidad perihelio
pp=147.10e9                #posicion perihelio
a=282.94*ca
i=0.00*ca
vic[3]=vector(-vp*sin(a),vp*cos(a),0)
xic[3]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[3]=xic[3]

#cuerpo5    Marte
vp=26.614e3                  #velocidad perihelio
pp=206.67e9                  #posicion perihelio
a=286.46*ca
i=1.85*ca
vic[4]=vector(-vp*sin(a),vp*cos(a),0)
xic[4]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[4]=xic[4]

#cuerpo6 Luna

vp=30.288e3-1000           #velocidad perihelio
pp=147.10e9+3.84e8         #posicion perihelio
a=282.94*ca
i=0.00*ca
vic[5]=vector(-vp*sin(a),vp*cos(a),0)
xic[5]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[5]=xic[5]

#cuerpo7    Jupiter
vp=13.73e3                  #velocidad perihelio
pp=740.57e9                 #posicion perihelio
a=274.19*ca
i=1.3*ca
vic[6]=vector(-vp*sin(a),vp*cos(a),0)
xic[6]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[6]=xic[4]

#cuerpo 8   Saturno
vp=10.195e3                  #velocidad perihelio
pp=1353.6e9                  #posicion perihelio
a=339.39*ca
i=2.48*ca
vic[7]=vector(-vp*sin(a),vp*cos(a),0)
xic[7]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[7]=xic[7]

#cuerpo 9   Urano
vp=7.135e3                  #velocidad perihelio
pp=2748.9e9                 #posicion perihelio
a=96.73*ca
i=0.769*ca
vic[8]=vector(-vp*sin(a),vp*cos(a),0)
xic[8]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[8]=xic[8]

#cuerpo 10  Neptuno
vp=5.48e3                 #velocidad perihelio
pp=4452.9e9                 #posicion perihelio
a=272.449*ca
i=1.769*ca
vic[9]=vector(-vp*sin(a),vp*cos(a),0)
xic[9]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[9]=xic[9]

#cuerpo 11  Pluton
vp=6.2618e3                  #velocidad perihelio
pp=4438.6e9                #posicion perihelio
a=113.76*ca
i=17.15*ca
vic[10]=vector(-vp*sin(a),vp*cos(a),0)
xic[10]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[10]=xic[10]

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

#SECCION 2-3
#Definicion de los cuerpos que se visualizarán
#A los cuales se les adjudica posición inicial igual a los vectores de posicion iniciales definidos anteriormente


ball1 = sphere(pos=xic[0],radius=5e9, color=color.yellow)
ball1.trail = curve(color=ball1.color)

ball2 = sphere(pos=xic[1],radius=3e9, color=(0.5,0.5,0.5))
ball2.trail = curve(color=ball2.color)

ball3 = sphere(pos=xic[2],radius=3e9, color=color.orange)
ball3.trail = curve(color=ball3.color)

ball4 = sphere(pos=xic[3],radius=3e9, color=color.blue)
ball4.trail = curve(color=ball4.color)

ball5 = sphere(pos=xic[4],radius=3e9, color=(0.5,0.2,0.))
ball5.trail = curve(color=ball5.color)

ball6 = sphere(pos=(xic[5]-xic[3])*es+xic[3],radius=1e9, color=(0.3,0.3,0.3))        #La posicion inicial de este cuerpo es escaladaa para poder visualizar
ball6.trail = curve(color=ball6.color)

ball7 = sphere(pos=xic[6],radius=4e9, color=(0.058,0.728,0.488))
ball7.trail = curve(color=ball7.color)

ball8 = sphere(pos=xic[7],radius=4e9, color=(0.8,0.4,0.4))
ball8.trail = curve(color=ball8.color)

ball9 = sphere(pos=xic[7],radius=4e9, color=(1,1,1))
ball9.trail = curve(color=ball9.color)

ball10 = sphere(pos=xic[9],radius=4e9, color=(0.,0.3,1.))
ball10.trail = curve(color=ball10.color)

ball11 = sphere(pos=xic[10],radius=4e9, color=(0.8,0.6,0.0))
ball11.trail = curve(color=ball11.color)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

#SECCION 2-4
#BLUCLE

n=0

while True:                         #El bucle se ejecutará siempre

    #SECCION 2-4-1
    #Calculo valores finales de posicion y velocidad de cada cuerpo
    n=n+1                           #Contador usado como opcion auxiliar para parar el programa, (se debe poner como condición ej: While n<1000)
    rate(1e8)
    
    tf=ti+dt
    
    xfc,vfc=rk4(fv1,fx,xic,vic,tf,ti,mc,nc)         #LLAma la funciòn Ruge-Kutta para hallar los valores finales de posicion y velociadad 

    #Se actualizan los valores iniciales de velocidad, posición y tiempo
    vic=vfc                                         
    xic=xfc
    ti=tf           
   #--------------------------------------------------------------------------------------------------------------------------------------------------------
    #SECCION 2-4-2
    #Actualización final de los cuerpos para visualizar nueva posición
    ball1.pos = xfc[0]
    ball1.trail.append(pos=ball1.pos)

    ball2.pos = xfc[1]
    ball2.trail.append(pos=ball2.pos,retain=1000)

    ball3.pos = xfc[2]
    ball3.trail.append(pos=ball3.pos,retain=1000)

    ball4.pos = xfc[3]
    ball4.trail.append(pos=ball4.pos,retain=1000)

    ball5.pos = xfc[4]
    ball5.trail.append(pos=ball5.pos,retain=5000)

    pl=(xfc[5]-xfc[3])*es+xfc[3]
    ball6.pos = pl
    ball6.trail.append(pos=ball6.pos,retain=400) #,retain=10000

    ball7.pos = xfc[6]
    ball7.trail.append(pos=ball7.pos)

    ball8.pos = xfc[7]
    ball8.trail.append(pos=ball8.pos)

    ball9.pos = xfc[8]
    ball9.trail.append(pos=ball9.pos)

    ball10.pos = xfc[9]
    ball10.trail.append(pos=ball10.pos)

    ball11.pos = xfc[10]
    ball11.trail.append(pos=ball11.pos)

#Fin del bucle    
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

#Fin del programa
##############################################################################################################################################################
##############################################################################################################################################################
    
    
