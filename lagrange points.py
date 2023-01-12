#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import math


# In[6]:


u=3.0039e-7
L1=0.995363
L2=1.004637
L3=-1.00001


# In[9]:


def Fx(x,y,vx,vy,t):
    x1,x2=-u,1-u
    y1,y2=0,0
    p1=np.sqrt((x-x1)**2+(y-y1)**2)
    p2=np.sqrt((x-x2)**2+(y-y2)**2)
    return 2*vy+x-(1-u)*(x-x1)/(p1^3)-u*(x-x2)/(p2^3)
def Fy(x,y,vx,vy,t):
    x1,x2=-u,1-u
    y1,y2=0,0
    p1=np.sqrt((x-x1)**2+(y-y1)**2)
    p2=np.sqrt((x-x2)**2+(y-y2)**2)
    return -2*vx+y-(1-u)*(y-y1)/(p1^3)-u*(y-y2)/(p2^3)
def EulerRichardson(x0,vx0,y0,vy0,start,end,steps):
    vx_out = []
    x_out = []
    vy_out = []
    y_out = []
    t = np.linspace(start,end,steps)
    h = t[1]-t[0]
    x_n = x0
    vx_n = vx0
    y_n = y0
    vy_n = vy0
    x_out.append(x0)
    vx_out.append(vx0)
    y_out.append(y0)
    vy_out.append(vy0)
    for i in range(1,len(t)):
        ax_n = Fx(x_n,y_n,vx_n,vy_n,t[i])/m
        ay_n = Fy(x_n,y_n,vx_n,vy_n,t[i])/m
        vx_mid = vx_n + 0.5*ax_n*h
        x_mid = x_n + 0.5*vx_n*h
        vy_mid = vy_n + 0.5*ay_n*h
        y_mid = y_n + 0.5*vy_n*h
        ax_mid = Fx(x_mid,y_mid,vx_mid,vy_mid,t[i]+0.5*h)/m
        ay_mid = Fy(x_mid,y_mid,vx_mid,vy_mid,t[i]+0.5*h)/m
        vx_n_plus_1 = vx_n + ax_mid*h
        x_n_plus_1 = x_n + vx_mid*h
        vy_n_plus_1 = vy_n + ay_mid*h
        y_n_plus_1 = y_n + vy_mid*h
        vx_out.append(vx_n_plus_1)
        x_out.append(x_n_plus_1)
        vy_out.append(vy_n_plus_1)
        y_out.append(y_n_plus_1)
        x_n = x_n_plus_1
        vx_n = vx_n_plus_1
        y_n = y_n_plus_1
        vy_n = vy_n_plus_1
    return x_out,vx_out,y_out,vy_out,t


# In[10]:


x_out,vx_out,y_out,vy_out,t = EulerRichardson(0.5,0,0.87,0,0,1000,100000)


# In[ ]:




