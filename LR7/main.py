from matrixMethodsV2 import *
from additionalFunc import *

from math import *

from var24 import *


#Функции K в методе Рунге-Кутта

def k1(F,x,y):
    res = []
    res.append(F(x,y))
    for i in range(len(y)-1):
        res.append(y[i])
    return res


def k2(F,x,y,h,k1_res):
    res = []
    res.append(F(x + h/2,Addition(y,Multiplicate(h/2,k1_res),True)))
    for i in range(len(y)-1):
        res.append(y[i]+k1_res[i+1]*h/2)
    return res

def k3(F,x,y,h,k2_res):
    res = []
    res.append(F(x + h/2,Addition(y,Multiplicate(h/2,k2_res),True)))
    for i in range(len(y)-1):
        res.append(y[i] + k2_res[i+1]*h/2)
    return res

def k4(F,x,y,h,k3_res):
    res = []
    res.append(F(x + h,Addition(y,Multiplicate(h,k3_res),True)))
    for i in range(len(y)-1):
        res.append(y[i] + k3_res[i+1]*h)
    return res


#Метод Рунге-Кутта

def next_step(F,x,y,h):
    k1_res = k1(F,x,y)
    k2_res = k2(F,x,y,h,k1_res)
    k3_res = k3(F,x,y,h,k2_res)
    k4_res = k4(F,x,y,h,k3_res)
    k2_res = Multiplicate(2,k2_res,True)
    k3_res = Multiplicate(2,k3_res,True)
    k_res = Addition(Addition(Addition(k1_res,k2_res),k3_res),k4_res,True)
    k_res = Multiplicate(h/6,k_res,True)
    res = Addition(y,k_res,True)
    return res

def check_step(F,x,y,h,eps):
    y_j = next_step(F,x,y,h)
    the_next = y_j
    y_j2 = next_step(F,x,y,h/2)
    y_j = Multiplicate(-1,y_j,True)
    delta = Addition(y_j,y_j2,True)
    norma = Vector_norm(delta)
    return norma < eps, the_next

def runge_kutta(F,a,b,y0,h0,eps,interval,Y_real):
    history = []
    x = a
    y = y0
    h = h0
    s = 0
    history.append((x,Y_real(x),y[1],y[0],abs(Y_real(x) - y[1])))
    last = x
    while x < b:
        if s%1000 == 0:
            s = 0
            flag, the_next = check_step(F,x,y,h,eps)
            while flag:
                h *= 2
                flag, the_next = check_step(F,x,y,h,eps)
            h /= 2

        flag, the_next = check_step(F,x,y,h,eps)
        if flag:
            prev_y = y
            prev_x = x
            y = the_next
            x += h
            s += 1
            # print(x)
        else:
            h /= 2
        if abs(x-last) > interval:
            temp = abs(abs(x-last) - interval)
            temp_h = abs(x - temp - prev_x)
            temp_y = next_step(F,prev_x,prev_y,temp_h)
            temp_x = prev_x + temp_h
            history.append((temp_x,Y_real(temp_x),temp_y[1],temp_y[0],abs(Y_real(temp_x) - temp_y[1])))
            last = temp_x
    history.append((x,Y_real(x),y[1],y[0],abs(Y_real(x) - y[1])))
    return y[len(y)-1], history


#Метод стрельб

def shooting_method(F,x0,x1,y0,y1,a1,a2,h0,eps_rk,eps_sm,Y_real):
    history = []
    B1, _ = runge_kutta(F,x0,x1,[a1,y0],h0,eps_rk,0.1,Y_real)
    B2, _ = runge_kutta(F,x0,x1,[a2,y0],h0,eps_rk,0.1,Y_real)
    if (B1 - y1)*(B2 - y1) > 0:
        raise ValueError("B1 и B2 на одной стороне от целевого значения")
    if B1 > y1:
        a1, a2 = a2, a1
    history.append((0,a1,B1,abs(B1 - y1)))
    history.append((0,a2,B2,abs(B2 - y1)))
    print(B1)
    print(B2)
    B = B1
    i = 1
    while abs(B - y1) > eps_sm:
        a = (a1 + a2)/2
        B, _ = runge_kutta(F,x0,x1,[a,y0],h0,eps_rk,0.1,Y_real)
        if B < y1:
            a1 = a
        else:
            a2 = a
        print(B)
        history.append((i,a,B,abs(y1 - B)))
        i += 1 
    return a, history
a, history = shooting_method(dZ,0,1,1,2,0,3,0.1,0.00001,0.0001,Y)
print(table_for_shooting_method(history))
_, history = runge_kutta(dZ,0,1,[a,1],0.1,0.00001,0.00625,Y)
print(table_for_runge_kutta(history))
print(" xI=%4.2f\n" % 0.2)