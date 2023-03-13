from matrixMethodsV2 import *
from additionalFunc import *

import sys
from math import *

from var18 import *

#write out into file
sys.stdout = open("out.txt", "w")

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


#Функции U(t,x) и F(t,x)

def utx(t, x, var):
    return x+0.1*var*t*sin(pi*x)

def ftx(t, x, var, Xi):
    return 0.1*var*sin(pi*x)*(1.0+t*pi*pi*Xi)


#Векторы X и T

def X_vector(n):
    h = 1 / n
    x = [j*h for j in range(n+1)]
    return x, h

def T_vector(Xi, h, explicit=True):
    if explicit:
        tau = (h * h) / (4 * Xi)
    else:
        tau = h
    t = []
    n = 0
    while (abs(1.0 - n * tau) > (abs(1.0 - (n+1) * tau))):
        t.append(n * tau)
        n += 1
    t.append(n * tau)
    return t, tau


#Решение трехдиагональной матрицы методом прогонки

def TriDiag_Matrix(m, p, n):
    r = [0] * n
    a = [0] * n
    y = [0] * n
    b = [0] * n
    y[0] = m[0][0]
    a[0] = -m[0][1] / y[0]
    b[0] = p[0] / y[0]
    for i in range(1, n-1):
        y[i] = m[i][i] + m[i][i-1] * a[i-1]
        a[i] = -m[i][i+1] / y[i]
        b[i] = (p[i] - m[i][i-1] * b[i-1]) / y[i]
    y[n-1] = m[n-1][n-1] + m[n-1][n-2] * a[n-2]
    b[n-1] = (p[n-1] - m[n-1][n-2] * b[n-2]) / y[n-1]
    r[n-1] = b[n-1]
    for i in range(n-2, -1, -1):
        r[i] = a[i] * r[i+1] + b[i]
    return r


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

print("Метод стрельб")
a, history = shooting_method(dZ,0,1,1,2,0,3,0.1,0.00001,0.0001,Y)
print(table_for_shooting_method(history))
_, history = runge_kutta(dZ,0,1,[a,1],0.1,0.00001,0.00625,Y)
print(table_for_runge_kutta(history))
print(" xI=%4.2f\n" % 0.2)


#Явная конечно-разностная схема

def explicitScheme(x, t, var, Xi, tau, h):
    answer = ""
    maxDelta = -1
    u0 = [0] * len(x)
    u1 = [0] * len(x)
    u0[len(x)-1] = 1
    u1[len(x)-1] = 1
    for i in range(1,len(x)-1,1):
        u0[i] = x[i]   
    for i in range(len(t)-1):
        for j in range(1, len(x)-1, 1):
            u1[j] = \
            u0[j] + tau * Xi * ((u0[j+1]-2*u0[j]+u0[j-1])/(h*h))\
            + tau * ftx(t[i],x[j],var,Xi)
        for p in range(len(x)):
            c = u0[p]
            u0[p] = u1[p]
            u1[p] = c
        maxDeltaj = -1
        for k in range(len(x)):
            if (abs(utx(t[i+1],x[k],var)-u0[k]) > maxDeltaj):
                maxDeltaj = abs(utx(t[i+1],x[k],var)-u0[k])
        if (maxDeltaj > maxDelta):
            maxDelta = maxDeltaj
        answer += "%7.3f" % t[i+1]
        answer += "  %18.14E  " % (maxDeltaj)
        lenght = len(u0)
        lenght1 = lenght//3
        lenght2 = lenght1*2
        for i,ui in enumerate(u0):
            if i == lenght1:
                answer += "\n                               "
            if i == lenght2:
                answer += "\n                               "
            answer += "%9.5f" % ui
        answer += "\n"
    answer += "  Del_T=%18.17f\n" % maxDelta
    return answer
strAnswer = ""
for N in [8,16,32]:
    strAnswer += " N={}\n".format(N)
    x, h = X_vector(N)
    strAnswer += getHeaderTable(x)
    t_e, tau_e = T_vector(0.2, h)
    strAnswer += explicitScheme(x, t_e, var, 0.2, tau_e, h)
print("Метод конечных разностей (явная схема)")
print(strAnswer)


#Неявная конечно-разностная схема

def implicitScheme(x, t, var, Xi, tau, h):
    answer = ""
    maxDelta = -1
    u0 = [0] * len(x)
    u0[len(x)-1] = 1
    d = (tau * Xi) / (h*h)
    for i in range(1,len(x)-1,1):
        u0[i] = x[i]
    n = len(x)-2
    m = [[0] * n for _ in range(n)]
    p = [0] * n
    m[0][0] = (1+2*d)
    m[0][1] = -d
    m[n-1][n-2] = -d
    m[n-1][n-1] = (1+2*d)
    for i in range(1, n-1):
        m[i][i-1] = -d
        m[i][i] = (1+2*d)
        m[i][i+1] = -d    
    for j in range(len(t)-1):
        for i in range(0, n):
            p[i] = u0[i+1] + tau * ftx(t[j+1],x[i+1],var,Xi)
        p[n-1] += d
        u0 = [0] + TriDiag_Matrix(m, p, n) + [1]
        maxDeltaj = -1
        for k in range(len(x)):
            if (abs(utx(t[j+1],x[k],var)-u0[k]) > maxDeltaj):
                maxDeltaj = abs(utx(t[j+1],x[k],var)-u0[k])
        if (maxDeltaj > maxDelta):
            maxDelta = maxDeltaj
        answer += "%7.3f" % t[j+1]
        answer += "  %18.14E  " % (maxDeltaj)
        lenght = len(u0)
        lenght1 = lenght//3
        lenght2 = lenght1*2
        for i,ui in enumerate(u0):
            if i == lenght1:
                answer += "\n                               "
            if i == lenght2:
                answer += "\n                               "
            answer += "%9.5f" % ui
        answer += "\n"
    answer += "  Del_T=%18.17f\n" % maxDelta
    return answer
strAnswer = ""
for N in [8,16,32]:
    strAnswer += " N={}\n".format(N)
    x, h = X_vector(N)
    strAnswer += getHeaderTable(x)
    t_e, tau_e = T_vector(0.2, h)
    strAnswer += implicitScheme(x, t_e, var, 0.2, tau_e, h)
print("Метод конечных разностей (неявная схема)")
print(strAnswer)