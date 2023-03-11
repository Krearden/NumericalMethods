from math import log, sin
def A(x):
    return 40*(-x+0.5)
def B(x):
    return -x**2+2
def C(x):
    return x+2
def Y(x):
    return 1+x+10*log(24+1)*x**3*(1-x)**3
def dY(x):
    return -30*x**3*(1-x)**2*log(25)+30*x**2*(1-x)**3*log(25) + 1
def ddY(x):
    return -30*x**3*(2*x-2)*log(25)-180*x**2*(1-x)**2*log(25)+60*x*(1-x)**3*log(25)
def F(x):
    return ddY(x) + A(x)*dY(x) - B(x)*Y(x) + C(x)*sin(Y(x))
def dZ(x,y):
    return F(x) - C(x)*sin(y[1]) + B(x)*y[1] - A(x)*y[0]
var = 24