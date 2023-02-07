import math as m
from matrixMethods import *
from main import getEucledianNorm

# Fi1(x,y)=0.8-cos(y-1)
# Fi2(x,y)=2+cos(x)

eps = 1e-4

def Jacobian(x, y):
    jak = [
    [0.0, m.sin(y-1)],
    [-m.sin(x),  0.0 ]
    ]
    return jak

def devMatrix(x, y):
    jak = [
    [1.0, -1*m.sin(y-1)],
    [m.sin(x),  1.0 ]
    ]
    return jak

def F1(x, y):
    return x-0.8+m.cos(y-1)

def F2(x, y):
    return y-2-m.cos(x)


def lenght(x, y):
    return m.sqrt(x * x + y * y)


def Iteration(x, y):
    round = 0

    mJac = Jacobian(x, y)
    norm = getEucledianNorm(mJac)
    error = 1 + eps
    print(error)

    while(lenght(F1(x, y), F2(x, y)) > eps):
        round+=1

        newX = 0.8-m.cos(y-1)
        newY = 2+m.cos(x)

        mJac = Jacobian(x, y)
        norm = getEucledianNorm(mJac)

        error = norm / (1 - norm) * lenght(newX - x, newY - y)
        
        x = newX
        y = newY

        print(round, x, y, norm, eps-lenght(F1(x, y), F2(x, y)))
    return x,y
if __name__ == "__main__":
    x0 = 0.900
    y0 = 2.600

    Iteration(x0, y0)
