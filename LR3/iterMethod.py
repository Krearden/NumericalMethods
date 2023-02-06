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


def Iteratiobackupn(x, y):
    round = 0
    error = eps+1

    while(True):
        round+=1

        newX = 0.8-m.cos(y-1)
        newY = 2+m.cos(x)

        x = newX
        y = newY

        mJac = Jacobian(x, y)
        norm = getEucledianNorm(mJac)
        error = norm
        print(round, x, y, norm)

def lenght(x, y):
    return m.sqrt(x * x + y * y)

def Iterationbacku0p(x, y):
    round = 0

    mJac = Jacobian(x, y)
    norm = getEucledianNorm(mJac)

    while(True):
        round+=1

        newX = 0.8-m.cos(y-1)
        newY = 2+m.cos(x)

        error = norm / (1 - norm) * lenght(newX - x, newY - y)

        if (error < eps):
            break

        mJac = Jacobian(x, y)
        norm = getEucledianNorm(mJac)

        print(round, x, y, norm,)

        x = newX
        y = newY
    return x,y

def Iteration(x, y):
    round = 0

    mJac = Jacobian(x, y)
    norm = getEucledianNorm(mJac)
    error = 1 + eps
    print(error)

    while(error > eps):
        round+=1

        newX = 0.8-m.cos(y-1)
        newY = 2+m.cos(x)

        error = norm / (1 - norm) * lenght(newX - x, newY - y)

        mJac = Jacobian(x, y)
        norm = getEucledianNorm(mJac)

        print(round, x, y, norm,)

        x = newX
        y = newY
    return x,y
if __name__ == "__main__":
    x0 = 0.900
    y0 = 2.600

    Iteration(x0, y0)
