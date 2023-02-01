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


def Iteration(x, y):
    round = 0
    error = eps+1

    while(error>eps):
        round+=1

        newX = 0.8-m.cos(y-1)
        newY = 2+m.cos(x)

        x = newX
        y = newY

        mJac = Jacobian(x, y)
        norm = getEucledianNorm(mJac)
        error = norm
        print(round, x, y, norm)

if __name__ == "__main__":
    x0 = 0.900
    y0 = 2.600

    Iteration(x0, y0)
