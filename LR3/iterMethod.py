import math
from matrixMethods import *


# Fi1(x,y)=0.8-cos(y-1)
# Fi2(x,y)=2+cos(x)



def getEucledianNorm(A):
    n = len(A)
    epsilon = 1e-06
    A = multiplyMatrix(getTransposedMatrix(A), A)
    max_non_diag = abs(A[0][1])
    A_copy = [[A[i][j] for j in range(n)] for i in range(n)]
    while (max_non_diag > epsilon):
        max_non_diag = abs(A[0][1])
        I = 0
        J = 1
        for i in range(n):
            for j in range(i + 1, n):
                if (abs(A[i][j]) > max_non_diag):
                    max_non_diag = abs(A[i][j])
                    I = i
                    J = j
        theta = math.atan(2 * A[I][J] / (A[I][I] - A[J][J])) / 2
        c = math.cos(theta)
        s = math.sin(theta)
        for m in range(n):
            if (m != I and m != J):
                A_copy[m][I] = c * A[I][m] + s * A[J][m]
                A_copy[m][J] = -s * A[I][m] + c * A[J][m]
                A_copy[I][m] = c * A[I][m] + s * A[J][m]
                A_copy[J][m] = -s * A[I][m] + c * A[J][m]
        A_copy[I][I] = c * c * A[I][I] + 2 * s * c * A[I][J] + s * s * A[J][J]
        A_copy[J][J] = s * s * A[I][I] - 2 * s * c * A[I][J] + c * c * A[J][J]
        A_copy[I][J] = (c * c - s * s) * A[I][J] + s * c * (A[J][J] - A[I][I])
        A_copy[J][I] = (c * c - s * s) * A[I][J] + s * c * (A[J][J] - A[I][I])
        A = [[A_copy[i][j] for j in range(n)] for i in range(n)]
    max_diagonal = 0
    for i in range(n):
        if (abs(A[i][i] > max_diagonal)):
            max_diagonal = abs(A[i][i])

    return math.sqrt(max_diagonal)

def Jacobian(x, y):
    jak = [
    [0.0, math.sin(y-1)],
    [-math.sin(x),  0.0 ]
    ]
    return jak

def devMatrix(x, y):
    jak = [
    [1.0, -1*math.sin(y-1)],
    [math.sin(x),  1.0 ]
    ]
    return jak

def F1(x, y):
    return x-0.8+math.cos(y-1)

def F2(x, y):
    return y-2-math.cos(x)


def lenght(x, y):
    return math.sqrt(x * x + y * y)


def Iteration(x, y):
    round = 0
    eps = 1e-4

    mJac = Jacobian(x, y)
    norm = getEucledianNorm(mJac)
    error = 1 + eps


    while(lenght(F1(x, y), F2(x, y)) > eps):
        round+=1

        newX = 0.8-math.cos(y-1)
        newY = 2+math.cos(x)

        mJac = Jacobian(x, y)
        norm = getEucledianNorm(mJac)

        error = norm / (1 - norm) * lenght(newX - x, newY - y)
        
        x = newX
        y = newY

        #print(round, x, y, norm, eps-lenght(F1(x, y), F2(x, y)))
        print("| {:3} |  {:3.8f} |  {:3.8f} |    {:1.8f} |   {: 3.8f} |   {: 3.8f} |".format(round, x, y, norm, F1(x, y), F2(x, y)))
    return x,y
if __name__ == "__main__":
    x0 = 0.900
    y0 = 2.600

    Iteration(x0, y0)
