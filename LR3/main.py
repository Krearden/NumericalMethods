#coding=utf-8
import math
from matrixMethods import *
from LU import *

# Лабораторная работа № 3 на тему «Решение нелинейных уравнений»
# Выполнили Запорожченко Кирилл и Педаев Михаил (ФЗ-11)
# 2023 г.


#FUNCTIONS
#print matrix on the screen
def printMatrix(matrix, accuracy = False):
    if matrix:
        for row in matrix:
            for element in row:
                if accuracy:
                    if (element == 0):
                        print("{:10.14f}  ".format(element), end="")
                    else:
                        print("{:.14}  ".format(element), end="")
                else:
                    print("{:10.6f} ".format(element), end="")
            print()
    else:
        print("Empty matrix given")

#get eucledian norm
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


#get eucledian norm of vector
def getEucledianVectorNorm(vector):
    pass


#Jacobi marix
def getJacobian():
    pass

def getFx(X, variant):
    fX = [0 for i in range(len(X))]
    if (variant == 0):
        X[0] = math.sin(X[0]) + 2 * X[1] - 1.6
        X[1] = math.cos(X[1] - 1) + X[0] - 1
    else:
        print("Такого еще нет варианта")

    return X


#Derivative matrix - матрица производных для систем в неявном виде
def getDerivativeMatrix(variant, x, y):
    N = 2
    dF = [[0 for i in range(N)] for i in range(N)]
    #test variant
    if (variant == 0):
        dF[0][0] = math.cos(x) #f1 over x
        dF[0][1] = 1 #f2 over x
        dF[1][0] = 2 #f1 over y
        dF[1][1] = -math.sin(y - 1) #f2 over y
    else:
        print("Пока здесь такого варианта нету")

    return dF


def length(X, variant):
    Fx = getFx(X, variant)
    sqrt = math.sqrt(Fx[0] * Fx[0] + Fx[1] * Fx[1])
    return sqrt


def getDeterminant(matrix):
    if (len(matrix) != 2):
        print("Wrong matrix size")
    else:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

#Newton's method
def methodNewton(X, variant):
    epsilon = 1e-04
    while (length(X, variant) > epsilon):
        #матрица производных
        dF = getDerivativeMatrix(variant, X[0], X[1])
        #определители
        determinant = getDeterminant(dF)
        #fx
        Fx = getFx(X, variant)
        #for x
        dF_x = [[dF[i] for i in range(len(dF))] for i in range(len(dF))]
        dF_x[0][0] = Fx[0]
        dF_x[1][0] = Fx[1]
        determinant_x = getDeterminant(dF_x)
        #for y
        dF_y = [[dF[i] for i in range(len(dF))] for i in range(len(dF))]
        dF_y[0][1] = Fx[0]
        dF_y[1][1] = Fx[1]
        determinant_y = getDeterminant(dF_y)
        #find new x and y
        newX = [0 for i in range(len(dF))]
        newX[0] = X[0] - determinant_x / determinant #find new X
        newX[1] = X[1] - determinant_y / determinant #find new Y
        #copy
        X = [newX[i] for i in range(len(dF))]
    print(X)
    return X

#MAIN
X = [0, 0.8] # [0][0] is X and [0][1] is Y (примерные значения)
methodNewton(X, 0)