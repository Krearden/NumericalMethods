#coding=utf-8
import math
from matrixMethods import *

# Лабораторная работа № 3 на тему «Решение нелинейных уравнений»
# Выполнили Запорожченко Кирилл и Педаев Михаил (ФЗ-11)
# 2022 г.


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

#MAIN
