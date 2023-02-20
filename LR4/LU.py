#coding=utf-8

#Для нахождения обратной матрицы

from matrixMethods import *

def swapRows(matrix, i, j):
    temp = matrix[i]
    matrix[i] = matrix[j]
    matrix[j] = temp
    return matrix

#get P matrix from p vector
def getPmatrix(p):
    n = len(p)
    P = list([[0  for j in range(n)] for i in range(n)])
    for i in range(n):
        P[i][p[i]] = 1
    return P


#forward substitution
def forwardSubstitution(L, b):
    n = len(L)
    y = [0 for i in range(n)]
    for i in range(n):
        summ = 0
        for j in range(i):
            summ += L[i][j] * y[j]
        y[i] = (b[i] - summ) / L[i][i]
    return y


#backward substitution
def backwardSubstitution(U, y):
    n = len(U)
    x = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        summ = 0
        for j in range(n - 1, i - 1, -1):
            summ += U[i][j] * x[j]
        x[i] = (y[i] - summ) / U[i][i]
    return x

#get inverse matirx
def getInverseMatrix(A):
    n = len(A)
    A, L, U, p, sign, rank = LU(A)
    A_inverse = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        e_i = [1 if j == i else 0 for j in range(n)]
        e_i = multipMatrixByVector(getPmatrix(p), e_i)
        y = forwardSubstitution(L, e_i)
        x = backwardSubstitution(U, y)
        for k in range(n):
            A_inverse[k][i] = x[k]
    return A_inverse

#LU
def LU(A):
    n = len(A)
    p = list([i for i in range(n)]) #вектор перестановок
    L = list([[0  for j in range(n)] for i in range(n)])
    U = list([[A[i][j]  for j in range(n)] for i in range(n)])
    sign = 1
    epsilon = 1e-06
    rank = 0
    for i in range(n):
        #находим главный элемент в столбце i и меняем строки местами
        max_el = 0
        max_el_row_index = -1
        for k in range(i, n):
            if(abs(U[k][i]) > max_el):
                max_el = abs(U[k][i])
                max_el_row_index = k
        if (max_el > epsilon):
            rank += 1
        if (max_el != 0 and max_el_row_index != i):
            p = swapRows(p, i, max_el_row_index)
            U = swapRows(U, i, max_el_row_index)
            sign *= -1
        #метод гаусса с перестановкой строк (формируем матрицу U)
        main_element = U[i][i]
        L[i][i] = main_element
        for row_element in range(n):
            U[i][row_element] /= main_element
        for k in range(i + 1, n):
            element_to_multiply = U[k][i] #- элемент на который умножать (первый элемент каждой строки ниже i-той)
            for j in range(i, n):
                U[k][j] = U[k][j] - U[i][j] * element_to_multiply
        #формируем матрицу L
        for j in range(i):
            summ = 0
            for k in range(j):
                summ += L[i][k] * U[k][j]
            L[i][j] = A[p[i]][j] - summ

    return A, L, U, p