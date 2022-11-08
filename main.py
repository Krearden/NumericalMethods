#coding=utf-8
from matrixMethods import *
import math

# Лабораторная работа № 2 на тему «Прямые методы решения систем линейных алгебраических уравнений»
# Выполнили Запорожченко Кирилл и Педаев Михаил (ФЗ-11)
# 2022 г.


#FUNCTIONS

#read matrix from file by given path/filename.txt
def readMatrixFromFile(filename):
    row = []
    matrix = []
    with open(filename) as file:
        for line in file:
            line = line.strip().split()
            for element in line:
                row.append(float(element))
            matrix.append(row)
            row = []
    return matrix

#print matrix on the screen
def printMatrix(matrix, accuracy = False):
    if matrix:
        for row in matrix:
            for element in row:
                if accuracy:
                    print("{:10.20f} ".format(element), end="")
                else:
                    print("{:10.6f} ".format(element), end="")
            print()
    else:
        print("Empty matrix given")

#единичная матрица
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

# def getRankMatrixFromU(U):
#     n = len(U)
#     temp = 0
#     rank = 0
#     for i in range(n):
#         for j in range(n):
#             temp += U[i][j]
#         if (temp != 0):
#             rank += 1
#         temp = 0
#
#     return rank

def computeDeterminantFromL(L, sign):
    n = len(L)
    det = 1
    for i in range(n):
        det *= L[i][i]
    det *= sign

    return det

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

        print("i = {}; max_el_row_index = {};".format(i, max_el_row_index))
        print("U: ")
        printMatrix(U)
        print("L: ")
        printMatrix(L)
        print("\n\n")

    return A, L, U, p, sign, rank

def checkIfCorrectLU(A, L, U, p):
    P = getPmatrix(p)
    LU = multiplyMatrix(L, U)
    PA = multiplyMatrix(P, A)

    return matrixSubtraction(LU, PA)


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
    n = len(L)
    x = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        summ = 0
        for j in range(n - 1, i - 1, -1):
            summ += U[i][j] * x[j]
        x[i] = (y[i] - summ) / U[i][i]
    return x

#get inverse matirx
def getInverseMatrix(A, p):
    n = len(A)
    A_inverse = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        e_i = [1 if j == i else 0 for j in range(n)]
        e_i = multipMatrixByVector(getPmatrix(p), e_i)
        y = forwardSubstitution(L, e_i)
        x = backwardSubstitution(U, y)
        for k in range(n):
            A_inverse[k][i] = x[k]
    return A_inverse

#getFirstNorma of matrix
def getFirstNorma(A):
    n = len(A)
    summ = 0
    max_summ = 0
    for i in range(n):
        for j in range(n):
            summ += abs(A[i][j])
        if summ > max_summ:
            max_summ = summ
        summ = 0
    return max_summ

#getSecondNorma of matrix
def getSecondNorma(A):
    n = len(A)
    summ = 0
    max_summ = 0
    for i in range(n):
        for j in range(n):
            summ += abs(A[j][i])
        if summ > max_summ:
            max_summ = summ
        summ = 0
    return max_summ

#getSummOfMultiplyOfNonDiagonalElements
def getSummOfMultiplyOfNonDiagonalElements(A):
    summ = 0.0
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            summ += A[i][j] * A[i][j]
    summ *= 2
    return summ

#find max matrix element (module) which is not on the diagonal
def getMaxNonDiagonalElement(A):
    n = len(A)
    max_el = 0.0
    max_el_i = max_el_j = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (abs(A[i][j] > max_el)):
                max_el = abs(A[i][j])
                max_el_i = i
                max_el_j = j
    return max_el_i, max_el_j

#get transposed matrix
def getTransposedMatrix(A):
    n = len(A)
    A_transposed = [[A[j][i] for j in range(n)] for i in range(n)]
    return A_transposed

#get eucledian norm
def getEucledianNorm(A):
    n = len(A)
    epsilon = 1e-06
    A = multiplyMatrix(getTransposedMatrix(A), A)
    max_a = abs(A[0][1])
    A_new = [[A[i][j] for j in range(n)] for i in range(n)]
    while (max_a > epsilon):
        ii = 0
        jj = 1
        for i in range(n):
            for j in range(i):
                if (abs(A[i][j] > max_a)):
                    max_a = abs(A[i][j])
                    ii = i
                    jj = j

        arctan_half = 0.5 * math.atan(2 * A[ii][jj] / (A[ii][ii] - A[jj][jj]))
        s = math.sin(arctan_half)
        c = math.cos(arctan_half)

        for m in range(n):
            if (m != ii and m != jj):
                A_new[m][ii] = A_new[ii][m] = c * A[ii][m] + s * A[jj][m]
                A_new[m][jj] = A_new[jj][m] = -s * A[ii][m] + c * A[jj][m]

        A_new[ii][ii] = c * c * A[ii][ii] + 2 * s * c * A[ii][jj] + s * s * A[jj][jj]
        A_new[jj][jj] = s * s * A[ii][ii] - 2 * s * c * A[ii][jj] + c * c * A[jj][jj]
        A_new[ii][jj] = A_new[jj][ii] = (c * c - s * s) * A[ii][jj] + s * c * (A[jj][jj] - A[ii][ii])

        A = [[A_new[i][j] for j in range(n)] for i in range(n)]
    return A

# #get eucledian norm
# def getEucledianNorm(A):
#     n = len(A)
#     epsilon = 1e-06
#     max_diagonal_element = 0.0
#     AAT = multiplyMatrix(getTransposedMatrix(A), A)
#     while (math.sqrt(getSummOfMultiplyOfNonDiagonalElements(AAT)) > epsilon):
#         print(getSummOfMultiplyOfNonDiagonalElements(AAT), math.sqrt(max_diagonal_element))
#         max_el_i, max_el_j = getMaxNonDiagonalElement(AAT)
#         if (AAT[max_el_i][max_el_i] == AAT[max_el_j][max_el_j]):
#             theta = math.pi / 4
#         else:
#             theta = math.atan(2 * A[max_el_i][max_el_j] / (A[max_el_i][max_el_i] - A[max_el_j][max_el_j])) / 2
#         C = math.cos(theta)
#         S = math.sin(theta)
#         B = [[AAT[i][j] for j in range(n)] for i in range(n)]
#
#         for k in range(n):
#             B[k][max_el_i] = C * AAT[k][max_el_i] + S * AAT[k][max_el_j]
#             B[k][max_el_j] = -S * AAT[k][max_el_i] + C * AAT[k][max_el_j]
#         for k in range(n):
#             AAT[max_el_i][k] = AAT[k][max_el_i] = C * B[max_el_i][k] + S * B[max_el_j][k]
#             AAT[max_el_j][k] = AAT[k][max_el_j] = -S * B[max_el_i][k] + C * B[max_el_j][k]
#
#         AAT[max_el_i][max_el_j] = AAT[max_el_j][max_el_i] = 0
#
#     for i in range(n):
#         if (abs(AAT[i][i]) > max_diagonal_element):
#             max_diagonal_element = abs(AAT[i][i])
#
#     return math.sqrt(max_diagonal_element)

#MAIN
filepath = "files/"
filename = "test_1.txt"
A = readMatrixFromFile(filepath + filename) #матрица A
A, L, U, p, sign, rank = LU(A)

print("Result of LU decomposition:" )
print("A: ")
printMatrix(A)
print("U: ")
printMatrix(U)
print("L: ")
printMatrix(L)
print("P: ")
printMatrix(getPmatrix(p))
print()
print("Vector p:")
print(p)

print()
print("Check if LU decomposition is correct (L*U-P*A): ")
should_be_zeros = checkIfCorrectLU(A, L, U, p)
printMatrix(should_be_zeros, True)

print()
print("Rank = {}".format(rank))

determinant = computeDeterminantFromL(L, sign)
print()
print("Determinant = {}".format(determinant))

#find x
print()
print("x (found): ")
b = multipMatrixByVector(A, [1, 2, 3, 4])
b = multipMatrixByVector(getPmatrix(p), b)
y = forwardSubstitution(L, b)
x = backwardSubstitution(U, y)
print(x)

#inverse matirx
print()
print("A^(-1):")
A_inverse = getInverseMatrix(A, p)
printMatrix(A_inverse)

#A*A^(-1)-E
print()
print("A*A^(-1)-E")
n = len(A)
E = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
whatever = matrixSubtraction(multiplyMatrix(A, A_inverse), E)
printMatrix(whatever, True)

#normas
print()
print("Число обусловленности в трех нормах:")
first_norma_A = getFirstNorma(A)
first_norma_A_inverse = getFirstNorma(A_inverse)
uno = first_norma_A * first_norma_A_inverse
print("Норма 1: {:10.6f}".format(uno))

second_norma_A = getSecondNorma(A)
second_norma_A_inverse = getSecondNorma(A_inverse)
dos = second_norma_A * second_norma_A_inverse
print("Норма 2: {:10.6f}".format(dos))


printMatrix(getEucledianNorm(A))
# third_norma_A = getEucledianNorm(A)
# third_norma_A_inverse = getEucledianNorm(A_inverse)
# tres = third_norma_A * third_norma_A_inverse
# print("Евклидова норма: {:10.6f}".format(tres))
