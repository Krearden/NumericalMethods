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
                    if (element == 0):
                        print("{:10.14f}  ".format(element), end="")
                    else:
                        print("{:.14}  ".format(element), end="")
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

# #getSummOfMultiplyOfNonDiagonalElements
# def getSummOfMultiplyOfNonDiagonalElements(A):
#     summ = 0.0
#     n = len(A)
#     for i in range(n):
#         for j in range(i + 1, n):
#             summ += A[i][j] * A[i][j]
#     summ *= 2
#     return summ

#find max matrix element (module) which is not on the diagonal
def getMaxNonDiagonalElement(A):
    n = len(A)
    max_el = abs(A[0][1])
    max_el_i = 0
    max_el_j = 1
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
    max_non_diag = abs(A[0][1])
    A_copy = [[A[i][j] for j in range(n)] for i in range(n)]
    while (max_non_diag > epsilon):
        max_non_diag = abs(A[0][1])
        I = 0
        J = 1
        for i in range(n):
            for j in range(i):
                if (abs(A[i][j] > max_non_diag)):
                    max_non_diag = abs(A[i][j])
                    I = i
                    J = j

        theta = math.atan(2 * A[I][J] / (A[I][I] - A[J][J])) / 2
        c = math.cos(theta)
        s = math.sin(theta)

        for m in range(n):
            if (m != I and m != J):
                A_copy[m][I] = A_copy[I][m] = c * A[I][m] + s * A[J][m]
                A_copy[m][J] = A_copy[J][m] = -s * A[I][m] + c * A[J][m]

        A_copy[I][I] = c * c * A[I][I] + 2 * s * c * A[I][J] + s * s * A[J][J]
        A_copy[J][J] = s * s * A[I][I] - 2 * s * c * A[I][J] + c * c * A[J][J]
        A_copy[I][J] = A_copy[J][I] = (c * c - s * s) * A[I][J] + s * c * (A[J][J] - A[I][I])

        A = [[A_copy[i][j] for j in range(n)] for i in range(n)]

    max_diagonal = 0
    for i in range(n):
        for j in range(n):
            if (abs(A[i][j] > max_diagonal)):
                max_diagonal = abs(A[i][j])

    return math.sqrt(max_diagonal)

# def vectorSubstraction(vector_one, vector_two):
#     n = len(vector_one)
#     if (n == len(vector_two)):
#         result = [0 for i in range(n)]
#         for i in range(n):
#             result[i] = vector_one[i] - vector_two[i]
#         return result
#     else:
#         return

#MAIN
filepath = "files/"
filename = "test_1.txt"
A = readMatrixFromFile(filepath + filename) #матрица A
A, L, U, p, sign, rank = LU(A)
n = len(A)

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

third_norma_A = getEucledianNorm(A)
third_norma_A_inverse = getEucledianNorm(A_inverse)
tres = third_norma_A * third_norma_A_inverse
print("Евклидова норма: {:10.6f}".format(tres))

# print()
# print("A*x-b:")
# whatever2 = vectorSubstraction(multipMatrixByVector(A, x), b)
# for i in range(n):
#     print("{:.14} ".format(whatever2[i]), end='')

print()
print("Погрешность нахождения x: ")
given_x = [1, 2, 3, 4]
temp = [given_x[i] - x[i] for i in range(n)]
for i in range(n):
    print("{:.14} ".format(temp[i]), end='')