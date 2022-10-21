#coding=utf-8


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
def printMatrix(matrix):
    if matrix:
        for row in matrix:
            for element in row:
                print("{:10.6f} ".format(element), end="")
            print()
    else:
        print("Empty matrix given")

#multiply Matrix By Vector (to find b)
def multipMatrix(m, v):
    rv = []

    for row in m:
        rowResult = 0
        for elemPos in range(len(row)):
            rowResult += row[elemPos] * v[elemPos]
        rv.append(rowResult)

    return rv

#get P matrix from p vector
def getPmatrix(p):
    n = len(p)
    P = list([[0  for j in range(n)] for i in range(n)])
    for i in range(n):
        P[i][p[i]] = 1
    return P

#LU
def LU(A):
    n = len(A)
    p = list([i for i in range(n)]) #вектор перестановок
    vector_p = [0, 1, 2, 3]
    L = list([[0  for j in range(n)] for i in range(n)])
    U = list([[A[i][j]  for j in range(n)] for i in range(n)])
    sign = 1

    for i in range(n):
        #находим главный элемент в столбце i и меняем строки местами
        max_el = 0
        max_el_row_index = -1
        for k in range(i, n):
            if(abs(U[k][i]) > max_el):
                max_el = abs(U[k][i])
                max_el_row_index = k
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
            L[i][j] = A[vector_p[i]][j] - summ



        print("i = {}".format(i))
        printMatrix(U)
        print()

    print("Result of LU decomposition:" )
    print("A: ")
    printMatrix(A)
    print("U: ")
    printMatrix(U)
    print("L: ")
    printMatrix(L)
    print("P: ")
    printMatrix(getPmatrix(p))
    print(p)
    return A

#еденичная матрица
def swapRows(matrix, i, j):
    temp = matrix[i]
    matrix[i] = matrix[j]
    matrix[j] = temp
    return matrix

#MAIN
filepath = "files/"
filename = "test_1.txt"
A = readMatrixFromFile(filepath + filename) #матрица A
LU(A)
