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

#LU
def LU(A):
    n = len(A)
    P = list([[1 if i == j else 0 for j in range(n)] for i in range(n)]) #матрица перестановок

    for i in range(n):

        #находим главный элемент в столбце i и меняем строки местами
        max_el = 0
        max_el_row_index = -1
        for k in range(i, n):
            if(abs(A[k][i]) > max_el):
                max_el = abs(A[k][i])
                max_el_row_index = k
        if (max_el != 0 and max_el_row_index != i):
            P = swapRows(P, i, max_el_row_index)
            A = swapRows(A, i, max_el_row_index)

        main_element = A[i][i]
        for row_element in range(n):
            A[i][row_element] /= main_element
        for k in range(i + 1, n):
            element_to_multiply = A[k][i] #- элемент на который умножать (первый элемент каждой строки ниже i-той)
            for j in range(i, n):
                A[k][j] = A[k][j] - A[i][j] * element_to_multiply
        print("i = {}".format(i))
        printMatrix(A)
        print()
    print("U matirx made from A: ")
    printMatrix(A)
    print("P: ")
    printMatrix(P)
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
print("A:")
printMatrix(A)
print()
LU(A)

