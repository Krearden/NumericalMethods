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
                print("{:10.5f} ".format(element), end="")
            print()
    else:
        print("Empty matrix given")


def getLUfromA(A):
    pass

def getUfromA(A, n, p):
    for k in range(n - 1):
        # находим индекс масксимального элемента в k-том столбце
        max_el = 0
        imax = -1
        for i in range(n):
            if(abs(A[i][k]) > abs(max_el)):
                max_el = A[i][k]
                imax = i
        #меняем строки местами
        if (imax != k):
            A = swapRows(A, k, imax)
            p = swapRows(p, k, imax)
        #проверка на вырожденность
        if (A[k][k] == 0):
            print("Вырожденная матрица")
            return
        #делаем волшебные действия (не работает) (нужно получить нули под элементом k, k)
        for i in range(k + 1):
            pass
    return A, p
def getLfromA(U):
    pass

def swapRows(A, i, j):
    temp = A[i]
    A[i] = A[j]
    A[j] = temp
    return A

#MAIN
filepath = "files/"
filename = "var_18_b.txt"
A = readMatrixFromFile(filepath + filename) #матрица A
n = len(A) #количество элементов матрицы
p = [0, 1, 2, 3] #вектор перестановок
printMatrix(A)
print()
printMatrix(getUfromA(A, n, p)[0])
print(getUfromA(A, n, p)[1])

