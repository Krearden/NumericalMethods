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

#get LU matrix from A matrix
def getLUfromA(A_matrix):
    n = len(A_matrix)
    for j in range(n):
        for i in range(n):
            number = min(i, j) - 1
            summ_of_multiply = getSummOfMultiply_ij(A_matrix, i, j, number)
            A_matrix[i][j]  = A_matrix[i][j] - summ_of_multiply
            if (i > j):
                A_matrix[i][j] = A_matrix[i][j] / A_matrix[j][j]
    return A_matrix



# сумма произведений элементов i-той j-того столбца матрицы
def getSummOfMultiply_ij(matrix, i, j, n):
    summ = 0.0
    for k in range(n):
        summ += matrix[i][k] * matrix[k][j]
        print("{} * {} + \n".format(matrix[i][k], matrix[k][j]))
    return summ


#MAIN
filepath = "files/"
filename = "var_18_b.txt"
A_matrix = readMatrixFromFile(filepath + filename)
LU_matrix = getLUfromA(A_matrix)
print("A:")
printMatrix(A_matrix)
print("LU:")
printMatrix(LU_matrix)


