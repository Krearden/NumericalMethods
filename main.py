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
                print("{} ".format(element), end="")
            print()
    else:
        print("Empty matrix given")


#MAIN
filepath = "files/"
filename = "var_18_a.txt"
matrix_18_a = readMatrixFromFile(filepath + filename)
printMatrix(matrix_18_a)


