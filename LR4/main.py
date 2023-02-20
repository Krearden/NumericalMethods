#coding=utf-8
from math import *
from LU import *


# Лабораторная работа № 4 на тему «Интерполирование. Среднеквадратичное приближение. Равномерное приближение»
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

#значение функции в точке x
def getFx(x, variant):
    if (variant == 0):
        return pow(3, x) + 3 * pow(x, 3) + 2
    else:
        print("Not done yet")

#Расстояние h между точками интерполяции по трем параметрам
def getH(a, b, n):
    return (b - a) / (n - 1)

#Факториал числа
def Factorial(n):
    if (n == 0):
        return 1
    else:
        return n * Factorial(n - 1)

#Таблица разделенных разностей
def createSplitDiffTable(variant, a, b, n):
    h = getH(a, b, n)
    xs = []
    differences_table = [[] for i in range(n)]
    while(a <= b):
        xs.append(a)
        differences_table[0].append(getFx(a, variant))
        a += h
        a = round(a, 2)

    for i in range(1, n):
        for j in range(n - i):
            differences_table[i].append((differences_table[i - 1][j + 1] - differences_table[i - 1][j]) / (xs[j + i] - xs[j]))

    return xs, differences_table

#Печать таблицы разделенных разностей
def printSplitDiffTable(xs, differences_table):
    print("Таблица разделенных разностей:")
    for i in range(len(xs)):
        print("{:2.1f}".format(xs[i]), end = '  ')
        for j in range(len(differences_table[i])):
            print(differences_table[j][i], end = ' ')
        print()

#Полином ньютона по таблице разделенных разностей
def getNewtonPolynom(x, xs, differences_table):
    n = len(xs)
    y = 0
    part_summ = 1
    for i in range(n):
        y += part_summ * differences_table[i][0]
        part_summ *= x - xs[i]
    return y

#Фукнция. Возвращает произведения разностей (x - xi): i = 0 ... n - 1; xi берется из таблицы разделенных разностей.
def getOmega(x, xs, n):
    result = 1
    for i in range(n):
        result *= x - xs[i]
    return abs(result)

#Фукнция - возвращает производную функции n-ного порядка
def Derivative(variant, x, n = 1):
    if (variant == 0):
        if (n >= 4):
            return pow(3, x) * pow(log(3), n)
        elif (n == 1):
            return log(3) * pow(3, x) + 9 * x * x
        else:
            return
            # return log(3) * pow(3, x) + 9 * x * x
    else:
        print("Not implemented yet")

def getMaxDerivatire(variant, xs, n):
    maximum_derivative = 1e-100
    for x in xs:
        first_derivative = abs(Derivative(variant, x, n))
        if (first_derivative > maximum_derivative):
            maximum_derivative = first_derivative
    return maximum_derivative

def getNewtonError(variant, x, xs, n):
    return getMaxDerivatire(variant, xs, n) * getOmega(x, xs, n)/ Factorial(n)


def printNewtonInterpolation(a, b, n, xs, differences_table):
    print("Интерполяционная формула Ньютона")
    print("M6 (derivative) = {}\n".format(getMaxDerivatire(variant, xs, 6)))
    h = getH(a, b, n)
    x = a + h / 2
    print("x       f(x)       Pn(x)       Delta       Оценка  ")
    while (x <= b):
        Fx = getFx(x, variant)
        newtonPx = getNewtonPolynom(x, xs, differences_table)
        print("{:1.2f}  {:2.6f}  {:2.6f}  {:2.6e}  {:2.6e}".format(x, Fx, newtonPx, abs(newtonPx - Fx), getNewtonError(variant, x, xs, n)))
        x += h



#solve SLAU
def solve_SLAU(matrix, b):
    A, L, U, p = LU(matrix)
    # find x
    b = multipMatrixByVector(getPmatrix(p), b)
    y = forwardSubstitution(L, b)
    x = backwardSubstitution(U, y)
    return x

def fi0(tau):
    return (1 + 2 * tau) * pow(1 - tau, 2)

def fi1(tau):
    return tau * pow(1 - tau, 2)

#S31 - кубический сплайн дефекта 1
def S31(x, xi, xnext, mi, mnext, h):
    tau = (x - xi) / h
    return fi0(tau) * getFx(xi, variant) + fi0(1 - tau) * getFx(xnext, variant) + h * (fi1(tau) * mi - fi1(1 - tau) * mnext)

def printCubicInterpolation(a, b, n, variant):
    print("Интерполяция кубическим сплайном деф. 1")
    #задаем трехдиагональную матрицу для выбранных точек интерполяции
    matrix = [[4, 1, 0, 0],[1, 4, 1, 0],[0, 1, 4, 1],[0, 0, 1, 4]]
    #задаем вектор 1-х пр-х
    m = [0 for i in range(n)]
    m[0] = Derivative(variant, a)
    m[5] = Derivative(variant, b)
    #задаем вектор правых частей
    h = getH(a, b, n)
    b_vector = [0 for i in range(n - 2)]
    for i in range(n - 2):
        b_vector[i] = 3 * ( getFx(xs[i + 2], variant) - getFx(xs[i], variant) ) / h
    #выполняем корректировку крайних элементов
    b_vector[0] -= m[0]
    b_vector[3] -= m[5]
    #решение СЛАУ
    mm = solve_SLAU(matrix, b_vector)
    for i in range(1, n - 1):
        m[i] = mm[i - 1]
    #print
    print()
    print("M5 (derivative) = {}\n".format(getMaxDerivatire(variant, xs, 5)))
    print("x[i]  df/dx(x[i])    m[i]      Delta       Оценка")
    for i in range(n):
        print("{:1.2f}  {:2.6f}  {:2.6f}  {:.6e}  {:.6e}".format(xs[i], Derivative(variant, xs[i]), m[i], abs(m[i] - Derivative(variant, xs[i])), getMaxDerivatire(variant, xs, 5) * pow(h, 4) / 60.0))

    print("\nM4 (derivative) = {}\n".format(getMaxDerivatire(variant, xs, 4)))
    print("x       f(x)     S31(f;x)     Abs(f - S31)        Оценка")
    x = a + h / 2
    i = 0
    while(x <= b):
        s31 = S31(x, xs[i], xs[i + 1], m[i], m[i + 1], h)
        fx = getFx(x, variant)
        error =  pow(h, 4) * (getMaxDerivatire(variant, xs, 4) / 384.0 + getMaxDerivatire(variant, xs, 5) * h / 240.0)
        print("{:1.2f}  {:2.6f}  {:2.6f}  {:1.12e}  {:1.12e}".format(x, fx, s31, abs(fx - s31), error))
        x += h
        i += 1




#MAIN
if __name__ == '__main__':
    variants = [0]
    for variant in variants:
        a = 1
        b = 2
        n = 6
        xs, differences_table = createSplitDiffTable(variant, a, b, n)
        printSplitDiffTable(xs, differences_table)
        print()
        printNewtonInterpolation(a, b, n, xs, differences_table)

        print()
        printCubicInterpolation(a, b, n, variant)