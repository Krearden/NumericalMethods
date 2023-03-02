#coding=utf-8
from math import *
from LU import *
import sys


# Лабораторная работа № 4 на тему «Интерполирование. Среднеквадратичное приближение. Равномерное приближение»
# Выполнили Запорожченко Кирилл и Педаев Михаил (ФЗ-11)
# 2023 г.


# delete comment below if need to write to file
# file_path = 'pedaev_zaporozhchenko_LR5.txt'
# sys.stdout = open(file_path, "w")


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




#базисные фукнкции для среднеквадратичного приближения
def getBasisFuncitons(x):
    basis = []
    basis.append(1) #b1
    basis.append(x) #b2
    basis.append(x * x) #b3

    return basis

#get (gi, gj) - скалярное произведение базисных фукнций
def get_Gi_Gj(i, j, xs):
    sum_multiply = 0
    for x in xs:
        sum_multiply += getBasisFuncitons(x)[i] * getBasisFuncitons(x)[j]

    return sum_multiply

#get (f, g) - скалярное произведение базис. фукнции. и ф-ии f(x)
def get_f_g(i, xs, variant):
    sum_multiply = 0
    for x in xs:
        sum_multiply += getBasisFuncitons(x)[i] * getFx(x, variant)
    return sum_multiply

#returns b_i * b_j
def getBasisFunctionForContiniousVariant(x, i, j):
    return getBasisFuncitons(x)[i] * getBasisFuncitons(x)[j]

#интерграл методом Симпсона для определения коэффициентов матрицы непрерывного варианта
def simpson_method(i, j, N = 10000, a = 1, b = 2):
    h = (b - a) / N
    f_x0 = getBasisFunctionForContiniousVariant(a, i, j)
    sum_even = 0
    sum_uneven = 0
    for k in range(1, N, 2):
        sum_even += getBasisFunctionForContiniousVariant((a + k * h), i, j)
        sum_uneven += getBasisFunctionForContiniousVariant(a + (k + 1) * h, i, j)

    return h / 3 * (f_x0 + 2 * sum_even + 4 * sum_uneven)

#вектор правых частей для непрервыного варианта, посчитано с помощью WolframAlpha
def get_b_vector_continious(variant):
    if (variant == 0):
        return [18.711435359761024361685441, 30.2823757012612226379606520, 50.3984832581529801104382368632604486]
    else:
        print("No variant like this yuet MotherFU@@'sdfer")

#ф-я F(X)^2 - P(X)^2 для определения нормы погр. непр. варианта
def FxPxSquaresDifference(x, coefficients_array):
    return getFx(x, variant) ** 2 - (coefficients_array[0] + coefficients_array[1] * x + coefficients_array[2] * x ** 2) ** 2

#вычисление нормы погрешности непрерывного интерграла.
#численно вычисляется интеграл методом трапеций
def getErrorContinious(coefficients_array):
    ans = 0
    eps = 1e-6
    x = a + eps
    while x <= b:
        vl1 = FxPxSquaresDifference(x, coefficients_array)
        vl2 = FxPxSquaresDifference(x - eps, coefficients_array)
        ans += (vl1 + vl2) / 2. * eps
        x += eps
    return abs(ans)

#среднеквадратичное приближение
def printMiddleSquareApproximation(xs):
    print("\nСреднеквадратичное приближение")
    #матрица коэффициентов n*n (n - amount of basis functions)
    #содержит коэффициенты среднеквадратичного приближения функции
    matrix = [[0 for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            matrix[i][j] = get_Gi_Gj(i, j, xs)
    print("\nМатрица")
    printMatrix(matrix)
    #вектор правых частей
    b = [0 for i in range(3)]
    for i in range(3):
        b[i] = get_f_g(i, xs, variant)
    print("\nВектор правых частей")
    print(b)
    #получаем массив искомых коэф. путем решения СЛАУ
    coefficients_array = solve_SLAU(matrix, b)
    #печать полинома второй степени на экран
    print(f"\nP2(x) = ({coefficients_array[0]}) + ({coefficients_array[1]}) * x + ({coefficients_array[2]}) * x^2")
    #норма погрешности
    F_square = 0
    P_square = 0
    for x in xs:
        F_square += pow(getFx(x, variant), 2)
        P_square += pow(coefficients_array[0] + coefficients_array[1] * x + coefficients_array[2] * x * x, 2)
    print(f"Норма погрешности: {sqrt(F_square - P_square)}")
    print("x       Погрешность")
    for x in xs:
        p2 = coefficients_array[0] + coefficients_array[1] * x + coefficients_array[2] * x * x
        Fx = getFx(x, variant)
        print("{:1.1f}   {:1.12f}".format(x, p2 - Fx))
    print()
    print("Непрерывный вариант")
    # #вычисляем матрицу коэффициентов
    # for i in range(3):
    #     for j in range(3):
    #         matrix[i][j] = simpson_method(i, j)
    matrix = [[1, 1.5, 2.3333333333333333333], [1.5, 2.333333333333333, 3.75], [2.3333333333333333333, 3.75, 6.2]]
    print("Матрица: ")
    printMatrix(matrix)
    print("Вектор правых частей: ")
    print(b)
    b = get_b_vector_continious(variant)
    coefficients_array = solve_SLAU(matrix, b)
    # печать полинома второй степени на экран
    print(f"\nP2(x) = ({coefficients_array[0]}) + ({coefficients_array[1]}) * x + ({coefficients_array[2]}) * x^2")
    error_norm = getErrorContinious(coefficients_array)
    print("Норма погрешности: {}".format(error_norm))
    print("x       Погрешность")
    for x in xs:
        p2 = coefficients_array[0] + coefficients_array[1] * x + coefficients_array[2] * x * x
        Fx = getFx(x, variant)
        print("{:1.1f}   {:1.12f}".format(x, p2 - Fx))


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

        print()
        printMiddleSquareApproximation(xs)