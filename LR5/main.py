#coding=utf-8
from math import *


# Лабораторная работа № 5 на тему «Численное интегрирование»
# Выполнили Запорожченко Кирилл и Педаев Михаил (ФЗ-11)
# 2023 г.


# delete comment below if need to write to file
# file_path = 'pedaev_zaporozhchenko_LR6(5).txt'
# sys.stdout = open(file_path, "w")


#GLOBAL VARIABLES
variant = 8
count_uses = 0


#FUNCTIONS

#значение фукнции
def f(x):
    global count_uses
    count_uses += 1
    if (variant == 8):
        return pow(5, x) - 6 * x + 3
    elif (variant == 18):
        pass
    elif (variant == 27):
        pass

#производная от фукнции
def dfdx(x):
    global count_uses
    count_uses += 1
    if (variant == 8):
        return log(5) * pow(5, x) - 6
    elif (variant == 18):
        pass
    elif (variant == 27):
        pass

#первообразная от фукнции
def F(x):
    if (variant == 8):
        return pow(5, x) / log(5) - 3 * x * x + 3 * x
    elif (variant == 18):
        pass
    elif (variant == 27):
        pass

#формула Ньютона-Лейбница
def NewtonLeybnic(a, b):
    return F(b) - F(a)

#шаг
def get_h(a, b, n):
    return (b - a) / n

#интеграл методом трапеций
def trapezoidalMethod(a, b, h):
    if a > b:
        a, b = b, a
    inner_sum = 0.0
    x_i = a + h
    while x_i < b:
        inner_sum += f(x_i)
        x_i += h

    return 0.5 * h * (f(a) + f(b) + 2 * inner_sum)

#Решение определенного интеграла с выбором метода, погрешность по Рунге
def solveIntegral(a, b, integral_sum_method, diff, alpha, epsilon):
    global count_uses
    n = 1
    h = get_h(a, b, n)
    I = 0
    I_before = 0
    I_2before = 0
    l = 1 / log(alpha)
    delta_r_coef = pow(alpha, diff) / (1.0 - pow(alpha, diff))
    print("|  N  |    h    | Интеграл     | Погрешность         | Оценка погр. по Рунге     |  k     |")
    while (True):
        #обнуляем счетчик обращений к фукнции
        count_uses = 0
        I = integral_sum_method(a, b, h)
        #Вычисляем оценку порядка сходимости. Должно приближаться к diff - порядку сходимости.
        if (n != 1 and ((I - I_2before) / (I_before - I_2before) - 1.0) > 0):
            k = l * log((I - I_2before) / (I_before - I_2before) - 1.0)
        else:
            k = nan
        #Оценка погрешности по правилу Рунге
        #Проверяем на первую итерацию, так как на первой итерации
        #Еще не существует прошлого значения интеграла
        if (n != 1):
            delta_r = delta_r_coef * (I - I_before)
        else:
            delta_r = nan
        #Относительная поргрешность
        delta = abs(NewtonLeybnic(a, b) - I) / NewtonLeybnic(a, b)
        #вывод на экран
        delta_r = "{:.12f}".format(delta_r) if not isnan(delta_r) else "nan"
        k = "{:.4f}".format(k) if not isnan(k) else "nan"
        print(f"|{n:5d}|  {h:.5f}|   {I:.9f}|       {delta:.12f}|            {delta_r:>15}|  {k:>6}|")
        #условие выхода из цикла
        if (delta < epsilon):
            break
        #заполняем информацию о предыдущих вычислениях
        I_2before = I_before
        I_before = I
        #обновляем шаг
        h *= alpha
        n = ceil((b - a) / h)
    print(f"Результат = {I:.12f}")
    print(f"Количество обращений к фукнции = {count_uses}")

#MAIN
if __name__ == '__main__':
    a = 1
    b = 2
    epsilon = 1e-08
    alpha = 0.5
    print(f"Вариант {variant}")
    print(f"Точное значение интеграла J = {NewtonLeybnic(a, b):.16f}")
    print()
    print("Метод трапеций")
    solveIntegral(a, b, trapezoidalMethod, 2, alpha, epsilon)



