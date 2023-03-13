#coding=utf-8
from math import *


# Лабораторная работа № 5 на тему «Численное интегрирование»
# Выполнили Запорожченко Кирилл и Педаев Михаил (ФЗ-11)
# 2023 г.


# delete comment below if need to write to file
# file_path = 'pedaev_zaporozhchenko_LR6(5).txt'
# sys.stdout = open(file_path, "w")


#GLOBAL VARIABLES
variant = 18
count_uses = 0
count_uses_dx = 0


#FUNCTIONS

#значение фукнции
def f(x):
    global count_uses
    count_uses += 1
    if (variant == 8):
        return pow(5, x) - 6 * x + 3
    elif (variant == 18):
        return pow(3, x) - 2 * x + 5
    elif (variant == 27):
        return pow(e, -2 * x) - 2 * x + 1

#производная от фукнции
def dfdx(x):
    global count_uses_dx
    count_uses_dx += 1
    if (variant == 8):
        return log(5) * pow(5, x) - 6
    elif (variant == 18):
        return log(3) * pow(3, x) - 2
    elif (variant == 27):
        return -2 * pow(e, -2 * x) - 2

#первообразная от фукнции
def F(x):
    if (variant == 8):
        return pow(5, x) / log(5) - 3 * x * x + 3 * x
    elif (variant == 18):
        return pow(3, x) / log(3) - (x * x) + 5 * x
    elif (variant == 27):
        return -(x * x) + x - pow(e, -2 * x) / 2

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

#интеграл методом трапеций с модификацией сплайном
def spline_trapezoidalMethod(a, b, h):
    if a > b:
        a, b = b, a
    inner_sum = 0.0
    x_i = a + h
    while x_i < b:
        inner_sum += f(x_i)
        x_i += h

    return 0.5 * h * (f(a) + f(b) + 2 * inner_sum) + h * h * (dfdx(a) - dfdx(b)) / 12.0

#интеграл методом Симпсона
def simpsonMethod(a, b, h):
    if a > b:
        a, b = b, a
    sum = 0.0
    x_i = a
    x_i1 = a + h
    f_i = f(x_i)
    while (x_i1 <= b):
        f_i1 = f(x_i1)
        sum += (x_i1 - x_i) * (f_i + 4 * f((x_i + x_i1) / 2) + f_i1) / 6
        x_i = x_i1
        x_i1 += h
        f_i = f_i1

    return sum

def map_x(a, b, x):
    return 0.5 * (a + b + x * (b - a))

def map_c(a, b, c):
    return 0.5 * c * (b - a)

#Интеграл трехточечным методом Гаусса
def gaussMethod(a, b, h):
    # Проверка порядка границ интервала
    if a > b:
        a, b = b, a
    # Значения коэффициентов и узлов для метода Гаусса
    x = sqrt(3.0 / 5.0)
    c1 = 5.0 / 9.0
    c2 = 8.0 / 9.0
    # Сумма внутри интервала
    sum = 0.0
    xi = a
    while xi < b:
        # Вычисление узлов и коэффициентов на текущем интервале
        c_1 = map_c(xi, xi + h, c1)
        c_2 = map_c(xi, xi + h, c2)
        x_1 = map_x(xi, xi + h, x)
        x_2 = map_x(xi, xi + h, -x)
        # Вычисление значения интеграла на текущем интервале
        sum += c_1 * f(x_1) + c_2 * f(0.5 * h + xi) + c_1 * f(x_2)
        # Переход к следующему интервалу
        xi += h

    return sum


#Решение определенного интеграла с выбором метода, погрешность по Рунге
def solveIntegral(a, b, integral_sum_method, diff, alpha, epsilon):
    global count_uses
    global count_uses_dx
    n = 1
    h = get_h(a, b, n)
    I = 0
    I_before = 0
    I_2before = 0
    l = 1 / log(alpha)
    delta_r_coef = pow(alpha, diff) / (1.0 - pow(alpha, diff))
    print("|  N  |    h    | Интеграл     | Погрешность (отн.)  | Оценка погр. по Рунге     |  k     |")
    while (True):
        #обнуляем счетчики обращений к фукнции и производной
        if (integral_sum_method != gaussMethod):
            count_uses = 0
        count_uses_dx = 0
        #вычисляем интеграл выбранным численным методом
        I = integral_sum_method(a, b, h)
        #Вычисляем оценку порядка сходимости. Должно приближаться к diff - порядку сходимости.
        if (n != 1 and n != 2):
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
        absolut_error = abs(NewtonLeybnic(a, b) - I)
        delta = absolut_error / NewtonLeybnic(a, b)
        #вывод на экран
        delta_r = "{:.12f}".format(delta_r) if not isnan(delta_r) else "nan"
        k = "{:.4f}".format(k) if not isnan(k) else "nan"
        print(f"|{n:5d}|  {h:.5f}|   {I:.9f}|       {delta:.12f}|            {delta_r:>15}|  {k:>6}|")
        #условие выхода из цикла
        if (abs(delta) < epsilon and n >= 4):
            break
        #заполняем информацию о предыдущих вычислениях
        I_2before = I_before
        I_before = I
        #обновляем шаг
        h *= alpha
        n = ceil((b - a) / h)

    print(f"Результат = {I:.12f}")
    print(f"Количество обращений к фукнции = {count_uses}", end = "")
    if (count_uses_dx != 0):
        print(f" + {count_uses_dx}", end = "")


#MAIN
if __name__ == '__main__':
    a = 1
    b = 2
    epsilon = 1e-08
    alpha = 0.5
    print(f"Вариант {variant}")
    print(f"Точное значение интеграла J = {NewtonLeybnic(a, b):.16f}")
    print("\n")
    print("Метод трапеций")
    solveIntegral(a, b, trapezoidalMethod, 2, alpha, epsilon)
    print("\n")
    print("Метод трапеций модифицированный сплайном")
    solveIntegral(a, b, spline_trapezoidalMethod, 4, alpha, epsilon)
    print("\n")
    print("Метод Симпсона")
    solveIntegral(a, b, simpsonMethod, 4, alpha, epsilon)
    print("\n")
    print("Метод Гаусса (трехточечный)")
    count_uses = 0
    solveIntegral(a, b, gaussMethod, 6, alpha, epsilon)



