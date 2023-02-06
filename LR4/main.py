#coding=utf-8
from math import *


# Лабораторная работа № 4 на тему «Интерполирование. Среднеквадратичное приближение. Равномерное приближение»
# Выполнили Запорожченко Кирилл и Педаев Михаил (ФЗ-11)
# 2023 г.


#FUNCTIONS

#значение функции в точке x
def getFx(x, variant):
    if (variant == 0):
        return pow(3, x) + 3 * pow(x, 3) + 2
    else:
        print("Not done yet")

# def createSplitDiffTable(variant):
#     a = 1
#     b = 2
#     n = 5
#     h = (b - a) / n
#     xs = [0 for i in range(n)]
#     differences_table = [[] for i in range(n)]
#     while(a <= b):
#         xs.append(a)
#         differences_table[0].append(getFx(a, variant))
#         a += h
#     for i in range(1, n):
#         for j in range(n - i):
#             differences_table[i].append((differences_table[i - 1][j + 1] - differences_table[i - 1][j]) / (xs[j + i] - xs[j]))
#
#     return xs, differences_table
#


#MAIN
if __name__ == '__main__':
    variants = [0]
    for variant in variants:
        xs, difft = createSplitDiffTable(variant)
        for i in range(len(xs)):
            print(xs[i], difft[i])