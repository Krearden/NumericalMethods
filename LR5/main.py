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
        pow(5, x) / log(5) - 3 * x * x + 3 * x
    elif (variant == 18):
        pass
    elif (variant == 27):
        pass

#формула Ньютона-Лейбница
def NewtonLeybnic(a, b):
    return f(b) - f(a)

#MAIN
if __name__ == '__main__':
    a = 1
    b = 2


