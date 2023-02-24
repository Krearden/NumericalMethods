#coding=utf-8
import math
from matrixMethods import *

# Лабораторная работа № 3 на тему «Решение нелинейных уравнений»
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

#get eucledian norm
def getEucledianNorm(A):
    n = len(A)
    epsilon = 1e-06
    A = multiplyMatrix(getTransposedMatrix(A), A)
    max_non_diag = abs(A[0][1])
    A_copy = [[A[i][j] for j in range(n)] for i in range(n)]
    while (max_non_diag > epsilon):
        max_non_diag = abs(A[0][1])
        I = 0
        J = 1
        for i in range(n):
            for j in range(i + 1, n):
                if (abs(A[i][j]) > max_non_diag):
                    max_non_diag = abs(A[i][j])
                    I = i
                    J = j
        theta = math.atan(2 * A[I][J] / (A[I][I] - A[J][J])) / 2
        c = math.cos(theta)
        s = math.sin(theta)
        for m in range(n):
            if (m != I and m != J):
                A_copy[m][I] = c * A[I][m] + s * A[J][m]
                A_copy[m][J] = -s * A[I][m] + c * A[J][m]
                A_copy[I][m] = c * A[I][m] + s * A[J][m]
                A_copy[J][m] = -s * A[I][m] + c * A[J][m]
        A_copy[I][I] = c * c * A[I][I] + 2 * s * c * A[I][J] + s * s * A[J][J]
        A_copy[J][J] = s * s * A[I][I] - 2 * s * c * A[I][J] + c * c * A[J][J]
        A_copy[I][J] = (c * c - s * s) * A[I][J] + s * c * (A[J][J] - A[I][I])
        A_copy[J][I] = (c * c - s * s) * A[I][J] + s * c * (A[J][J] - A[I][I])
        A = [[A_copy[i][j] for j in range(n)] for i in range(n)]
    max_diagonal = 0
    for i in range(n):
        if (abs(A[i][i] > max_diagonal)):
            max_diagonal = abs(A[i][i])

    return math.sqrt(max_diagonal)



#Jacobi marix
def getJacobian(x, y, variant):
    if (variant == 18):
        jak = [
        [0.0, -1*math.sin(y+1)],
        [-1*math.cos(x), 0.0 ]
        ]
    elif (variant == 27):
        jak = [
        [0.0, -1*math.sin(y+0.5)],
        [math.cos(x)/2,  0.0 ]
        ]
    else:
        jak == [[0, 0], [0, 0]]
    return jak



def getFxy(XY, variant):
    fX = [0 for i in range(len(XY))]
    if (variant == 18):
        fX[0] = 2 * XY[0] - math.cos(XY[1] + 1)
        fX[1] = XY[1] + math.sin(XY[0]) + 0.4
    elif (variant == 27):
        fX[0] = math.cos(XY[1] + 0.5) - XY[0] - 0.8
        fX[1] = math.sin(XY[0]) - 2 * XY[1] - 1.6
    else:
        print("Такого еще нет варианта")

    return fX


#Derivative matrix - матрица производных для систем в неявном виде
def getDerivativeMatrix(variant, x, y):
    N = 2
    dF = [[0 for i in range(N)] for i in range(N)]
    if (variant == 18):
        dF[0][0] = 2  # f1 over x
        dF[0][1] = math.cos(XY[0])  # f2 over x
        dF[1][0] = math.sin(XY[1] + 1)  # f1 over y
        dF[1][1] = 1  # f2 over y
    elif (variant == 27):
        dF[0][0] = -1 # f1 over x
        dF[0][1] = math.cos(XY[0]) # f2 over x
        dF[1][0] = -math.sin(XY[1] + 0.5) # f1 over y
        dF[1][1] = -2 # f2 over y
    else:
        print("Пока здесь такого варианта нету")

    return dF


#length() for newton method to stop
def length(XY, variant):
    Fxy = getFxy(XY, variant)
    sqrt = math.sqrt(Fxy[0] * Fxy[0] + Fxy[1] * Fxy[1])
    return sqrt

#len() for others to calculate pogreshnost'
def leng(XY):
    return math.sqrt(XY[0] * XY[0] + XY[1] * XY[1])


def getDeterminant(matrix):
    if (len(matrix) != 2):
        print("Wrong matrix size")
    else:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def printHeader(method_name):
    if (method_name == "newton"):
        print("")
        print("+-NEWTON'S METHOD-----------------------------------------------------------------------------------+")
        print("| Itr |      X      |      Y      |    Норма невязки    |          F1         |          F2         |")
        print("+-----+-------------+-------------+---------------------+---------------------+---------------------+")

    elif (method_name == "iteration"):
        print("")
        print("+-SIMPLE ITERATION METHOD-------------------------------------------------------------------------------+")
        print("| Itr |      X      |      Y      | Норма невязки | Погрешность | Оценка погрешности |  Норма якобиана  |")
        print("+-----+-------------+-------------+---------------+-------------+--------------------+------------------+")


    elif (method_name == "gradient"):
        print("")
        print("+-GRADIENT METHOD---+-------------+---------+--------------------+--------------------+-----+")
        print("| Itr |      X      |      Y      |  Alpha  |    Норма невязки   |     Погрешность    |  k  |")
        print("+-----+-------------+-------------+---------+--------------------+--------------------+-----+")

def writeHeader(method_name, output_file):
    if (method_name == "newton"):
        output_file.write("\n")
        output_file.write("\n+-NEWTON'S METHOD-----------------------------------------------------------------------------------+")
        output_file.write("\n| Itr |      X      |      Y      |    Норма невязки    |          F1         |          F2         |")
        output_file.write("\n+-----+-------------+-------------+---------------------+---------------------+---------------------+")

    elif (method_name == "iteration"):
        output_file.write("\n")
        output_file.write("\n+-SIMPLE ITERATION METHOD-------------------------------------------------------------------------------+")
        output_file.write("\n| Itr |      X      |      Y      | Норма невязки | Погрешность | Оценка погрешности |  Норма якобиана  |")
        output_file.write("\n+-----+-------------+-------------+---------------+-------------+--------------------+------------------+")
        pass
    elif (method_name == "gradient"):
        output_file.write("\n")
        output_file.write("\n+-GRADIENT METHOD-------------------------------------------------------------------+-----+")
        output_file.write("\n| Itr |      X     |      Y     |  Alpha  |    Норма невязки   |     Погрешность    |  k  |")
        output_file.write("\n+-----+------------+------------+---------+--------------------+--------------------+-----+")

#Newton's method
def methodNewton(XY, variant, output_file):
    epsilon = 1e-12
    iteration_counter = 0
    while (length(XY, variant) > epsilon):
        iteration_counter += 1
        #матрица производных
        dF = getDerivativeMatrix(variant, XY[0], XY[1])
        #определители
        determinant = getDeterminant(dF)
        #fx
        Fx = getFxy(XY, variant)
        n = len(dF)
        #for x
        dF_x = [[dF[j][i] for j in range(n)] for i in range(n)]
        dF_x[0][0] = Fx[0]
        dF_x[1][0] = Fx[1]
        determinant_x = getDeterminant(dF_x)
        #for y
        dF_y = [[dF[j][i] for j in range(n)] for i in range(n)]
        dF_y[0][1] = Fx[0]
        dF_y[1][1] = Fx[1]
        determinant_y = getDeterminant(dF_y)
        #find new x and y
        newXY = [0 for i in range(len(dF))]
        newXY[0] = XY[0] - determinant_x / determinant #find new XY
        newXY[1] = XY[1] - determinant_y / determinant #find new Y
        #copy
        XY = [newXY[i] for i in range(len(dF))]
        #printinfo
        print("| {:3} |  {:3.8f} | {: 3.8f} |          {:1.8f} |         {: 3.8f} |         {: 3.8f} |".format(iteration_counter, XY[0], XY[1], length(XY, variant), getFxy(XY, variant)[0], getFxy(XY, variant)[1]))
        output_file.write("\n| {:3} |  {:3.8f} |  {:3.8f} |          {:1.8f} |         {:3.8f} |         {:3.8f} |".format(iteration_counter, XY[0], XY[1], length(XY, variant), getFxy(XY, variant)[0], getFxy(XY, variant)[1]))
    return XY


#func to find min
def minimumFunction(XY, variant):
    Fxy = getFxy(XY, variant) # Fx[0] - F1, Fx[1] - F2
    return math.pow(Fxy[0], 2) + math.pow(Fxy[1], 2)


def minimumFunctionDXDY(XY, variant):
    Fxy = getFxy(XY, variant) # Fxy[0] - F1, Fxy[1] - F2
    dF = getDerivativeMatrix(variant, XY[0], XY[1])
    dx = 2 * Fxy[0] * dF[0][0] + 2 * Fxy[1] * dF[0][1]
    dy = 2 * Fxy[0] * dF[1][0] + 2 * Fxy[1] * dF[1][1]
    return [dx, dy]



#get eucledian norm of vector
def getEucledianVectorNorm(vector):
    result = 0
    for element in vector:
        result += element * element
    return math.sqrt(result)


#метод Градиентного Спуска не работает - проблема в условии завершения
def methodGradient(XY, variant, newtonXY, output_file):
    epsilon = 1e-04
    timetostop = 1 + epsilon
    iteration_counter = 0
    while (timetostop > epsilon):
        iteration_counter += 1
        alpha = 1
        lambdaa = 0.5
        k = -1
        min_func_dxdy = minimumFunctionDXDY(XY,variant)
        #вычисляем значение переменной alpha адаптивным методом
        while (minimumFunction([XY[0] - alpha * min_func_dxdy[0], XY[1] - alpha * min_func_dxdy[1]], variant) > minimumFunction(XY, variant)):
            alpha *= lambdaa
            k += 1
        # find new x and y
        newXY = [0 for i in range(len(XY))]
        newXY[0] = XY[0] - alpha * min_func_dxdy[0] #new X
        newXY[1] = XY[1] - alpha * min_func_dxdy[1] #new Y
        #copy
        XY = [newXY[i] for i in range(len(XY))]
        timetostop = getEucledianVectorNorm(getFxy(XY, variant))
        #printinfo
        print("| {:3} | {: 3.8f} | {:3.8f} | {:3.5f} | {:3.16f} | {:3.16f} | {:3} |".format(iteration_counter, XY[0], XY[1], alpha, timetostop, leng([XY[0] - newtonXY[0], XY[1] - newtonXY[1]]), k))
        output_file.write("\n| {:3} | {:3.8f} | {:3.8f} | {:3.5f} | {:3.16f} | {:3.16f} | {:3} |".format(iteration_counter, XY[0], XY[1], alpha, timetostop, leng([XY[0] - newtonXY[0], XY[1] - newtonXY[1]]), k))

    return XY

def Iteration(x, y, newtonXY, variant):
    round = 0
    eps = 1e-4

    while(length([x, y], variant) > eps):
        round+=1

        if (variant == 18):
            newX = math.cos(y+1)/2
            newY = -math.sin(x) - 0.4
        elif (variant == 27):
            newX = -0.8 + math.cos(y + 0.5)
            newY = (math.sin(x) - 1.6)/2

        mJac = getJacobian(x, y, variant)
        normJac = getEucledianNorm(mJac)
        errorRate = normJac / (1 - normJac) * leng([x - newtonXY[0], y - newtonXY[1]])

        
        x = newX
        y = newY
        Fxy = getFxy([x, y], variant)

        print("| {:3} | {: 3.8f} | {: 3.8f} |    {:1.8f} |  {:1.8f} |         {:1.8f} |      {: 3.8f} |".format(round, x, y, length([x, y], variant), leng([x - newtonXY[0], y - newtonXY[1]]), errorRate, normJac))
        output_file.write("\n| {:3} | {: 3.8f} | {: 3.8f} |    {:1.8f} |  {:1.8f} |         {:1.8f} |      {: 3.8f} |".format(round, x, y, length([x, y], variant), leng([x - newtonXY[0], y - newtonXY[1]]), errorRate, normJac))
    return x,y


#MAIN
if __name__ == "__main__":
    variants = [18, 27]
    output_filename = "lr3_output.txt"
    output_file = open(output_filename, "w")
    output_file.write("Запорожченко, Педаев. ЛР3.")

    for i in range(len(variants)):
        if (variants[i] == 18):
            XY = [0.5, -0.9]
        elif (variants[i] == 27):
            XY = [0.2, -0.7]

        print("\n\n\nVARIANT {}, {}".format(variants[i], XY))
        output_file.write("\n\nVARIANT {}".format(variants[i]))

        print("\nx0={} y0={}".format(XY[0], XY[1]))
        output_file.write("\n\nx0={} y0={}".format(XY[0], XY[1]))
        print("Alfa0= 1.000 Lambda= 0.500")
        output_file.write("\nAlfa0= 1.000 Lambda= 0.500\n")

        print("\nМатрица производных")
        output_file.write("\n\nМатрица производных\n")

        if (variants[i] == 18):
            print("2;  cos(x);\nsin(y+1);  1;")
            output_file.write("2;  cos(x);\nsin(y+1);  1;")
        elif (variants[i] == 27):
            print("-1;  cos(x);\n-sin(y+0.5);  -2;")
            output_file.write("-1;  cos(x);\n-sin(y+0.5);  -2;")

        printHeader("newton")
        writeHeader("newton", output_file)
        newtonXY = methodNewton(XY, variants[i], output_file)
        print("+-----+-------------+-------------+-----------------------------------------------------------------+")
        print(f"x = {newtonXY[0]}, y = {newtonXY[1]};")
        output_file.write("\n+---------------------------------------------------------------------------------------------------+")
        output_file.write(f"\nx = {newtonXY[0]}, y = {newtonXY[1]};")


        print("\nМетод простой итерации")
        output_file.write("\n\nМетод простой итерации")

        if (variants[i]==18):
            print("Fi1(x,y)=cos(y+1)/2")
            print("Fi2(x,y)=-sin(x)-0.4")
            output_file.write("\nFi1(x,y)=cos(y+1)/2\nFi2(x,y)=-sin(x)-0.4")
        elif(variants[i]==27):
            print("Fi1(x,y)=-0.8+cos(y+0.5)")
            print("Fi2(x,y)=(sin(x)-1.6)/2")
            output_file.write("\nFi1(x,y)=-0.8+cos(y+0.5)\nFi2(x,y)=(sin(x)-1.6)/2")
        
        print("\nЯкобиан")
        output_file.write("\n\nЯкобиан\n")
        if (variants[i]==18):
            print("0;  -sin(y+1);")
            print("-cos(x);  0;")
            output_file.write("0;  -sin(y+0.5);\ncos(x);  0;")
        elif(variants[i]==27):
            print("0;  -sin(y+0.5);")
            print("cos(x);  0;")
            output_file.write("0;  -sin(y+0.5);\ncos(x);  0;")
        
        print("\nЗначение")
        output_file.write("\n\nЗначение\n")
        jak = getJacobian(XY[0], XY[1], variants[i])
        print("{: 2.4f}  {: 2.4f}".format(jak[0][1], jak[0][1]))
        print("{: 2.4f}  {: 2.4f}".format(jak[1][1], jak[1][1]))
        output_file.write("{: 2.4f}  {: 2.4f}\n{: 2.4f}  {: 2.4f}".format(jak[0][1], jak[0][1], jak[1][1], jak[1][1]))

        print("Норма={:2.4f}".format(getEucledianNorm(jak)))
        output_file.write("\nНорма={:2.4f}".format(getEucledianNorm(jak)))

        printHeader("iteration")
        writeHeader("iteration", output_file)
        itX, itY = Iteration(XY[0], XY[1], newtonXY, variants[i])
        print("+-------------------------------------------------------------------------------------------------------+")
        print(f"x = {itX}, y = {itY};")
        output_file.write("\n+-------------------------------------------------------------------------------------------------------+")
        output_file.write(f"\nx = {itX}, y = {itY};")

        printHeader("gradient")
        writeHeader("gradient", output_file)
        gXY = methodGradient(XY, variants[i], newtonXY, output_file)
        print("+-----+-------------+-------------+---------+--------------------+--------------------+-----+")
        output_file.write("\n+-----+------------+------------+---------+--------------------+--------------------+-----+")
        print(f"x = {gXY[0]}, y = {gXY[1]};")
        output_file.write(f"\nx = {gXY[0]}, y = {gXY[1]};")
        #write to file

