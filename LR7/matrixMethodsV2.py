from math import sqrt

#Вектор в матрицу размерностью Nx1 и обратно

def Vector_to_List(vector):
    n = len(vector)
    new_vector = [[0] for i in range(n)]
    for i in range(n):
        new_vector[i][0] = vector[i]
    return new_vector

def Vector_from_List(vector):
    n = len(vector)
    new_vector = [0]*n
    for i in range(n):
        new_vector[i] = vector[i][0]
    return new_vector


#Преобразование вектора в матрицу и обратно

def Vector_to_Matrix(vector):
    if type(vector[0]) != list:
        return Vector_to_List(vector)
    else:
        return vector

def Vector_from_Matrix(vector):
    if len(vector[0]) == 1:
        return Vector_from_List(vector)
    else:
        print("Данный объект не является вектором")
        return vector


#Норма вектора

def Vector_norm(vector):
    if type(vector[0]) != list:
        vector = Vector_to_List(vector)
    if len(vector[0]) != 1:
        raise ValueError("Данный объект не является вектором")
    square = 0
    for i in vector:
        square += i[0]*i[0]
    result = sqrt(square)
    return result


#Сложение матриц и произведение матриц на константы

def Multiplicate(A,B,is_res_vec = False):
    if type(B[0]) != list:
        B = Vector_to_List(B)
    h = len(B)
    w = len(B[0])
    C = [[0]*w for i in range(h)]
    if type(A) != list:
        for i in range(h):
            for j in range(w):
                C[i][j] = A*B[i][j]
    else:
        if type(A[0]) != list:
            A = Vector_to_List(A)
        if h != len(A) or w != len(A[0]):
            raise ValueError("Матрицы разного размера")
        for i in range(h):
            for j in range(w):
                C[i][j] = A[i][j]*B[i][j]
    if is_res_vec:
        C = Vector_from_List(C)
    return C

def Addition(A,B,is_res_vec = False):
    if type(A[0]) != list:
        A = Vector_to_List(A)
    if type(B[0]) != list:
        B = Vector_to_List(B)
    h = len(A)
    w = len(A[0])
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise ValueError("Матрицы разного размера")
    C = [[0]*w for i in range(h)]
    for i in range(h):
        for j in range(w):
            C[i][j] = A[i][j] + B[i][j]
    if is_res_vec:
        C = Vector_from_List(C)
    return C