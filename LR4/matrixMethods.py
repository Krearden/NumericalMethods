#multiply Matrix By Vector (to find b)
def multipMatrixByVector(m, v):
    rv = []
    for row in m:
        rowResult = 0
        for elemPos in range(len(row)):
            rowResult += row[elemPos] * v[elemPos]
        rv.append(rowResult)

    return rv


#multiply Matrix by Matrix
def multiplyMatrix(first_matrix, second_matrix):
    length = len(first_matrix) 
    result_matrix = [[0 for i in range(length)] for i in range(length)]
    for i in range(length):
        for j in range(length):
            for k in range(length):
                result_matrix[i][j] += first_matrix[i][k] * second_matrix[k][j]

    return result_matrix


#substract Matrix from Matrix
def matrixSubtraction(first_matrix, second_matrix):
    length = len(first_matrix)
    result_matrix = [[0 for i in range(length)] for i in range(length)]
    for i in range(length):
        for j in range(length):
            result_matrix[i][j] = first_matrix[i][j] - second_matrix[i][j]
    
    return result_matrix


def vectorSubstraction(vector_one, vector_two):
    n = len(vector_one)
    if (len(vector_one) == len(vector_two)):
        result = [0 for i in range(n)]
        for i in range(len(vector_one)):
            result[i] = vector_one[i] - vector_two[i]
        return result
    else:
        return


def multiplyMatrixByNumber(first_matrix, number):
    pass


def MatrixPlusMatrix(first_matrix, second_matrix):
    pass


#get transposed matrix
def getTransposedMatrix(A):
    n = len(A)
    A_transposed = [[A[j][i] for j in range(n)] for i in range(n)]
    return A_transposed


#compute determinant for 2x2 matrix
def computeDeterminant(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


if __name__ == "__main__":
    first_matrix = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ]
    second_matrix = [
        [1, 4, 7],
        [2, 5, 8],
        [3, 6, 9]
    ]

    print(multiplyMatrix(first_matrix, second_matrix))