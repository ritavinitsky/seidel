
matrixA = [[3,-1,1],[0,1,-1],[1,1,-2]]
vectorB = [4,-1,-3]

epsilon = 0.00001
numOfIterations = 100

def isMatrix(m):
    """
    Returns True if m is matrix and False otherwise
    :param m: parameter to be checked if it is matrix
    :return: True if m is matrix and False otherwise
    """
    try:
        # check if input are matrixes
        iter(m)
        if len(m) < 1:
            return False  # matrix is empty
        matcolcount = len(m[0])
        for i in m:
            iter(i)
            if len(i) != matcolcount:
                return False  # columns of matrix have different number of elements
    except TypeError as e:
        return False  # input is not matrixes
    return True


def makeDominantMatrix(matrix,resvec):
    """
    Tryes to make dominant matrix by swapping lines, returns the dominant matrix and vector, if success, if not success return 2 empty lists
    :param matrix: matrix to make dominant
    :param resvec: solution vector of matrix
    :return: the dominant matrix and vector, if success, if not success return 2 empty lists
    """
    matrix = matrix.copy()
    resvec = resvec.copy()
    i = 0
    while i < len(matrix):
        if not isDominantRow(matrix[i],i):
            j = 0
            flag = 0
            while j < len(matrix):
                if isDominantRow(matrix[j],i):
                    temp = matrix[j].copy()
                    matrix[j] = matrix[i].copy()
                    matrix[i] = temp
                    temp = resvec[j]
                    resvec[j] = resvec[i]
                    resvec[i] = temp
                    flag = 1
                    break
                j += 1
            if flag == 0:
                return [],[] #can't make dominant matrix
        i += 1
    if not isDominantMatrix(matrix):
        return [],[]
    return matrix,resvec

def buildIterationMatrix(matrix,resvec):
    """
    Returns iteration matrix built from matrix and vector resvec
    :param matrix: matrix from which to build iteration matrix
    :param resvec: vector from which to build iteration matrix
    :return: iteration matrix built from matrix and vector resvec
    """
    matrix = matrix.copy()
    resvec = resvec.copy()
    i = 0
    iterationMat = []
    while i < len(matrix):
        if matrix[i][i] == 0:
            s = i + 1
            flag = 0
            while s < len(matrix):
                if matrix[s][i] != 0:
                    # swap lines
                    temp = matrix[s].copy()
                    matrix[s] = matrix[i].copy()
                    matrix[i] = temp
                    flag = 1
                    temp = resvec[s]
                    resvec[s] = resvec[i]
                    resvec[i] = temp
                    break
                s += 1
            if flag == 0:
                return []  # can't build iterationMatrix
        iterationMat.append([])
        iterationMat[i].append(resvec[i] / matrix[i][i])
        j = 0
        while j < len(matrix):
            if i != j:
                iterationMat[i].append(-matrix[i][j] / matrix[i][i])
            else:
                iterationMat[i].append(0)
            j += 1
        i += 1
    return iterationMat

def solutionOfMatrixByYacobi(matrix, resvec):
    """
    Prints the solution vector of square matrix that has result vector resvec by Yacobi
    :param matrix: matrix
    :param resvec: result vector
    """
    print("Solution by Yacobi")
    if not isMatrix(matrix):
        return
    try:
        iter(resvec)
    except TypeError as e:
        return  # resvec is not vector
    matrix1,resvec1 = makeDominantMatrix(matrix,resvec)
    if len(matrix1) > 0:
        matrix = matrix1
        resvec = resvec1
    else:
        print("Matrix don't have dominant diagonal, trying iteration")
    #iteration
    iterationMat = buildIterationMatrix(matrix,resvec)
    if len(iterationMat) < 1:
        print("Can't build iteration matrix")
        return
    i = 0
    variables = []
    for j in matrix:
        variables.append(0) #start all variables to 0
    tempvariables = variables.copy()
    stop = False
    while i < numOfIterations:
        printStr = str(i) + "   "
        j = 0
        while j < len(variables):
            printStr += "X" + str(j) + " = " + "%.10f" % variables[j] + " "
            j += 1
        print(printStr)
        if stop:
            print("Number of iterations: " + str(i + 1))
            return
        j = 0
        while j < len(variables):
            result = iterationMat[j][0]
            k = 0
            while k < len(variables):
                result += iterationMat[j][k+1]*variables[k]
                k += 1
            tempvariables[j] = result
            j += 1
        #check if X0r+1 - X0r <=epsilon
        if abs(variables[0]-tempvariables[0]) <= epsilon:
            stop = True
        variables = tempvariables.copy()
        i += 1
    if not stop:
        print("The system does not converge")
    else:
        printStr = str(i) + "   "
        j = 0
        while j < len(variables):
            printStr += "X" + str(j) + " = " + "%.10f" % variables[j] + " "
            j += 1
        print(printStr)

def solutionOfMatrixByGaussZidel(matrix, resvec):
    """
    Prints the solution vector of square matrix that has result vector resvec by Gauss Zidel
    :param matrix: matrix
    :param resvec: result vector
    """
    print("Solution by Gauss Zidel")
    if not isMatrix(matrix):
        return
    try:
        iter(resvec)
    except TypeError as e:
        return  # resvec is not vector
    matrix1, resvec1 = makeDominantMatrix(matrix, resvec)
    if len(matrix1) > 0:
        matrix = matrix1
        resvec = resvec1
    else:
        print("Matrix don't have dominant diagonal, trying iteration")
    # iteration
    iterationMat = buildIterationMatrix(matrix, resvec)
    if len(iterationMat) < 1:
        print("Can't build iteration matrix")
        return
    i = 0
    variables = []
    for j in matrix:
        variables.append(0)  # start all variables to 0
    stop = False
    while i < numOfIterations:
        tempvariables = variables.copy()
        printStr = str(i) + "   "
        j = 0
        while j < len(variables):
            printStr += "X" + str(j) + " = " + "%.10f" % variables[j] + " "
            j += 1
        print(printStr)
        if stop:
            print("Number of iterations: " + str(i + 1))
            return
        j = 0
        while j < len(variables):
            result = iterationMat[j][0]
            k = 0
            while k < len(variables):
                result += iterationMat[j][k + 1] * variables[k]
                k += 1
            variables[j] = result
            j += 1
        # check if X0r+1 - X0r <=epsilon
        if abs(variables[0] - tempvariables[0]) <= epsilon:
            stop = True
        i += 1
    if not stop:
        print("The system does not converge")
    else:
        printStr = str(i) + "   "
        j = 0
        while j < len(variables):
            printStr += "X" + str(j) + " = " + "%.10f" % variables[j] + " "
            j += 1
        print(printStr)

choice = -1
while choice != '0':
    print("Enter 1 to Yacobi, 2 to Gauss Zidel or 0 to exit: ")
    choice = input()
    if choice == '1':
        solutionOfMatrixByYacobi(matrixA, vectorB)
    elif choice == '2':
        solutionOfMatrixByGaussZidel(matrixA, vectorB)
    elif choice != '0':
        print("Invalid input")