from giaiDT import *        # Gói dùng để giải đa thức đặc trưng
import numpy as np

def simpleCaseA(A, row):
    '''
       Trường hợp lí tưởng : Phần tử A[row][row-1] != 0
       :argument
        A: Ma trận 2 chiều A
        row: chỉ số dòng đang xét
       :return
       similarA: Ma trận A sau khi đã biến đổi
       inverseM: ma trận nghịch đảo của M
    '''

    # tạo ra ma trận M cho phép biến đổi trên A
    n= A.shape[0]
    M = np.eye(n)
    M[row - 1, :] = A[row, :]

    # tạo ra ma trận nghịch đảo M
    inverseM = np.eye(n)
    inverseM[row - 1, :] = -A[row, :] / A[row, row - 1]
    inverseM[row - 1, row - 1] = 1 / A[row, row - 1]

    # similarA = M*A*inverseM
    similarA = M.dot(A).dot(inverseM)
    return similarA,inverseM

def specialCase1(A, row, col):
    '''
       Trường hợp đặc biệt 1: Nếu A[row][row-1] == 0 và A[row][col] != 0 với col < row -1
       :argument 
       Ma trận A ,
       row: chỉ số dòng đang xét
       col: chỉ số cột đang xét
       :return 
       copA: Ma trận mới sau khi biến đổi từ A
    '''
    copA= A.copy()
    copA[[row - 1,col], :] = copA[[col, row - 1], :]
    copA[:, [row - 1, col]] = copA[:, [col, row - 1]]
    return copA

def getPolynomialCoefficients(matrixFrobenius):
    '''
       Lấy đa thức đặc trưng từ ma trận Frobenius
       :argument
       Ma trận arrayFrobenius
       :return
       polynomialCoefficients: Ma trận 1 chiều chứa hệ số của đa thức đặc trưng
    '''

    size = matrixFrobenius.shape[0]
    polynomialCoefficients = (-1) ** size * np.ones(size + 1)
    polynomialCoefficients[1:] = polynomialCoefficients[1:] * (-1) * matrixFrobenius[0, :]
    return polynomialCoefficients

def findEigenValue(matrixFrobenius):
    """
    Tìm trị riêng từ ma trận Frobenius
    :argument
    arrayFrobenius: Ma trận Frobenius
    :return
    eigenValues: Danh sách các trị riêng đã được sắp xếp
    """

    p = getPolynomialCoefficients(matrixFrobenius)
    p = list(p)[::-1]
    eigenValues = findRoots(p)
    return sorted(eigenValues)

def Danilevski(A):
    '''
       :argument
       A : Ma trận 2 chiều
       :return
       list_eigenvalues : Danh sách các trị riêng
       list_eigenvectors : Danh sách các vector riêng của từng trị riêng
    '''
    sizeOfA = A.shape[0]
    back = np.eye(sizeOfA)
    similar = A.copy()
    list_eigenvalues = []
    list_eigenvectors = []
    row = sizeOfA - 1
    tmpSize = sizeOfA
    while row > 0:
        if similar[row, row - 1] != 0:                           # Trường hợp lí tưởng
            similar, inverseM = simpleCaseA(similar, row)
            back = back.dot(inverseM)
            print('Case lí tưởng, vị trí đang xét tại dòng thứ: ',row)
            print(similar)

        else:
            case1 = False
            for col in range(0, row - 1):
                if similar[row, col] != 0:                       # Trường hợp đặc biệt thứ 1:
                    similar = specialCase1(similar, row, col)
                    back[:, [col, row - 1]] = back[:, [row - 1, col]]
                    case1 = True
                    print('Case 1,vị trí đang xét tại dòng thứ: ', row)
                    print(similar)
                    row = row + 1
                    break

            if not case1:
                # Trường hợp đặc biệt thứ 2
                for col in range(row, tmpSize - 1):
                    S = np.eye(tmpSize)
                    S[:row, col + 1] = -similar[:row, col]
                    inverseS = np.eye(tmpSize)
                    inverseS[:row, col + 1] = similar[:row, col]
                    similar = S.dot(similar).dot(inverseS)
                    back = back.dot(inverseS)
                    print('Case 2, vị trí đang xét tại dòng thứ: ', row)
                    print(similar)

                case2_1 = False
                # Trường hợp đặc biệt thứ 2.1
                for _row in range(row - 1, -1, -1):
                    if similar[_row, tmpSize - 1] != 0:
                        U = np.eye(tmpSize)
                        x = U[row:tmpSize, :]
                        y = U[row - 1, :]
                        U = np.vstack([U[0:row - 1, :], x, y])
                        inverseU = U.T
                        similar = U.dot(similar).dot(inverseU)
                        back = back.dot(inverseU)
                        print('Case 2.1, vị trí đang xét tại dòng thứ: ', row)
                        print(similar)
                        row = tmpSize
                        case2_1 = True
                        break

                if not case2_1:
                    # Trường hợp đặc biệt thứ 2.2
                    print('Case 2.2, vị trí đang xét tại dòng thứ: ', row)
                    print(similar)
                    X = similar[row:, row:]      # X đã là ma trận Frobenius
                    sizeOfX = X.shape[0]
                    eigenValues = findEigenValue(X)
                    for j in range(len(eigenValues)):
                        print("Vector riêng ứng với trị riêng đây(case2.2): ", eigenValues[j])
                        list_eigenvalues.append(eigenValues[j])
                        eigenVectorOfX = np.power(eigenValues[j], np.arange(sizeOfX))[::-1].reshape((sizeOfX, 1))
                        eigenVectorOfA = np.zeros((sizeOfA, 1))
                        p = np.zeros((tmpSize, 1))
                        p[row:tmpSize] = eigenVectorOfX
                        p = back.dot(p)
                        eigenVectorOfA[:tmpSize] = p
                        list_eigenvectors.append(eigenVectorOfA)
                        print(eigenVectorOfA)
                    tmpSize = row
                    similar = similar[:row, :row]
                    back = np.eye(tmpSize)
        row = row - 1
    X = similar        # X đã là ma trận Frobenius
    sizeOfX = X.shape[0]
    print('Dòng thứ: ', row)
    print(similar)
    eigenValues = findEigenValue(X)
    for j in range(len(eigenValues)):
        print("Vector riêng ứng với trị riêng: ", eigenValues[j])
        list_eigenvalues.append(eigenValues[j])
        eigenVectorOfX = np.power(eigenValues[j], np.arange(sizeOfX))[::-1].reshape((sizeOfX, 1))
        p = np.zeros((tmpSize, 1))
        eigenVectorOfA = np.zeros((sizeOfA, 1))
        p[row:tmpSize] = eigenVectorOfX
        p = back.dot(p)
        eigenVectorOfA[:tmpSize] = p
        list_eigenvectors.append(eigenVectorOfA)
        print(eigenVectorOfA)
    return list_eigenvalues, list_eigenvectors
def readMatrix(path):
    file = open("A.txt", 'r')
    lt = list([int(number) for number in line.split()] for line in file)
    A = np.array(lt)
    return A

if __name__ == '__main__':
    A = readMatrix('E://Test1.txt')
    if A.shape[0] != A.shape[1]:
        print('Ma trận A không phải là ma trận vuông!')
    else:
        eigenvalues, eigenvectors = Danilevski(A)
        print('Kiểm tra kết quả: ')
        for i in range(len(eigenvalues)):
            B = A.dot(eigenvectors[i])
            C = eigenvalues[i] * eigenvectors[i]
            print('landa = ', eigenvalues[i])
            print('A.X - landa.X = 0:\n', (B - C))