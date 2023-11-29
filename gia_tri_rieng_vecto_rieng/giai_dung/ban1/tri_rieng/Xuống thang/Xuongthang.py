# Khai báo thư viện
import math
import numpy as np
import numpy.linalg as la

# Sai số
eps = 1e-6

# Nhận ma trận
A = np.loadtxt('Xuống thang\XuongthangInput.txt')

# Ma trận đơn vị
I = np.eye(A.shape[1], dtype=float)

# Chuẩn hóa
def normalize(M):
    maxi = abs(M[0][0])
    c = M[0][0]
    for i in range(len(M)):
        if maxi < abs(M[i][0]):
            maxi = abs(M[i][0])
            c = M[i][0]
    matrix = M / c
    return matrix

# Tính A^m*x và lưu vào cột m của B


def calc(A, m, B):
    M = np.dot(A, B[m-1])
    M = normalize(M)
    B.append(M)
    return

# Kiểm tra


def check(A, x):
    B = []
    B.append(x)
    calc(A, 1, B)
    calc(A, 2, B)
    m = 2
    TH = 3
    while True:
        if m == 100:
            break
        if (la.norm(B[m]-B[m-1], np.inf) < eps):
            TH = 1
            break
        if (la.norm(B[m]-B[m-2], np.inf) < eps):
            TH = 2
            break
        m += 1
        calc(A, m, B)
    return m, TH, B

# Xuống thang đơn, thực


def deflate1(A, v):
    for i in range(len(v)):
        if (v[i] == 1):
            index = i
            break
    theta = I
    for i in range(theta.shape[1]):
        theta[i, index] = theta[i, index] - v[i]
    A = np.dot(theta, A)
    return A, theta

# Xuống thang trái dấu, thực


def deflate2(A, v1, v2):
    A, theta1 = deflate1(A, v1)
    v2 = normalize(np.dot(theta1, v2))
    A, theta2 = deflate1(A, v2)
    return A, theta2

# Xử lý TH1


def TH_1(x, i, A, B):
    v = B[-1]
    Am1 = np.dot(A, B[-1])
    for j in range(A.shape[1]):
        if (Am1[j, 0] != 0):
            lamda = Am1[j, 0]/B[-1][j, 0]
            break
    for j in range(len(B)):
        print("\nA{}Y = ".format(j))
        print(B[j])
    print('\nGia tri rieng thu {} la: '.format(i+1), lamda)
    print('\nVector rieng tuong ung la v[{}] la '.format(i+1))
    print(v)
    print("\n")
    print('*' * 100)
    A,theta = deflate1(A,v)
    print ("\nMa tran theta")
    print(theta)
    print("\nMa tran A moi")
    print(A)
    print("\n")
    print('*' * 100)
    i += 1
    return A,i

# Xử lý TH2
def TH_2(x,i,A,m,B):
    Am1 = np.dot(A,B[-1])
    Am2 = np.dot(A,Am1)
    Am3 = np.dot(A,Am2)
    for j in range(A.shape[1]):
        if (Am1[j,0] != 0):
            lamda1 = math.sqrt(Am3[j,0]/Am1[j,0])
            lamda2 = -math.sqrt(Am3[j,0]/Am1[j,0])
            break
    if (m%2 == 0):
        v1 = normalize(Am2 + lamda1*Am1)
        v2 = normalize(Am2 + lamda2*Am1)
    else:
        v1 = normalize(Am3 + lamda1*Am2)
        v2 = normalize(Am3 + lamda2*Am2)

    print('\nGia tri rieng thu {} la: '.format(i+1), lamda1)
    print('\nVector rieng tuong ung la v[{}] la '.format(i+1))
    print(v1)

    print('\nGia tri rieng thu {} la: '.format(i+2), lamda2)
    print('\nVector rieng tuong ung la v[{}] la '.format(i+2))
    print(v2) 
    print("\n")
    print('*' * 100)
    A, theta2 = deflate2(A,v1,v2)
    print (" \n Ma tran theta")
    print(theta2)
    print("\nMa tran A moi")
    print(A)
    print("\n")
    print('*' * 100)
    i += 2
    return A,i

# Xử lý TH3
def TH_3(x,i,A,B):
    Am = B[-1]
    Am1 = np.dot(A,Am)
    Am2 = np.dot(A,Am1)
    a1 = Am2[0][0]
    b1 = Am1[0][0]
    c1 = Am[0][0]
        
    a2 = Am2[1][0]
    b2 = Am1[1][0]
    c2 = Am[1][0]

    a = c1 * b2 - b1 * c2
    b = (a1 * c2 - c1 * a2) 
    c = (b1 * a2 - a1 * b2)
    delta = b ** 2 - 4 * a * c
    print("\nGia pt tim lam da la: a*Z^2 + b*Z + c = 0")
    print("\na = ", a)
    print("\nb = ", b)
    print("\nc = ", c)
    print("\ndelta = ", delta)

    if (delta < 0):
        lamda1 = complex(-b / (2 * a), -math.sqrt(abs(delta)) / (2 * a))
        lamda2 = complex(-b / (2 * a), math.sqrt(abs(delta)) / (2 * a))
        v1 = normalize(Am2 - lamda2 * Am1)
        v2 = normalize(Am2 - lamda1 * Am1)
        
        print('\nGia tri rieng thu {} la: '.format(i+1), lamda1)
        print('\nVector rieng tuong ung la v[{}] la '.format(i+1))
        print(v1)        

        print('\nGia tri rieng thu {} la: '.format(i+2), lamda2)
        print('\nVector rieng tuong ung la v[{}] la '.format(i+2))
        print(v2) 
        print("\n")
        print('*' * 100)        
    else: print('delta khong am')

# Chương trình chính
def power_iteration(A):
    print(A)
    x = np.random.rand(A.shape[1],1)
    i=0
    while (i<A.shape[1]):
        m,TH,B = check(A,x)
        print('Sau {} lần lặp'.format(m))
        print('Đây là TH ',TH)
        if (TH == 1): 
            A,i = TH_1(x,i,A,B)
        if (TH == 2): 
            A,i = TH_2(x,i,A,m,B)
        if (TH == 3): 
            TH_3(x,i,A,B)
            break

power_iteration(A)