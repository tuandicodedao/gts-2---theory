#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#define N 50
#define ep 1e-10
#define step 0.25
using namespace std;

// hàm gán nội dung của ma trận b cho ma trận a
void assign (double a[][N], double b[][N], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = b[i][j];
        }
    }
}

// nhân 2 ma trận vuông a và b, kết quả lưu vào c
void matrixMultiply (double a[][N], double b[][N], double c[][N], int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            c[i][j] = 0;
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k]*b[k][j];
            }
            if (c[i][j] < ep && c[i][j] > -ep)
                c[i][j] = 0;
        }
}

// in ma trận cỡ m x n
void print (double matrix[][N], int m, int n){
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            cout << setw(12) << matrix[i][j];
        cout << endl;
    }
    cout << endl;
}

// đổi chỗ hàng, cột i và j
void swapRowColumn (double a[][N], int n, int i, int j) {
    // đổi chỗ cột i và j
    for (int k = 0; k < n; k++) {
        double temp = a[k][i];
        a[k][i] = a[k][j];
        a[k][j] = temp;
    }
    // đổi chỗ hàng i và j
    for (int k = 0; k < n; k++)
	{
		double temp = a[i][k];
		a[i][k] = a[j][k];
		a[j][k] = temp;
	}
}

// tính giá trị của đa thức tại x theo sơ đồ Horner
double horner (double a[], int n, double x) {
    double result = a[0];
    for (int i = 1; i < n; i++) {
        result = a[i] + x*result;
    }
    return result;
}

// hàm tìm khoảng chứa tất cả các nghiệm thực
vector<double> distance (double a[], int n) {
    // vecto chứa khoảng nghiệm
    vector <double> distance;
    // cận trên và dưới của khoảng chứa nghiệm thực
    double upper, lower;
    int i;
    // tìm hệ số âm đầu tiên của đa thức
    for (i = 1; i < n; i++) {
        if (a[i] < 0)
            break;
    }
    // các hệ số dương hết thì pt k có nghiệm dương
    if (i == n) {
        upper = 0;
    }
    else {
        int max = i;
        // tìm hệ số có trị tuyệt đối lớn nhất trong các hệ số âm
        for (int j = i + 1; j < n - 1; j++) {
            if (a[j] < 0 && -a[j] > -a[max]) {
                max = j;
            }
        }
        // cận trên của khoảng chứa nghiệm thực
        upper = 1 + pow(-a[max] / a[0], 1.0 / i);
    }
    
    double b[n];    // ma trận chứa hệ số của đa thức sau khi thay x thành -x
    // trường hợp n lẻ, tức bậc lớn nhất của đa thức là chẵn
    if (n % 2 == 1) {
        for (int i = 0; i < n - 1; i++) {
            // tại các hệ số của x có số mũ chẵn thì k đổi
            if (i % 2 == 0)
                b[i] = a[i];
            // tại các hệ số của x có số mũ lẻ thì đổi dấu
            else b[i] = -a[i];
        }
    }
    // trường hợp n chẵn, tức bận lớn nhất của đa thức là lẻ
    else {
        for (int i = 0; i < n - 1; i++) {
            if (i % 2 == 0) {
                b[i] = a[i];
            }
            else {
                b[i] = -a[i];
                b[n-1] = -a[n-1];
            }
        }
    }

    // tìm cận dưới của khoảng chứa nghiệm thực
    for (i = 1; i < n; i++) {
        if (b[i] < 0)
            break;
    }
    if (i == n)
        lower = 0;
    else {
        int max = i;
        // tìm hệ số có trị tuyệt đối lớn nhất trong các hệ số âm
        for (int j = 2; j < n - 1; j++) {
            if (b[j] < 0 && -b[j] > -b[i]) {
                max = j;
            }
        }
        lower = -1 - pow(-b[max] / b[0], 1.0/i);
    }
    distance.push_back(lower);
    distance.push_back(upper);
    return distance;
}

// hàm giải pt bằng chia đôi
double bisection(double arr[], int n, double a, double b) {
    double c;
    do {
        c = (a + b) / 2;
        if(horner(arr, n, c) == 0)
            return c;
        
        if (horner(arr, n, c) * horner(arr, n, a) < 0) {
            b = c;
        }
        else a = c;
    } while (abs(b-a) > ep);
	return c;
}

// hàm giải pt đa thức chính
// các nghiệm dc lưu vào vecto eigenvalue
void solve (double a[], int n, vector<double> &eigenvalue) {
    vector <double> ab = distance(a, n);
    double c = *ab.begin();
    //cout << c << endl;
    double b = *(ab.end() - 1);
    //cout << b << endl;
    // kiểm tra xem 2 cận trên và dưới có là nghiệm hay k
    if (horner(a, n, b) == 0)
        eigenvalue.push_back(b);
    // giải nghiệm
    for (double i = c; i < b; i = i + step) {
        if (horner(a, n, i) == 0) 
        eigenvalue.push_back(i);
        if (horner(a, n, i) * horner(a, n, i + step) < 0) {
            double x = bisection(a, n, i, i + step);
            eigenvalue.push_back(x);
            // cout << x << endl;
        }
    }
    return;
}

// hàm lấy hệ số của pt đa thức sau khi biến ma trận input về dạng Frobenius
// và lưu các hệ số vào ma trận b
void getcoefficient (double a[][N], int n, double b[]) {
    b[0] = 1;
    for (int i = 1; i <= n; i++) {
        b[i] = -a[0][i-1];
    }
}

// đổi chỗ hàng i và j
void swap (double matrix[][N], int n, int i, int j) {
    for (int k=0; k<=n; k++)
	{
		double temp = matrix[i][k];
		matrix[i][k] = matrix[j][k];
		matrix[j][k] = temp;
	}
}

//tính nghiệm và in nghiệm
void nghiem (double matrix[][N], int n){
    double x[n];
    //vị trí cột có giá trị khác 0 đầu tiên ở từng hàng
    int index[n];
    for (int i = 0; i<n; i++) {
        int j = 0;
        while (matrix[i][j] == 0 && j < n) j++;
        if (j == n) index[i] = -1;
        else index[i] = j;
    }

    //giải nghiệm
    for (int i = 0; i < n; i++) {
        if (index[i] != i)
            x[i] = 1;
    }
    for (int i = n-1; i >= 0; i--) {
        if (index[i] != -1) {
            x[index[i]] = 0;
            for (int j=index[i]+1; j < n; j++)
                x[index[i]] -= x[j]*matrix[i][j];
            x[index[i]] = x[index[i]]/matrix[i][index[i]];
        }
    } 
    cout << "Vecto rieng:" << endl;
    for (int i = 0; i < n; i++)
        cout << setw(12) << x[i] << endl;
}

void guass(double matrix[][N], int n) {
    //ưu tiên giá trị lớn nhất
    for (int i = 0; i < n-1; i++) {
        int key = i;
        for (int j = i; j < n; j++) {
            if (abs(matrix[j][i]) == 1) {
                key = j; 
                break;
            }
            if (abs(matrix[j][i])>abs(matrix[key][i]))
                key = j;
        }
        double khu = matrix[key][i];
        if (key != i) swap (matrix, n, key, i);
    //lưu vị trí cột có giá trị khác 0 đầu tiên của từng hàng
    int index[n];
    for (int i = 0; i<n; i++) {
        int j = 0;
        while (matrix[i][j] == 0  && j < n) j++;
        if (j == n) index[i] = -1;
        else index[i] = j;
    } 
    //bậc thang
    //cout << "\ni = " << i << endl;
    for (int j = i+1; j<n; j++) {
       // cout << index[j] << " ";
        if (index[j] == index[i])  {
            double h = matrix[j][index[i]]/matrix[i][index[i]];
        for (int k = i+1; k<=n; k++) {
            matrix[j][k] -= matrix[i][k]*h;
            if ((matrix[j][k] < 5*ep && matrix[j][k] >-5*ep) 
                || (matrix[j][k] < 50*ep && matrix[j][k] > -50*ep))
                matrix[j][k] = 0;
        }
        matrix[j][i] = 0;
        }
    }
    //print(matrix, n, n+1);
    }
    nghiem(matrix, n);
}

// biến đổi ma trận về dạng Frobenius
void danilevskii (double a[][N], int n, double p[][N], vector<double> &eigenvalue) {
    double c[n][N];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                c[i][j] = 1;
            else c[i][j] = 0;
        }
    }
    // số lần lặp
    int i = 0;
    // duyệt các phần tử xét từ dưới lên, từ phải sang
    for (int k = n-2; k >= 0; k--) {
        // trường hợp phần tử xét ở hàng cuối cùng = 0
        if (a[n-1][n-2] == 0) {
            int x = k;
            for (int h = n-3; h >= 0; h--) {
                if (a[n-1][h] != 0) {
                    x = h;
                    break;
                }
            }
            // tìm dc phần tử khác 0 đứng trc phần tử đang xét
            if (x != k) {
                cout << "Doi cho hang va cot " << x+1 << " " << k+1 << endl;
                swapRowColumn(a, n, x, k);
                // đổi chỗ cột x và k của ma trận c
                for (int r = 0; r < n; r++) { // r stands for row
                    double temp = c[r][x];
                    c[r][x] = c[r][k];
                    c[r][k] = temp;
                }
                print(a, n, n);
            }
            // các phần tử đứng trc phần tử đang xét đều = 0
            else {
                cout << "Co 1 tri rieng la: " << a[n-1][n-1] << endl;
                eigenvalue.push_back(a[n-1][n-1]);
                cout << "Xet ma tran con:" << endl;
                n--;
                print(a, n, n);
                // if (n == 1) {
                //     cout << "Co 1 tri rieng la: " << a[0][0] << endl;
                //     eigenvalue.push_back(a[0][0]);

                // }
                continue;
            }
        }
        // trường hợp phần tử xét không ở hàng cuối = 0
        if (a[k+1][k] == 0 && k != n-2) {
            int x = k;
            for (int h = k-1; h >= 0; h--) {
                if (a[k+1][h] != 0) {
                    x = h;
                    break;
                }
            }
            // tìm dc phần tử khác 0 đứng trc phần tử đang xét
            if (x != k) {
                cout << "Doi cho hang va cot " << x+1 << " " << k+1 << endl;
                swapRowColumn(a, n, x, k);
                // đổi chỗ cột x và k của ma trận c
                for (int r = 0; r < n; r++) { // r stands for row
                    double temp = c[r][x];
                    c[r][x] = c[r][k];
                    c[r][k] = temp;
                }
                print(a, n, n);
            }
            // các phần tử đứng trc phần tử đang xét đều = 0
            else {
                // lấy hệ số của pt đa thức
                double b[n - k];
                b[0] = 1;
                for (int i = 1; i < n-k; i++) {
                    b[i] = -a[k+1][k+i];
                }
                // in hệ số
                // for (int i = 0; i < n-k; i++) {
                //     cout << b[i] << " ";
                // }
                // cout << endl;
                
                // giải pt đa thức bên phải
                solve(b, n-k, eigenvalue);

                n = k + 1;
                cout << "Xet ma tran con:" << endl;
                print (a, n, n);
                if (n == 1) {
                    cout << "Co 1 tri rieng la: " << a[0][0] << endl;
                    eigenvalue.push_back(a[0][0]);
                }
                continue;
            }
        }
        // tìm m1 và m
        double m1[n][N];
        double m[n][N];
        double b[n][N];
       
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != k) {
                    if (i == j) {
                        m[i][j] = 1;
                        m1[i][j] = 1;
                    }
                    else {
                        m[i][j] = 0;
                        m1[i][j] = 0;
                    }
                }
                else {
                    m1[i][j] = a[k+1][j];
                    if (j == k) m[i][j] = 1.0/a[k+1][k];
                    else m[i][j] = -a[k+1][j]/a[k+1][k];
                }
            }
        }

        // cout << "Ma tran M1:" << endl;
        // print(m1, n, n);
        // cout << "Ma tran M:" << endl;
        // print(m, n, n);
        
        matrixMultiply(a, m, b, n);
        matrixMultiply(m1, b, a, n);
        i++;
        
        matrixMultiply(c, m, p, n);
        assign(c, p, n);

        cout << "A_" << i << ":" << endl;
        print(a, n, n);
        // cout << "Ma tran P:" << endl;
        // print(p, n, n);
    }
    // đã về dạng Frobenius, ta tìm trị riêng bằng giải pt đa thức
    double b[n+1];
    getcoefficient(a, n, b);
    solve(b, n+1, eigenvalue);

    cout << "Cac tri rieng: ";
    int size = eigenvalue.size();
    for (int i = 0; i < size; i++) {
        cout << setw(12) << eigenvalue[i];
    }

    return;
}

// hàm tìm vecto riêng sau khi có các giá trị riêng
void eigenvector(double A[][N], double a[][N], int n, double p[][N], vector<double> eigenvalue) {
    int size = eigenvalue.size();

    cout << setprecision(12);
    
    for (int i = 0; i < size; i++) {
        cout << "\nVoi tri rieng " << eigenvalue[i] << " giai he phuong trinh theo Gauss:" << endl;
        double b[n][N];
        for (int i = 0; i < n; i++) {
            b[i][n] = 0;
        }
        for (int h = 0; h < n; h++) {
            for (int k = 0; k < n; k++) {
                if (h == k) {
                    b[h][k] = A[h][k] - eigenvalue[i];
                }
                else b[h][k] = A[h][k];
            }
        }
        //print(b, n, n+1);
        guass(b, n);
    }
}

int main() {
    // đọc ma trận từ file Nhap.txt
    ifstream read ("Nhap.txt");
    int n = 0;
    string s;
    while (getline(read, s))
        n++;
    cout << "Kich co ma tran: " << n << endl;
    read.close();

    // ma trận Input
    double a[n][N];
    read.open ("Nhap.txt");
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            read >> a[i][j];
    read.close();
    cout << "Ma tran input:" << endl;
    print(a, n, n);
    double A[n][N];
    assign(A, a, n);

    // vecto chứa các giá trị riêng
    vector<double>eigenvalue;
    // ma trận P
    double p[n][N];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                p[i][j] = 1;
            else p[i][j] = 0;
        }
    }

    danilevskii(a, n, p, eigenvalue);
    eigenvector(A, a, n, p, eigenvalue);
    return 0;
}