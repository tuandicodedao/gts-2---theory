#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstdlib>


using namespace std;

	double A[20][20], Av[20][50], v[20][50], eps = 1e-06; 
	int m, n;

// Dữ liệu đầu vào
void input()
{
    cin >> m;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            cin >> A[i][j];
	cout << "Ma tran A" << endl;
	for(int i = 0; i < m; i++){
		for(int j = 0; j < m; j++){
			cout << setw(15) << right << setprecision(10) << fixed << A[i][j];
			Av[i][j] = 0;}
		cout << endl;}
}

//Nhân ma trận A với 1 vecto
void Nhan(double A[][20], double v[][50], double Av[][50], int m, int k)
{
	for(int i = 0; i < m; i++)
		for(int j = 0; j < m; j++)
			Av[i][k] += A[i][j]*v[j][k-1];
}

//Tìm max của vecto Av
double max(double Av[][50], int m, int k)
{
	double max = Av[0][k];
	for(int i = 1; i < m; i++)
		if(Av[i][k] > max) {
			max = Av[i][k];
		}
	if (max == 0) return 1;
	return max;
}

//Điều kiện dừng cho trường hợp 1
double dkDung1(double v[][50], int m, int k)
{
	double max = 0, tmp[m];
	for(int i = 0; i < m; i++){
		tmp[i] = fabs(v[i][k-1] - v[i][k-2]);
		if(tmp[i] > max) max = tmp[i];
	}
	return max;
}

//Điều kiện dừng cho trường hợp 2
double dkDung2(double v[][50], int m, int k)
{
	double max = 0, tmp;
	for(int i = 0; i < m; i++){
		tmp = fabs(v[i][k] - v[i][k-2]);
		if(tmp > max) max = tmp;
	}
	return max;
}

//Tính giá trị riêng trường hợp 1
double lamda1(double Av[][50], double v[][50], int m, int k)
{
	double sum1 = 0, sum2 = 0;
	for(int i = 0; i < m; i++){
		sum1 += Av[i][k]*v[i][k-1];
		sum2 += v[i][k-1]*v[i][k-1];}
	return sum1/sum2;
}

//Tính giá trị riêng trường hợp 3
void lamda3(double Av[][50], double real[], double img[], int m, int k, int count)
{
	double a, b, c, delta;
	int i, j;
	do{
			srand((unsigned int) time(NULL));
			i = rand() % m;
			j = rand() % m;
	}while((i == j) || (Av[i][k] == 0 || Av[j][k] == 0));
	a = Av[i][k-1]*Av[j][k-2] - Av[j][k-1]*Av[i][k-2];
	b = Av[i][k-2]*Av[j][k] - Av[j][k-2]*Av[i][k];
	c = Av[i][k]*Av[j][k-1] - Av[j][k]*Av[i][k-1];
	delta = b*b - 4*a*c;
	real[count] = -b/(2*a);
	img[count] = sqrt(fabs(delta))/(2*a);
}

void Kiemtra(double A[][20], double lamda, double x[], int m)
{
	double I[m][m], c[m][m], b[m];
	for(int i = 0; i < m; i++)
		for(int j = 0; j < m; j++)
			if(i != j) I[i][j] = 0;
			else I[i][j] = lamda;
	for(int i = 0; i < m; i++)
		for(int j = 0; j < m; j++)
			c[i][j] = A[i][j] - I[i][j];
	for(int i = 0; i < m; i++){
		double sum = 0;
		for(int j = 0; j < m; j++)
			sum += c[i][j]*x[j];
		b[i] = sum;
	}
	cout << "______________________________" << endl;
	cout << "Ma tran kiem tra xap xi ma tran 0 => true" << endl;
	for(int i = 0; i < m; i++)
		cout << setw(15) << right << setprecision(10) << fixed << b[i] << endl;
	cout << "______________________________" << endl;
}

int main()
{
//Nhap file va in ra ma tran A
	freopen("input.txt", "r", stdin);
    freopen("LuyThuaOutput.txt", "w", stdout);
    input();
	
//Truong hop 1: |lamda1| > |lamda2|
{
	double lamda[50] = {0};
	cout << "Xet truong hop 1 tri rieng troi" << endl;
	int k = 0;
	for(int i = 0; i < m; i++)
		v[i][k] = 1;
	k++;
	Nhan(A, v, Av, m, k);
	double p = max(Av, m, k);
	for(int i = 0; i < m; i++)
		v[i][k] = Av[i][k]/p;
	k++;
	while(dkDung1(v, m, k) > eps) {
		Nhan(A, v, Av, m, k);
		p = max(Av, m, k);
		for(int i = 0; i < m; i++)
			v[i][k] = Av[i][k]/p;
		lamda[k] = lamda1(Av, v, m, k);
		k++;
		if (k > 50) break;
	}
	if(abs(lamda[k-1] - lamda[k-2]) > 1e-03) cout << "Khong thoa man" << endl;
	else {
		for(int i = 0; i < k; i++){
			cout << "==========================================" << endl;
			cout << "k = " << i << endl;
			cout << "Ma tran Av" << endl;
			for(int j = 0; j < m; j++)
				cout << setw(15) << right << setprecision(10) << fixed << Av[j][i] << endl;
			cout << "---------------------" << endl;
			cout << "lamda[" << i << "] = " << lamda[i] << endl;
		}
		cout << "==========================================" << endl;
		cout << "Ket thuc sau " << k+1 << " lan lap" << endl;
		cout << "---------------------" << endl;
		cout << "lamda = " << setprecision(10) << fixed << lamda[k-1] << endl;
		cout << "---------------------" << endl;
		cout << "Vecto rieng x tuong ung" << endl << endl;
		double x[m];
		for(int i = 0; i < m; i++){
			x[i] = v[i][k-1];
			cout << setprecision(10) << fixed << v[i][k-1] << endl;
		}
		Kiemtra(A, lamda[k-1], x, m);
		return 0;}
}

//Truong hop 2: |lamda1| = |lamda2| > |lamda3|, lamda1 = - lamda2
{
	cout << "====================================" << endl;
	cout << "Xet truong hop 2 tri rieng troi thuc" << endl;
	for(int i = 0; i < m; i++)
		for(int j = 0; j < 50; j++)
			Av[i][j] = 0;
	int k = 0;
	for(int i = 0; i < m; i++)
		v[i][k] = 1;
	k++;
	Nhan(A, v, Av, m, k);
	double p = max(Av, m, k);
	for(int i = 0; i < m; i++)
		v[i][k] = Av[i][k]/p;
	k++;
	double lamdabinhphuong[50] = {0};
	do {
		Nhan(A, Av, Av, m, k);
		p = max(Av, m, k);
		int index;
		for(int i = 0; i < m; i++)
			if(Av[i][k] == p){
				index = i;
				break;
			}
		for(int i = 0; i < m; i++)
			v[i][k] = Av[i][k]/p;
		if(k > 2) lamdabinhphuong[k] = Av[index][k]/Av[index][k-2];
		k++;
		if(k > 50) break;
	} while(dkDung2(v, m, k) > eps);
	if(fabs(lamdabinhphuong[k-1] - lamdabinhphuong[k-2]) > 1e-03)
		cout << "Khong thoa man" << endl;
	else {
		for(int i = 0; i < k; i++){
			cout << "==========================================" << endl;
			cout << "k = " << i << endl;
			cout << "Ma tran Av" << endl;
			for(int j = 0; j < m; j++)
				cout << setw(15) << right << setprecision(10) << fixed << Av[j][i] << endl;
			cout << "---------------------" << endl;
			cout << "lamdabinhphuong[" << i << "] = " << lamdabinhphuong[i] << endl;
		}
		cout << "==========================================" << endl;
		cout << "Ket thuc sau " << k << " lan lap" << endl;
		double tmp = sqrt(lamdabinhphuong[k-1]);
		double tmp1 = sqrt(lamdabinhphuong[k-2]);
		int index1, index2;
		for(int i = 0; i < m; i++) {
			if(v[i][k-1] == 1) index1 = i;
			if(v[i][k-2] == 1) index2 = i;
		}
		cout << "---------------------" << endl;
		cout << "lamda1 = " << setprecision(10) << fixed << tmp << endl;
		double x1[m], x2[m];
		for(int i = 0; i < m; i++){
			Av[i][k-1] = Av[i][k-1]/pow(tmp, k-1);
			Av[i][k-2] = Av[i][k-2]/pow(tmp1, k-2);
			x1[i] = (Av[i][k-1] + Av[i][k-2])/2;
			x2[i] = (Av[i][k-1] - Av[i][k-2])/2;
		}
		double p1 = x1[index1], p2 = x2[index2];
		for(int i = 0; i < m; i++){
			x1[i] = x1[i]/p1;
			x2[i] = x2[i]/p2;
		}
		cout << "---------------------" << endl;
		cout << endl << "Vecto rieng x1 tuong ung" << endl;
		for(int i = 0; i < m; i++)
			cout << setprecision(10) << fixed << x1[i] << endl;
		Kiemtra(A, tmp, x1, m);
		tmp = -tmp;
		cout << "---------------------" << endl;
		cout << "lamda2 = " << setprecision(10) << fixed << tmp << endl;
		cout << "---------------------" << endl;
		cout << endl << "Vecto rieng x2 tuong ung" << endl;
		for(int i = 0; i < m; i++)
			cout << setprecision(10) << fixed << x2[i] << endl;
		Kiemtra(A, tmp, x2, m);
		return 0;
	}
}

//Truong hop 3: |lamda1| = |lamda2| > |lamda3|, lamda1, lamda2 phuc lien hop
{
	cout << "====================================" << endl;
	cout << "Xet truong hop 2 tri rieng troi phuc" << endl;
	for(int i = 0; i < m; i++)
		for(int j = 0; j < 50; j++)
			Av[i][j] = 0;
	int k = 0;
	double real[50], img[50];
	for(int i = 0; i < m; i++)
		v[i][k] = 1;
	cout << "==========================================" << endl;
	cout << "k = " << k << endl;
	cout << "Ma tran Av" << endl;
	for(int j = 0; j < m; j++)
		cout << setw(15) << right << setprecision(10) << fixed << Av[j][k] << endl;
	k++;
	Nhan(A, v, Av, m, k);
	cout << "==========================================" << endl;
	cout << "k = " << k << endl;
	cout << "Ma tran Av" << endl;
	for(int j = 0; j < m; j++)
		cout << setw(15) << right << setprecision(10) << fixed << Av[j][k] << endl;
	double p = max(Av, m, k);
	for(int i = 0; i < m; i++)
		v[i][k] = Av[i][k]/p;
	int count = 1;
	for(int i = 0; i < 50; i++)
		real[i] = img[i] = i;
	do {
		k++;
		if(k > 50) break;
		Nhan(A, v, Av, m, k);
		cout << "==========================================" << endl;
		cout << "k = " << k << endl;
		cout << "Ma tran Av" << endl;
		for(int j = 0; j < m; j++)
			cout << setw(15) << right << setprecision(10) << fixed << Av[j][k] << endl;
		p = max(Av, m, k);
		if(p == 1) continue;
		k++;
		if(k > 50) break;
		Nhan(A, Av, Av, m, k);
		cout << "==========================================" << endl;
		cout << "k = " << k << endl;
		cout << "Ma tran Av" << endl;
		for(int j = 0; j < m; j++)
			cout << setw(15) << right << setprecision(10) << fixed << Av[j][k] << endl;
		k++;
		if(k > 50) break;
		Nhan(A, Av, Av, m, k);
		cout << "==========================================" << endl;
		cout << "k = " << k << endl;
		cout << "Ma tran Av" << endl;
		for(int j = 0; j < m; j++)
			cout << setw(15) << right << setprecision(10) << fixed << Av[j][k] << endl;
		p = max(Av, m, k);
		for(int i = 0; i < m; i++)
			v[i][k] = Av[i][k]/p;
		lamda3(Av, real, img, m, k, count);
		cout << "--------------------" << endl;
		cout << "real[" << count << "] = " << real[count] << endl;
		cout << "img[" << count << "] = " << img[count] << endl;
		cout << "--------------------" << endl;
		count++;
	} while(fabs(img[count-1] - img[count-2]) > eps);
	if(fabs(img[count-1] - img[count-2]) > 1e-03) cout << "Khong thoa man" << endl;
	else {
		cout << "==========================================" << endl;
		cout << "Ket thuc sau " << k+1 << " lan lap" << endl;
		cout << "--------------------" << endl;
		cout << "lamda1 = " << setprecision(10) << fixed << real[count-1] << " + " << setprecision(10) << fixed << abs(img[count-1]) << "i" << endl;
		cout << "lamda2 = " << setprecision(10) << fixed << real[count-1] << " - " << setprecision(10) << fixed << abs(img[count-1]) << "i" << endl;
	}
}
 	return 0;
}
