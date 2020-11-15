#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;

double *y = nullptr, * b = nullptr, * res = nullptr;
double **a = nullptr, * *l = nullptr;
void solve(int n);
void lu(int n);

int main() {
	int n = 0;
	cout << "输入方程阶数：";
	cin >> n;
	a = new double* [n + 5];
	l = new double* [n + 5];
	b = new double[n + 5] { 0 };
	y = new double[n + 5] { 0 };
	res = new double[n + 5] { 0 };
	for (int i = 0; i <= n; ++i) {
		a[i] = new double[n + 5] { 0 };
		l[i] = new double[n + 5] { 0 };
	}
	cout << "输入矩阵：" << endl;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			cin >> a[i][j];
			if (i == 0)
				l[i][j] = a[i][j];
		}
	cout << "输入常向量：" << endl;
	for (int i = 0; i < n; ++i)
		cin >> b[i];
	for (int i = 1; i < n; ++i)
		l[i][0] = a[i][0] / l[0][0];
	lu(n);
	cout << "LU矩阵: " << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j)
			cout << l[i][j] << " ";
		cout << endl;
	}
	solve(n);
	cout << "Y向量: ";
	for (int i = 0; i < n; ++i)
		cout << y[i] << " ";
	cout << endl;
	cout << "Res向量: ";
	for (int i = 0; i < n; ++i)
		cout << res[i] << " ";
	cout << endl;
	for (int i = 0; i <= n; ++i) {
		delete[] a[i];
		delete[] l[i];
	}
	delete[] a;
	delete[] l;
	delete[] b;
	delete[] y;
	delete[] res;
	return 0;
}

void solve(int n) {
	y[0] = b[0];
	for (int i = 1; i < n; ++i) {
		double mid = 0;
		for (int j = 0; j < i; ++j)
			mid += (y[j] * l[i][j]);
		y[i] = b[i] - mid;
	}
	res[n - 1] = (y[n - 1] / l[n - 1][n - 1]);
	for (int i = n - 2; i >= 0; --i) {
		double mid = 0;
		for (int j = i; j < n; ++j)
			mid += (l[i][j] * res[j]);
		res[i] = (((y[i] - mid)) / l[i][i]);
	}
}

void lu(int n) {
	for (int i = 1; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			double mid = 0;
			for (int k = 0; k < i; ++k)
				mid += (l[i][k] * l[k][j]);
			l[i][j] = a[i][j] - mid;
		}
		for (int j = i + 1; j < n; ++j) {
			double mid = 0;
			for (int k = 0; k < i; ++k)
				mid += (l[j][k] * l[k][i]);
			l[j][i] = (a[j][i] - mid) / l[i][i];
		}
	}
}