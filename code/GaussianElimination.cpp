#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

double **a = nullptr, * res = nullptr;
void Select_c(int k, int n);
bool Guass(int n);

int main() {
	int n = 0;
	cout << "输入方程阶数：";
	cin >> n;
	a = new double* [n + 5];
	res = new double[n + 5] { 0 };
	for (int i = 0; i <= n; ++i)
		a[i] = new double[n + 5] { 0 };
	cout << "输入矩阵：" << endl;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j <= n; ++j)
			cin >> a[i][j];
	if (Guass(n)) {
		for (int i = 0; i < n; ++i)
			cout << "x" << i + 1 << " = " << res[i] << endl;
		cout << endl;
	} else
		cout << "无解！" << endl;
	for (int i = 0; i <= n; ++i)
		delete[] a[i];
	delete[] a;
	delete[] res;
	return 0;
}

bool Guass(int n) {
	for (int i = 0; i < n - 1; ++i) {
		Select_c(i, n);
		for (int j = i + 1; j < n; ++j) {
			double mid = (a[i][i] / a[j][i]);
			for (int k = i; k <= n; ++k) {
				if (k == i)
					a[j][k] -= ((a[i][k] / a[i][i]) * a[j][i]);
				else
					a[j][k] -= (a[i][k] / mid);
			}
		}
	}
	if (count(a[n - 1], a[n - 1] + n, 0) == n && a[n])
		return false;
	int loc = n - 1;
	for (int i = n - 1; i >= 0; --i) {
		res[i] = (a[i][n] / a[i][loc]);
		for (int j = i; j >= 0; --j) {
			a[j][n] -= (a[j][loc] * res[i]);
			a[j][loc] = 0;
		}
		loc--;
	}
	return true;
}

void Select_c(int k, int n) {
	double emax = a[k][k], loc = k;
	for (int i = k + 1; i < n; ++i)
		if (a[i][k] > emax) {
			emax = a[i][k];
			loc = i;
		}
	if (loc > k)
		for (int i = 0; i <= n; ++i)
			swap(a[k][i], a[(int)loc][i]);
}