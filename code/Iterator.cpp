#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

double **a = nullptr, * d = nullptr, * res = nullptr, maxk = 0;
int n = 0, num = 0;
bool Gauss_Seidel();
bool Jacobi();

int main() {
	cout << "���뷽��������������";
	cin >> num;
	cout << "����������ޣ�";
	cin >> maxk;
	cout << "���뷽�̽�����";
	cin >> n;
	a = new double* [n + 5];
	d = new double[n + 5] { 0 };
	res = new double[n + 5] { 0 };
	for (int i = 0; i < n + 5; ++i)
		a[i] = new double[n + 5] { 0 };
	cout << "���뷽�̾���" << endl;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j <= n; ++j) {
			cin >> a[i][j];
			if (i != j && j != n)
				a[i][j] = -a[i][j];
		}
	for (int i = 0; i < n; ++i)
		swap(d[i], a[i][i]);
	cout << "���������ʼ������";
	for (int i = 0; i < n; ++i)
		cin >> res[i];
	if (Gauss_Seidel())
		for (int i = 0; i < n; ++i)
			cout << "x" << i + 1 << " = " << res[i] << endl;
	else
		cout << "�����������ޣ�";
	cout << endl;
	delete[] d;
	delete[] res;
	for (int i = 0; i < n + 5; ++i)
		delete[] a[i];
	delete[] a;
	return 0;
}

bool Gauss_Seidel() {
	for (int i = 0; i < num; ++i) {
		double mid = 0, esp = 0;
		for (int j = 0; j < n; ++j) {
			mid = 0;
			for (int k = 0; k < n; ++k)
				mid += (a[j][k] * res[k]);
			mid += a[j][n];
			mid /= d[j];
			if (esp < fabs(res[j] - mid))
				esp = fabs(res[j] - mid);
			res[j] = mid;
		}
		if (esp <= maxk)
			return true;
	}
}

bool Jacobi() {
	for (int i = 0; i < num; ++i) {
		double *mid = nullptr, esp = 0;
		mid = new double[n + 5] { 0 };
		for (int j = 0; j < n; ++j) {
			for (int k = 0; k < n; ++k)
				mid[j] += (a[j][k] * res[k]);
			mid[j] += a[j][n];
			mid[j] /= d[j];
			if (esp < fabs(res[j] - mid[j]))
				esp = fabs(res[j] - mid[j]);
		}
		swap(res, mid);
		delete[] mid;
		if (esp <= maxk)
			return true;
	}
	return false;
}