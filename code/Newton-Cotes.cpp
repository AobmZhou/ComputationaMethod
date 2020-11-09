#include<iostream>
#include<cmath>
using namespace std;

double c[8][8 + 2] = { 
	{0},
	{2, 1, 1}, 
	{6, 1, 4, 1}, 
	{8, 1, 3, 3, 1}, 
	{90, 7, 32, 12, 32, 7}, 
	{288, 19, 75, 50, 50, 75, 19}, 
	{840, 41, 216, 27, 272, 27, 216, 41},
	{17280, 751, 3577, 1323, 2989, 2989, 1323, 3577, 751} 
};

int main()
{
	double a = 0, b = 0, n = 0;
	double res[3] = { 0 }, y[3][8] = { 0 };
	cout << "ÊäÈëa, b£º ";
	cin >> a >> b;
	cout << "ÊäÈën£º";
	cin >> n;
	double h = (b - a) / n, x = a;
	for (int i = 1; i <= (int)n + 1; ++i) {
		if (!x)
			y[0][i] = 1;
		else
			y[0][i] = sin(x) / x;
		y[1][i] = exp(x) / (4 + (x * x));
		y[2][i] = log(1 + x) / (1 + (x * x));
		x += h;
	}
	for (int i = 1; i <= (int)n + 1; ++i) {
		res[0] += y[0][i] * c[(int)n][i];
		res[1] += y[1][i] * c[(int)n][i];
		res[2] += y[2][i] * c[(int)n][i];
	}
	for (int i = 0; i < 3; ++i)
		res[i] *= (b - a) / c[(int)n][0];
	for (int i = 0; i < 3; ++i) {
		x = a;
		for (int j = 1; j <= n + 1; ++j, x += h)
			cout << "f(" << x << ") = " << y[i][j] << " " << endl;
		cout << "res = " << res[i] << endl;
		cout << endl;
	}
	return 0;
}
