#include<iostream>
#include<cmath>
#include<algorithm>
#define p(x) (x)*(x)
using namespace std;

double func(double x, double y) { return (3 * x) + (2 * p(y)) - 12; }

int main()
{
	double a = 0, b = 0, n = 0, y = 2;
	cout << "ÊäÈëxµÄ·¶Î§: ";
	cin >> a >> b;
	cout << "ÊäÈën: ";
	cin >> n;
	double h = (b - a) / n, x = a;
	for (double i = 1; i <= n; ++i) {
		double mid1 = y + h * func(x, y);
		double mid2 = y + h * func(x + h, mid1);
		y = (mid1 + mid2) / 2;
		cout << "x[" << i << "] = " << x + (h * i) << "	   " << "y[" << i << "] = " << y << endl;
		x += h;
	}
	return 0;
}