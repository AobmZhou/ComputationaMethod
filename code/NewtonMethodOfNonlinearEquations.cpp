#include<iostream>
#include<cmath>
using namespace std;

bool func(double & x, double num, double eps);
double f1(double x) { return exp(x) + x * exp(x); }
double f(double x) { return x * exp(x) - 1; }

int main()
{
	double x0 = 0, num = 0, eps = 0;
	cout << "输入迭代初始值：";
	cin >> x0;
	cout << "输入迭代最大次数：";
	cin >> num;
	cout << "输入精度要求：";
	cin >> eps;
	if (!func(x0, num, eps))
		cout << "超过迭代上限！" << endl;
	else
		cout << "x = " << x0 << endl;
	return 0;
}

bool func(double & x, double num, double eps)
{
	double xi = 0;
	for (int i = 0; i < num; ++i) {
		xi = x;
		x = xi - (f(xi) / f1(xi));
		if (fabs(xi - x) <= eps)
			return true;
	}
	return false;
}
