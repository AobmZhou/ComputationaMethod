# 实验题目(五)：$ 解线性方程组的迭代法$

## 1. 实验目的

(1) 掌握雅可比迭代和Seidel迭代来求解方程组；
(2) 掌握常用的几种迭代格式；
(3) 编写程序实现上述迭代方法；
(4) 分析实验结果，并估计误差；

## 2. 实验环境

Windows 10 + Visual Studio

## 3. 实验内容

有如下线性方程组 Ax = b 如下：

$$
\left\{
\begin{matrix}
a_{11} & a_{12} &\cdots& a_{1n}\\
a_{21} & a_{22} &\cdots& a_{2n}\\
\vdots & \vdots & \ddots & \vdots\\
a_{n1} & a_{n2} &\cdots& a_{nn}
\end{matrix}
\right\}
\left\{
\begin{matrix}
x_1 \\
x_2\\
\vdots\\
x_n
\end{matrix}
\right\}
=\left\{
\begin{matrix}
b_1 \\
b_2\\
\vdots\\
b_n
\end{matrix}
\right\}
$$

使用迭代法进行求解，主要迭代方法为雅可比迭代和 Gauss-Seidel 迭代。

## 4. 实验程序

```c++
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

double **a = nullptr, * d = nullptr, * res = nullptr, maxk = 0;
int n = 0, num = 0;
bool Gauss_Seidel();
bool Jacobi();

int main() {
	cout << "输入方程最大迭代次数：";
	cin >> num;
	cout << "输入误差上限：";
	cin >> maxk;
	cout << "输入方程阶数：";
	cin >> n;
	a = new double* [n + 5];
	d = new double[n + 5] { 0 };
	res = new double[n + 5] { 0 };
	for (int i = 0; i < n + 5; ++i)
		a[i] = new double[n + 5] { 0 };
	cout << "输入方程矩阵：" << endl;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j <= n; ++j) {
			cin >> a[i][j];
			if (i != j && j != n)
				a[i][j] = -a[i][j];
		}
	for (int i = 0; i < n; ++i)
		swap(d[i], a[i][i]);
	cout << "输入迭代初始向量：";
	for (int i = 0; i < n; ++i)
		cin >> res[i];
	if (Gauss_Seidel())
		for (int i = 0; i < n; ++i)
			cout << "x" << i + 1 << " = " << res[i] << endl;
	else
		cout << "超出迭代上限！";
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
```

## 5. 实验结果

![Image text](https://gitee.com/aska-z/computationa-method/raw/master/picture/Iterator.PNG)

## 6. 实验总结

jacobi迭代Gauss_Seidel迭代没有更多优势，特别是在速度方面。但是良好的并行性质使得它在解方程组当中有了一席之地。Gauss_Seidel迭代其实是jacobi迭代的一种改进，在计算x的每一分量时，有新值不用旧值，使得迭代的速度有了一个质的飞跃。

# 实验题目(六)：$线性方程组的高斯消元法$

## 1. 实验目的

(1) 掌握高斯消元法求解方程组；
(2) 掌握列主元高斯消元法求解方程组；
(3) 分析实验结果，并估计误差；

## 2. 实验环境

Windows 10 + Visual Studio

## 3. 实验内容

使用高斯消元法求解下列方程组：

$$
\begin{cases}
10x_1-x_2-2x_3=7.2\\
-x_1+10x_2-2x_3=8.3\\
-x_1-x_2+5x_3=4.2
\end{cases}
$$

## 4. 实验程序

```c++
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
```

## 5. 实验结果

![Image text](https://gitee.com/aska-z/computationa-method/raw/master/picture/GaussianElimination.PNG)

## 6. 实验总结

高斯消元法主要用于求解线性方程组，也可以求矩阵的秩、矩阵的逆等，是一个重要的数学方法。其时间复杂度主要与方程组个数、方程组未知数个数有关，一般来说，时间复杂度为 O(n^3)。

# 实验题目(七)：$线性方程组的矩阵分解法$

## 1. 实验目的

(1) 掌握采用矩阵 LU 分解方法来求解线性方程组；
(2) 编程实现矩阵 LU 分解算法；

## 2. 实验环境

Windows 10 + Visual Studio

## 3. 实验内容

(1) 写出矩阵 LU 分解法解线性方程组算法，编一程序上机调试出结果，要求所编
程序适用于任何一解线性方程组问题，即能解决这一类问题，而不是某一个问
题。
(2) 使用矩阵 Doolittle 分解法求解下列方程组

$$
\begin{cases}
10x_1-x_2-2x_3=7.2\\
-x_1+10x_2-2x_3=8.3\\
-x_1-x_2+5x_3=4.2
\end{cases}
$$

## 4. 实验程序

```c++
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
```

## 5. 实验结果

![Image text](https://gitee.com/aska-z/computationa-method/raw/master/picture/Doolittle.PNG)

## 6. 实验总结

LU分解是矩阵分解的一种，可以将一个矩阵分解为一个单位下三角矩阵和一个上三角矩阵的乘积。LU分解是高斯消元的一种表达方式。首先，对矩阵A通过初等行变换将其变为一个上三角矩阵。然后，将原始矩阵A变为上三角矩阵的过程，对应的变换矩阵为一个下三角矩阵。这中间的过程，就是Doolittle算法。

# 实验题目(八)：$常微分方程求解算法$

## 1. 实验目的

(1) 掌握采用欧拉法来求解常微分方程；
(2) 掌握采用改进的欧拉法来求解常微分方程；
(3) 编程实现上述两个算法；

## 2. 实验环境

Windows 10 + Visual Studio

## 3. 实验内容

(1) 写出欧拉法或改进的欧拉法来求解常微分方程，编程序上机调试出结果。
(2) 使用常微分方程例子如下：

$$
\begin{cases}
y'=3x-2y^2-12&(0 < x < 5)\\
y(0)=2
\end{cases}
$$

## 4. 实验程序

```c++
#include<iostream>
#include<cmath>
#include<algorithm>
#define p(x) (x)*(x)
using namespace std;

double func(double x, double y) { return (3 * x) + (2 * p(y)) - 12; }

int main()
{
	double a = 0, b = 0, n = 0, y = 2;
	cout << "输入x的范围: ";
	cin >> a >> b;
	cout << "输入n: ";
	cin >> n;
	double h = (b - a) / n, x = a;
	for (double i = 1; i <= n; ++i) {
		double mid1 = y + h * func(x, y);
		double mid2 = y + h * func(x + h, mid1);
		y = (mid1 + mid2) / 2;
		cout << "x[" << i << "] = " << x + (h * i) << "	   " <<;
		cout << "y[" << i << "] = " << y << endl;
		x += h;
	}
	return 0;
}
```

## 5. 实验结果

![Image text](https://gitee.com/aska-z/computationa-method/raw/master/picture/EulerMethod.PNG)

## 6. 实验总结

可以将区间[a,b]分成n段，那么方程在第xi点有y'(xi)=f(xi,y(xi))，再用向前差商近似代替导数则为：(y(xi+1)-y(xi))/h=f(xi,y(xi))，在这里，h是步长，即相邻两个结点间的距离。因此可以根据xi点和yi点的数值计算出yi+1来。先用欧拉法求得一个初步的近似值，称为预报值，然后用它替代梯形法右端的yi+1再直接计算fi+1，得到校正值yi+1，这样称为改进的欧拉格式。