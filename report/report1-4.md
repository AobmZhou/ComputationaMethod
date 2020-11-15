# 实验题目(一)：$Lagrange插值方法$

## 1. 实验目的

(1) 熟悉简单的一阶和二阶Lagrange插值方法；
(2) 学会计算Lagrange基函数；
(3) 正确构造插值多项式；
(4) 对插值结果进行合理分析；

## 2. 实验环境

Windows 10 + Visual Studio

## 3. 实验内容

设函数为 $f(x) = 5\sqrt{x}$，已知
| &nbsp;$x$&nbsp;| $f(x)$     |
| :---:  | :-----:  |
| 24     | 1.888175 |
| 26     | 1.918645 |
| 28     | 1.947294 |
| 30     | 1.961009 |

使用 $Lagrange$ 插值多项式计算 $f(25), f(27), f(29)$，并给出插值多项式。修改程序
直至运行成功，查看运行结果，并和如下真实值进行比较。
|&nbsp;&nbsp;&nbsp;&nbsp; $x$ &nbsp;&nbsp;&nbsp;&nbsp; |    $f(x)$   |
| :------: |  :-----:  |
|    25    | 1.903653938715879 |
|    27    | 1.933182044931763 |
|    29    | 1.961009057454548 |

## 4. 实验程序

```c++
#include<iostream>
#include<cmath>
using namespace std;
double a[35] = { 0 };

int main()
{
	int loc[4] = { 24, 26, 28, 30 };
	a[24] = 1.888175 / (-48);
	a[26] = 1.918645 / (16);
	a[28] = 1.947294 / (-16);
	a[30] = 1.961009 / (48);
	double num = 0, res = 0;
	for (int i = 0; i < 3; ++i) {
		res = 0;
		cout << "请输入x的值：";
		cin >> num;
		for (int j = 0; j < 4; ++j) {
			double mid = 1;
			for (int k = 0; k < 4; ++k) {
				if (j != k)
					mid *= (num - loc[k]);
			}
			res += (a[loc[j]] * mid);
		}
		cout << res << endl;
	}
	return 0;
}
```

## 5. 实验结果

![Image text](https://gitee.com/aska-z/computationa-method/raw/master/picture/Lagrange.PNG)

## 6. 实验总结

通过这次实验，使我对Lagrange插值方法有了更深刻的了解。拉格朗日插值法的公式结构整齐紧凑，在理论分析中十分方便，然而在计算中，当插值点增加或减少一个时，所对应的基本多项式就需要全部重新计算，于是整个公式都会变化，非常繁琐。这时可以用重心拉格朗日插值法或牛顿插值法来代替。


# 实验题目(二)：$牛顿插值方法$

## 1. 实验目的

(1) 理解牛顿插值方法；
(2) 学会计算差商；
(3) 正确构造插值多项式；
(4) 设计程序并调试得到正确结果；

## 2. 实验环境

Windows 10 + Visual Studio

## 3. 实验内容

样本数据值如下：
| &nbsp;&nbsp;&nbsp;$x$&nbsp;&nbsp;&nbsp; |0.4|0.55|0.65|0.8|0.9|
|:---:|:---:|:---:|:---:|:---:|:---:|
| $f(x)$|0.41075|0.57815|0.69675|0.88811|1.02652|

## 4. 实验程序

```c++
#include<iostream>
#include<cmath>
#include<vector>
using namespace std;

class loc {
public:
	double x = 0, y = 0;
	loc(double x = 0, double y = 0) { this->x = x; this->y = y; }
};

int main()
{
	double num = 0.895, res = 0;
	vector<loc> v;
	v.push_back(loc(0.4, 0.41075));
	v.push_back(loc(0.55, 0.57815));
	v.push_back(loc(0.65, 0.69675));
	v.push_back(loc(0.8, 0.88811));
	v.push_back(loc(0.9, 1.02652));
	vector<vector<double>> dc(v.size());
	for (int i = 0; i < dc.size(); ++i)
		dc[i].push_back(v[i].y);
	for (int i = 1; i < dc.size(); ++i) {
		for (int j = 1; j <= i; ++j) {
			double mid = 0;
			mid = dc[i - 1][j - 1] - dc[i][j - 1];
			mid /= (v[i-j].x-v[i].x);
			dc[i].push_back(mid);
		}
	}
	for (int i = 0; i < dc.size(); ++i) {
		double mid = 1;
		for (int j = 0; j < i; ++j) {
			mid *= (num - v[j].x);
		}
		res += (mid * dc[i][i]);
	}
	cout << "f(x) = " << res << endl;
	return 0;
}
```

## 5. 实验结果

![Image text](https://gitee.com/aska-z/computationa-method/raw/master/picture/NewtonInterpolation.PNG)

## 6. 实验总结

通过实验结果来看，牛顿插值法的特点在于：每增加一个点，不会导致之前的重新计算，只需要算和新增点有关的就可以了。牛顿插值法相对于拉格朗日插值法具有承袭性的优势，即在增加额外的插值点时，可以利用之前的运算结果以降低运算量。

# 实验题目(三)：$Newton-Cotes方法$

## 1. 实验目的

(1) 掌握Newton-Cotes算法；
(2) 要求程序不断加密对积分区间的等分，自动地控制Newton-Cotes算法中的加
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;速收敛过程；
(3) 编写程序，分析实验结果；

## 2. 实验环境

Windows 10 + Visual Studio

## 3. 实验内容

计算以下积分值：
(1)&nbsp;$I$&nbsp;$=$&nbsp;$\int^0_1\frac{sinx}{x}{\rm d}x$&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;($f(0)$&nbsp;$=$&nbsp;$1$,&nbsp;&nbsp;$I$&nbsp;$\approx$ 0.9460831)
(2)&nbsp;$I$&nbsp;$=$&nbsp;$\int^0_1\frac{e^x}{4+x^2}{\rm d}x$
(3)&nbsp;$I$&nbsp;$=$&nbsp;$\int^0_1\frac{\ln(1+x)}{1+x^2}{\rm d}x$

## 4. 实验程序

```c++
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
	cout << "输入a, b： ";
	cin >> a >> b;
	cout << "输入n：";
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
```

## 5. 实验结果

![Image text](https://gitee.com/aska-z/computationa-method/raw/master/picture/Newton-Cotes.PNG)

## 6. 实验总结

Newton-Cotes方法通过以函数于等距n+1点的值，取得一个n次的多项式来近似原来的函数，再行求积。不过从结果来看，对于次数较高的多项式而有很大误差。

# 实验题目(四)：$求非线性方程根的牛顿法$

## 1. 实验目的

(1) 掌握求非线性方程根的牛顿法；
(2) 进一步了解牛顿法的改进算法；
(3) 编写程序，分析实验结果；

## 2. 实验环境

Windows 10 + Visual Studio

## 3. 实验内容

用牛顿迭代法求$xe^x-1=0$的根，迭代初始值为$x_0=0.5$。

## 4. 实验程序

```c++
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
```

## 5. 实验结果

![Image text](https://gitee.com/aska-z/computationa-method/raw/master/picture/NewtonMethodOfNonlinearEquations.PNG)

## 6. 实验总结

牛顿法亦称作切线法。一般非线性方程组的根通常无法直接求解，因此需要使用如牛顿法一类的迭代算法求近似解。通过实验来看，牛顿法只能逼近解，不能计算精确解。不过实际应用中，我们都不要求绝对精确的解，只要精度足够高就好了。