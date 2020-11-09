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
		cout << "ÇëÊäÈëxµÄÖµ£º";
		cin >> num;
		for (int j = 0; j < 4; ++j) {
			double mid = 1;
			for (int k = 0; k < 4; ++k) {
				if (j != k)
					mid *= (num - loc[k]);
			}
			res += (a[loc[j]] * mid);
		}
		cout << "f(x) = " << res << endl;
	}
	return 0;
}
