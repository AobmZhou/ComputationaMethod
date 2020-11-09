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
