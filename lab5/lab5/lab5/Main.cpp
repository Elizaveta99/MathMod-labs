#include<iostream>
#include<vector>
#include<map>
#include<algorithm>
#include<cmath>
#include<cstdlib>
#include<climits>
#include<iomanip>

using namespace std;

const unsigned int a_start = 65539, a_start_2 = 16387, M = 2147483648, GLOB_N = 1000, GLOB_M = 1000;
const double PI = 3.14159265358979323846, x_ok = 4.84732824427, y_ok = -1.29770992366, z_ok = -2.06106870229;

double Accuracy(vector<double> x, int N)
{
	double E2 = 0;
	for (int i = 0; i < N; i++)
		E2 += x[i];
	E2 /= N;

	double D2 = 0;
	for (int i = 0; i < N; i++)
		D2 += ((x[i] - E2) * (x[i] - E2));
	D2 /= (N - 1);

	return 0.6745 * sqrt(D2 / N);
}

int main()
{
	setlocale(LC_ALL, "Russian");

	for (int N = 10; N <= GLOB_N; N *= 10)
		for (int M = 10; M <= GLOB_M; M *= 10)
		{
			int n = 3;
			double x_ans = 0, y_ans = 0, z_ans = 0;
			vector<vector<double> > a(n + 1), p(n + 1);
			for (int i = 0; i < n; i++) 
			{
				a[i].resize(n + 1);
				p[i].resize(n + 1, 1.0 / double(n));
			}
			vector<double> pi(n + 1, 1.0 / double(n)), f(n + 1), h1(n + 1, 0),
				h2(n + 1, 0), h3(n + 1, 0), Q1(N + 1), Q2(N + 1, 0), Q3(N + 1, 0), 
				ksi1(M + 1, 0), ksi2(M + 1, 0), ksi3(M + 1, 0);
			vector<int> cp(N + 1);

			a[0][0] = 0.2, a[0][1] = -0.2, a[0][2] = -0.3, f[0] = 3.0;
			a[1][0] = 0.2, a[1][1] = 0.5, a[1][2] = 0.3, f[1] = -1.0;
			a[2][0] = -0.4, a[2][1] = -0.2, a[2][2] = -0.3, f[2] = -1.0;
			h1[0] = 1; h2[1] = 1; h3[2] = 1;

			for (int i = 0; i < M; i++)
			{
				double alpha = (double)rand() / double(RAND_MAX);
				if (alpha < pi[0])
					cp[0] = 0;
				else
					if (alpha >= pi[0] && alpha < 2 * pi[0]) cp[0] = 1;
						else cp[0] = 2;
				for (int j = 1; j <= N; j++)
				{
					double alpha = (double)rand() / double(RAND_MAX);
					if (alpha < 1.0 / double(n))
						cp[j] = 0;
					else
						if (alpha >= 1.0 / double(n) && alpha < 2.0 * 1.0 / double(n)) cp[j] = 1;
							else cp[j] = 2;

				}

				if (pi[cp[0]] > 0) {
					Q1[0] = h1[cp[0]] / pi[cp[0]];
					Q2[0] = h2[cp[0]] / pi[cp[0]];
					Q3[0] = h3[cp[0]] / pi[cp[0]];
				}
				else {
					Q1[0] = 0; Q2[0] = 0; Q3[0] = 0;
				}
				for (int j = 1; j <= N; j++)
				{
					if (p[cp[j - 1]][cp[j]] > 0) {
						Q1[j] = Q1[j - 1] * a[cp[j - 1]][cp[j]] / p[cp[j - 1]][cp[j]];
						Q2[j] = Q2[j - 1] * a[cp[j - 1]][cp[j]] / p[cp[j - 1]][cp[j]];
						Q3[j] = Q3[j - 1] * a[cp[j - 1]][cp[j]] / p[cp[j - 1]][cp[j]];
					}
					else {
						Q1[j] = 0; Q2[j] = 0; Q3[j] = 0;
					}
				}		
				
				for (int j = 0; j <= N; j++) {
					ksi1[i] += (Q1[j] * f[cp[j]]); 
					ksi2[i] += (Q2[j] * f[cp[j]]);
					ksi3[i] += (Q3[j] * f[cp[j]]);
				}
			}
			double x = 0, y = 0, z = 0;
			for (int i = 0; i < M; i++) {
				x += ksi1[i]; y += ksi2[i]; z += ksi3[i];
			}
			x /= M;  y /= M;  z /= M;

			cout << "N = " << N << " M = " << M << " " << fixed << setprecision(2) << x << " " << y << " " << z << "\n";
			cout << "Right value : " << x_ok << " " << y_ok << " " << z_ok << "\n";
			cout << "Accuracy of calculating : " << Accuracy(ksi1, M) << "\n";
		}

	system("pause");
	return 0;
}