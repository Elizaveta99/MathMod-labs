#include<iostream>
#include<vector>
#include<map>
#include<algorithm>
#include<cmath>

using namespace std;

const unsigned int a_start = 65539, M = 2147483648, GLOB_N = 1000000;
const double PI = 3.14159265358979323846, E = 0, D = 1, integral_math = 0.22721;

//R[0, 1]
void mkm(vector<double> &a, unsigned int beta, int N)
{
	vector<unsigned int> az(N * 12);
	az[0] = beta * beta;
	a[0] = az[0] / (double)M;
	for (int i = 1; i < a.size(); i++)
	{
		az[i] = (beta * az[i - 1]) % M;
		a[i] = az[i] / (double)M;
	}
}

double func_Gauss(double x)
{
	return 0.5 * (1.0 + erf((x - E) / sqrt(2 * D)));
}

double func(double x)
{
	double res = 1.0;
	res /= (x*x*x*x + 3.0*x*x + 17.0);
	return res;
}


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

	return 3.0 * sqrt(D2 / N);
}

int main()
{
	setlocale(LC_ALL, "Russian");

	for (int N = 10; N <= GLOB_N; N *= 10)
	{
		vector<double> b(N * 12, 0), x(N, 0);
		//R[0, 1]
		mkm(b, a_start, N);

		int cnt = 0;
		for (int i = 0; i < N * 12; i += 12)
		{
			for (int j = i; j < i + 12; j++)
				x[cnt] += b[j];
			x[cnt] -= 6;
			x[cnt] *= sqrt(D);
			x[cnt] += E;
			cnt++;
		}

		double integral = sqrt(2.0 * PI) / (double)N;
		double sum = 0;
		for (int i = 0; i < N; i++)
		{
			sum += (func(x[i]) * exp(x[i] * x[i] / 2.0));
		}
		integral *= sum;

		cout << "N = " << N << "\n";
		cout << "Approximate value of integral : " << integral << "\n";
		cout << "Value of integral from Mathematica : " << integral_math << "\n";
		cout << "Accuracy of calculating : " << Accuracy(x, N) << "\n";
	}

	system("pause");
	return 0;
}