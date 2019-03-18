#include<iostream>
#include<vector>
#include<map>
#include<algorithm>
#include<cmath>

using namespace std;

const unsigned int a_start = 65539, M = 2147483648, N = 1000;
const double PI = 3.14159265358979323846, E = -4, D = 4;
const double e = 0.05, k = 10;
double critical_value_kolmogorov = 1.36, critical_value_pirson = 16.92, lmbd = 0.5, A = 0, B = 1.5;
const double hb[5] = {0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};
const double p = 0.2316419;

//R[0, 1]
void mkm(vector<double> &a, unsigned int beta)
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

double func_exp(double x)
{
	return 1.0 - exp(-lmbd * x);
}

double func_logistic(double x)
{
	return 0.5 * (1.0 + tanh((x - A) / (2 * B)));
}

void Kolmogorov(vector<double> a, int type)
{
	sort(a.begin(), a.end());
	double D = 0;
	for (int i = 0; i < N; i++) 
	{
		if (type == 1)
			D = max(D, fabs(((double)(i + 1) / (double)N) - func_Gauss(a[i])));
		else 
			if (type == 2)
				D = max(D, fabs(((double)(i + 1) / (double)N) - func_exp(a[i])));
			else
				if (type == 3)
					D = max(D, fabs(((double)(i + 1) / (double)N) - func_logistic(a[i])));
		//cout << "D = " << D << ' ';
	}
	//cout << endl;

	//3.1
	double res_check = D * sqrt(N);
	double prob = abs((critical_value_kolmogorov - res_check) / critical_value_kolmogorov * 100.0);
	cout << "Вероятность ошибки I рода : " << prob << "\n";

	//5.3.1
	cout << "Check Kolmogorov\n" << res_check << "\n";
	if (res_check < critical_value_kolmogorov)
		cout << "ok\n";
	else cout << "failed\n";
}

void Pirson_Gauss(vector<double> a, int type)
{
	sort(a.begin(), a.end());
	int i = 0;
	int cnt = 0;
	double X = 0;
	int temp = 0, h = N / k - 1;
	double hh = (a[N - 1] - a[0]) / (N / k);
	for (int j = 1; j <= k; j++)
	{
		cnt = 0;
		while (a[i] <= a[h] && i < N)
		{
			i++;
			cnt++;
		}
		double pk;
		if (type == 1)
			pk = func_Gauss(a[h]) - func_Gauss(a[h - (N / k - 1)]);
		else 
			if (type == 2)
				pk = func_exp(a[h]) - func_exp(a[h - (N / k - 1)]);
			else 
				if (type == 3)
					pk = func_logistic(a[h]) - func_logistic(a[h - (N / k - 1)]);
		h += (N / k - 1);
		X += (((double)cnt - (double)N * pk) * ((double)cnt - (double)N * pk) / ((double)N * pk));
	}

	//4.1
	double res_check = X;
	double prob = abs((critical_value_pirson - res_check) / critical_value_pirson * 100.0);
	cout << "Вероятность ошибки I рода : " << prob << "\n";

	// 5.4.1
	cout << "Check Pirson\n" << res_check << "\n";
	if (res_check < critical_value_pirson)
		cout << "ok\n";
	else cout << "failed\n";
}

void E_D(vector<double> x, double E_ok, double D_ok)
{
	double E2 = 0;
	for (int i = 0; i < N; i++)
		E2 += x[i];
	E2 /= N;

	double D2 = 0;
	for (int i = 0; i < N; i++)
		D2 += ((x[i] - E2) * (x[i] - E2));
	D2 /= (N - 1);

	if (E2 < E_ok)
		cout << "Несмещённая оценка матожидания " << E2 << " меньше истинной " << E_ok << "\n";
	else cout << "Несмещённая оценка матожидания " << E2 << " больше истинной " << E_ok << "\n";
	if (D2 < D_ok)
		cout << "Несмещённая оценка дисперсии " << D2 << " меньше истинной " << D_ok << "\n";
	else cout << "Несмещённая оценка дисперсии " << D2 << " больше истинной " << D_ok << "\n";
}


int main()
{
	setlocale(LC_ALL, "Russian");
	vector<double> b(N * 12, 0), y(N, 0), x(N, 0), expc(N, 0), logt(N, 0);

	//R[0, 1]
	mkm(b, a_start);

	// 1 Gauss
	cout << "Нормальное распределение \n";

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

	E_D(x, E, D);

	Kolmogorov(x, 1);
	Pirson_Gauss(x, 1);

	//2 Exp
	cout << "Экспоненциальное распределение \n";

	for (int i = 0; i < N; i++)
		expc[i] = -1.0 / lmbd * log(b[i]);
	for (int i = 0; i < N; i++)
		cout << func_exp(expc[i]) << ' ';
	cout << "\n";

	double E_ok = 1.0 / lmbd, D_ok = 1.0 / (lmbd * lmbd);
	E_D(expc, E_ok, D_ok);

	Kolmogorov(expc, 2);
	Pirson_Gauss(expc, 2);

	//2 Logistic
	cout << "Логистическое распределение \n";

	for (int i = 0; i < N; i++)
		logt[i] = A + B * log((1.0 - b[i]) / b[i]) ; 
	for (int i = 0; i < N; i++)
		cout << func_logistic(logt[i]) << ' ';
	cout << "\n";

	E_ok = A, D_ok = 3.2899 * B * B;
	E_D(logt, E_ok, D_ok);

	Kolmogorov(logt, 3);
	Pirson_Gauss(logt, 3);

	system("pause");
	return 0;
}