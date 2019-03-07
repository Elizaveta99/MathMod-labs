#include<iostream>
#include<vector>
#include<map>
#include<algorithm>

using namespace std;

const unsigned int a_start = 65539, M = 2147483648, N = 1000, K = 48, m = 6, r = 4;
const double e = 0.05,
				delta_pirson1 = 16.919 /*11.071*/, // ?? кол-во степ свободы = 6 (6 - 1 = 5) ??
				delta_pirson2 = 16.919 /*7.815*/, // 
				k = 10, // change
				p_binom = 0.3333333, p_pascal = 0.2; // delta change

vector<int> kol;
int c[100][100];

void mkm(vector<double> &a, unsigned int beta, int size)
{
	vector<unsigned int> az(size + 1);
	az[0] = beta * beta;
	a[0] = az[0] / (double)M;
	for (int i = 1; i < a.size(); i++)
	{
		az[i] = (beta * az[i - 1]) % M;
		a[i] = az[i] / (double)M;
	}
}


double Pirson_binom(vector<int> a) // Binom
{
	for (int i = 0; i < a.size(); i++)
		kol[a[i]]++; // или это количество до i-ого, равных i-тому
	//for (int i = 0; i < a.size(); i++)
		//cout << kol[a[i]] << ' '; cout << endl;
	double X = 0;
	for (int ii = 0; ii < N; ii++)
	{
		double pi = (double)c[m][a[ii]] * pow(p_binom, (double)a[ii]) * pow(1 - p_binom, (double)(m - a[ii]));
		double deg1 = 1.0, deg2 = 1.0;
		for (int i = 0; i < a[ii]; i++)
			deg1 *= p_binom;
		for (int i = 0; i < (m - a[ii]); i++)
			deg2 *= (1.0 - p_binom);
		//cout << "pow = " << pi << ' ' << (int)c[m][a[ii]] << ' ' << a[ii] << ' ' << deg1 << ' ' << deg2 << endl;
		X += (((((double)kol[a[ii]] - (double)N * pi) * ((double)kol[a[ii]] - (double)N * pi))) / (double)N / pi);
		//cout << "!!! X = " << X << ' ' << (((((double)kol[a[ii]] - (double)N * pi) * ((double)kol[a[ii]] - (double)N * pi))) / (double)N / pi) << endl;
	}
	return X;
}

double Pirson_neg_binom(vector<int> a) // Pascal
{
	

	return 0;
}

int main()
{
	setlocale(LC_ALL, "Russian");

	// сочетания
	for (int i = 0; i <= 25; i++) // <= ???
		c[i][0] = 1,
		c[i][i] = 1;
	for (int i = 2; i <= 25; i++)
		for (int j = 1; j <= 25; j++)
			c[i][j] = c[i - 1][j] + c[i - 1][j - 1]; // ?? правильность проверить


	for (int i = 0; i <= 25; i++)
	{
		for (int j = 0; j <= 25; j++)
			cout << c[i][j] << ' ';
		cout << endl;
	}
	cout << endl;

	kol.resize(m + 1, 0);
	//system("pause");
	//return 0;


	cout << "Биномиальное распределение : \n";
	// 1
	vector<double> b(m * N + 1, 0);
	vector<int> x(N + 1, 0);
	mkm(b, a_start, m * N);
	int X = 0, cnt = 0;
	for (int i = 0; i < m * N; i += m)
	{
		for (int j = i; j < i + m; j++)
			x[cnt] += (p_binom > b[j] ? 1 : 0);
		cnt++;
	}

	for (int i = 0; i < N; i++)
		cout << x[i] << ' '; cout << endl;
	// 2
	double E = (double)(m) * p_binom;
	double D = (double)(m) * p_binom * (1.0 - p_binom);

	double E1 = 0;
	for (int i = 0; i < N; i++)
		E1 += x[i];
	E1 /= N;

	double D1 = 0;
	for (int i = 0; i < N; i++)
		D1 += ((x[i] - E1) * (x[i] - E1));
	D1 /= N;

	if (E1 < E)
		cout << "Несмещённая оценка матожидания " << E1 << " меньше истинной " << E << "\n";
	else cout << "Несмещённая оценка матожидания " << E1 << " больше истинной " << E << "\n";
	if (D1 < D)
		cout << "Несмещённая оценка матожидания " << D1 << " меньше истинной " << D << "\n";
	else cout << "Несмещённая оценка матожидания " << D1 << " больше истинной " << D << "\n";

	// 3
	double res_check;
	res_check = Pirson_binom(x);
	double prob = abs((delta_pirson1 - res_check) / delta_pirson1 * 100.0); // ??
	cout << "Вероятность ошибки I рода : " << prob << "\n";

	// 4
	cout << "Check Pirson first\n" << res_check << "\n";
	if (res_check < delta_pirson1)
		cout << "ok\n";
	else cout << "failed\n";




	cout << "Отрицательное биномиальное распределение : \n";
	// 1
	b.resize(r * 10 * N + 1, 0);
	x.resize(N + 1, 0);
	mkm(b, a_start, r * 10 * N);
	X = 0, cnt = 0;
	int i = 0, nu = 0, u = 0;
	while (i < r * 10 * N && cnt < N)
	{
		u = 0;
		while (i < r * 10 * N && u <= r)
		{
			if (p_pascal < b[i])
				x[cnt]++;
			else u++;
			i++;
		}
		cnt++;
	}

	cout << "i = " << i << endl;

	/*system("pause");
	return 0;*/

	
	cout << cnt - 1 << "\n";
	for (int i = 0; i < N; i++)
		cout << x[i] << ' '; cout << endl;

	// 2
	E = (double)(r) * (1.0 - p_pascal) / p_pascal;
	D = E / p_pascal;

	E1 = 0;
	for (int i = 0; i < N; i++)
		E1 += x[i];
	E1 /= N;

	D1 = 0;
	for (int i = 0; i < N; i++)
		D1 += ((x[i] - E1) * (x[i] - E1));
	D1 /= N;

	if (E1 < E)
		cout << "Несмещённая оценка матожидания " << E1 << " меньше истинной " << E << "\n";
	else cout << "Несмещённая оценка матожидания " << E1 << " больше истинной " << E << "\n";
	if (D1 < D)
		cout << "Несмещённая оценка матожидания " << D1 << " меньше истинной " << D << "\n";
	else cout << "Несмещённая оценка матожидания " << D1 << " больше истинной " << D << "\n";

	// 3
	res_check;
	res_check = Pirson_neg_binom(x);
	prob = (delta_pirson2 - res_check) / delta_pirson2 * 100; // ??
	cout << "Вероятность ошибки I рода : " << prob << "\n";

	// 4
	cout << "Check Pirson first\n" << res_check << "\n";
	if (res_check < delta_pirson2)
		cout << "ok\n";
	else cout << "failed\n";

	system("pause");
	return 0;
}