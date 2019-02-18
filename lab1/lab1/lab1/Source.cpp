#include<iostream>
#include<vector>
#include<algorithm>

using namespace std;

const unsigned int a_start = 16387, M = 2147483648, N = 1000, K = 48;
const double e = 0.05, critical_value_kolmogorov = 1.36, critical_value_pirson = 18.307, k = 10;

void mkm(vector<double> &a, unsigned int beta)
{
	vector<unsigned int> az(N);
	az[0] = beta * beta; 
	a[0] = az[0] / (double)M;
	for (int i = 1; i < N; i++)
	{
		az[i] = (beta * az[i - 1]) % M;
		a[i] = az[i] / (double)M;
	}
}

double Kolmogorov(vector<double> a)
{
	sort(a.begin(), a.end());
	double D = 0;
	for (int i = 0; i < N; i++) {
		D = max(D, fabs(((double)(i + 1) / (double)N) - a[i]));
		cout << "!!!! = " << fabs(((double)(i + 1) / (double)N) - a[i]) << "\n";
	}
	return D /** (double)N*/;
}

double Pirson(vector<double> a)
{
	sort(a.begin(), a.end());
	double i = 1.0 / k;
	int j = 0, cnt = 0;
	double X = 0;
	while (i <= 1)
	{
		while (j < N && a[j] < i)
		{
			j++;
			cnt++;
		}
		X += (((double)cnt - (double)N / k) * ((double)cnt - (double)N / k) / ((double)N / k));
		cnt = 0;
		i += (1.0 / (double)k);
	}
	return X;
}

int main()
{
	vector<double> b(N, 0);
	mkm(b, a_start);
	/*cout << "First implementations of base random variables using a multiplicative congruential method\n";
	for (int i = 0; i < N; i++)
		cout << b[i] << ' ';
	cout << "\n";*/

	double res_check = Kolmogorov(b);
	cout << "Check Kolmogorov first\n" << res_check << "\n";
	if (res_check < critical_value_kolmogorov)
		cout << "ok\n";
	else cout << "failed\n";

	res_check = Pirson(b);
	cout << "Check Pirson first\n" << res_check << "\n";
	if (res_check < critical_value_pirson)
		cout << "ok\n";
	else cout << "failed\n";

	vector<double> c(N, 0);
	mkm(c, a_start * 4 + 1); 
	/*cout << "Second implementations of base random variables using a multiplicative congruential method\n";
	for (int i = 0; i < N; i++)
		cout << c[i] << ' ';
	cout << "\n";*/

	vector<double> v(K);
	for (int i = 0; i < K; i++)
		v[i] = b[i]; 

	//McLaren - Marsaly
	vector<double> a(N, 0);
	for (int i = 0; i < N; i++)
	{
		unsigned int s = c[i] * K / M;
		a[i] = v[s]; 
		if (i + K < N)
			v[s] = b[i + K];
		else v[s] = b[i];
	}
	/*cout << "Third implementations of base random variables using the method of McLaren - Marsaly\n";
	for (int i = 0; i < N; i++)
		cout << a[i] << ' ';
	cout << "\n";*/

	res_check = Kolmogorov(a);
	cout << "Check Kolmogorov second\n" << res_check << "\n";
	if (res_check < critical_value_kolmogorov)
		cout << "ok\n";
	else cout << "failed\n";

	res_check = Pirson(a);
	cout << "Check Pirson second\n" << res_check << "\n";
	if (res_check < critical_value_pirson)
		cout << "ok\n";
	else cout << "failed\n";
	return 0;
}