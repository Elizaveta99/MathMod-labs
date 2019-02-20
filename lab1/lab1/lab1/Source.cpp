#include<iostream>
#include<vector>
#include<algorithm>

using namespace std;

const unsigned int a_start = 16387, M = 2147483648, N = 1000, K = 48;
const double e = 0.05, critical_value_kolmogorov = 1.36, critical_value_pirson = 18.307, k = 11;

void mkm(vector<double> &a, unsigned int beta)
{
	vector<unsigned int> az(N + K);
	az[0] = beta * beta; 
	a[0] = az[0] / (double)M;
	for (int i = 1; i < a.size(); i++)
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
	}
	return D /** sqrt(N)*/;
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
	vector<double> b(N + K, 0);
	mkm(b, a_start);
	
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
	mkm(c, a_start * 2 + 1); 
	res_check = Pirson(c);
	cout << "Check Pirson first\n" << res_check << "\n";
	if (res_check < critical_value_pirson)
		cout << "ok\n";
	else cout << "failed\n";

	vector<double> v(K);
	for (int i = 0; i < K; i++)
		v[i] = b[i]; 

	//McLaren - Marsaly
	vector<double> a(N + K, 0);
	for (int i = 0; i < N; i++)
	{
		int s = (int)(c[i] * (double)K);
		a[i] = v[s]; 
		v[s] = b[(i + K) % M];
	}
	
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