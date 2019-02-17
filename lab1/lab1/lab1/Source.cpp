#include<iostream>
#include<vector>

using namespace std;

const unsigned int a_start = 16387, M = 2147483648, N = 1000, K = 48;

void mkm(vector<double> &a, unsigned int beta)
{
	vector<unsigned int> az(N);
	az[0] = beta * beta; // ??? // beta * a_start
	a[0] = az[0] / (double)M;
	for (int i = 1; i < N; i++)
	{
		az[i] = (beta * az[i - 1]) % M;
		a[i] = az[i] / (double)M;
	}
}

int main()
{
	vector<double> b(N, 0);
	mkm(b, a_start);
	//cout << "Реализации базовых случайных величин с помощью мультипликативного конгруэнтного метода\n";
	cout << "First implementations of base random variables using a multiplicative congruential method\n";
	for (int i = 0; i < N; i++)
		cout << b[i] << ' ';
	cout << "\n";
	//check
	vector<double> c(N, 0);
	mkm(c, a_start * 4 + 1); // ???
	cout << "Second implementations of base random variables using a multiplicative congruential method\n";
	for (int i = 0; i < N; i++)
		cout << c[i] << ' ';
	cout << "\n";

	vector<unsigned int> v(K);
	for (int i = 0; i < K; i++)
		v[i] = b[i]; // ??

	//mmm
	vector<double> a(N, 0);
	for (int i = 0; i < N; i++)
	{
		unsigned int s = c[i] * K;
		a[i] = v[s]; // ??
		v[s] = b[i + K];
	}
	cout << "Third implementations of base random variables using the method of McLaren - Marsaly\n";
	for (int i = 0; i < N; i++)
		cout << a[i] << ' ';
	cout << "\n";
	//check
	return 0;
}