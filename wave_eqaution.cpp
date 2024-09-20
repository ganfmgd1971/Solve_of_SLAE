// wave_eqaution.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#define _USE_MATH_DEFINES // для C++
#include <iostream>
#include <vector>
#include <cmath>
#include<vector>
#include <fstream>

//std::vector<double> operator * (const std::vector<std::vector<double>>& A, const std::vector<double>& v) {
//
//	std::vector<double> tmp;
//	tmp.resize(v.size());
//
//	for (int i = 0; i < v.size(); i++) {
//
//		tmp[i] = 0;
//
//		for (int j = 0; j < v.size(); j++) {
//			tmp[i] += A[i][j] * v[j];
//		}
//
//	}
//
//	return tmp;
//}
//
//std::vector<double> operator * (const double c, const std::vector<double>& v) {
//
//	std::vector<double> tmp;
//
//	tmp.resize(v.size());
//
//	for (int i = 0; i < v.size(); i++) {
//		tmp[i] = c * v[i];
//	}
//
//	return tmp;
//}
//
//std::vector<double> operator + (const std::vector<double>& v1, const std::vector<double>& v2) {
//
//	std::vector<double> tmp;
//
//	tmp.resize(v1.size());
//
//	for (int i = 0; i < v1.size(); i++) {
//		tmp[i] = v1[i] + v2[i];
//	}
//
//	return tmp;
//}
//
//std::vector<double> operator - (const std::vector<double>& v1, const std::vector<double>& v2) {
//
//	std::vector<double> tmp;
//
//	tmp.resize(v1.size());
//
//	for (int i = 0; i < v1.size(); i++) {
//		tmp[i] = v1[i] - v2[i];
//	}
//
//	return tmp;
//}

double sqr(const std::vector<double>& v1, const std::vector<double>& v2) {

	double tmp = 0;

	for (int i = 0; i < v1.size(); i++) {
		tmp += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}

	return sqrt(tmp);
}

double maxnorm(const std::vector<double>& v1, const std::vector<double>& v2) {

	const int m = v1.size();
	double norm = 0;
	double tmp;

	for (int i = 0; i < m; i++) {
		tmp = abs(v1[i] - v2[i]);
		if (tmp > norm) {
			norm = tmp;
		}
	}

	return norm;
}

void printvector(std::vector<double>& a) {

	int n = a.size();

	for (int i = 0; i < n; i++) {
		std::cout << a[i] << "\n";
	}
	std::cout << "\n";
}

double phi(double t, const int key = 1) {

	switch (key) {

	case 1:
		return 0;
	
	case 2:
		return 2 * t * t;
	
	case 3:
		return 2 * t + 1;

	case 4:
		return sin(t);
	}
}

double psi(double t, const int key = 1) {

	switch (key) {

	case 1:
		return 0;

	case 2:
		return 1;
	}
}

double f(double x, const int key = 1) {

	switch (key) {

	case 1:
		return sin(M_PI * x);

	case 2:
		return x * (1 - x);

	case 3:
		return 0.5 * x * (x + 1);

	case 4:
		return (1 - x) * cos(M_PI * x * 0.5);

	case 5:
		if ((x > -1 / 3) && (x < 1 / 3)) {
			return 1;
		}

		else {
			return 0;
		}

	case 6:
		return 0;

	}
}

double ddf(double x, const int key = 1) {

	switch (key) {

	case 1:
		return -M_PI * M_PI * sin(M_PI * x);

	case 2:
		return -2;

	case 3:
		return 1;
	
	case 4:
		return -M_PI * M_PI * 0.25 * (1 - x) * cos(M_PI * x * 0.5) + M_PI * sin(M_PI * x * 0.5);
	}

}

double g(double x, const int key = 1) {

	switch (key) {

	case 1:
		return 0;

	case 2:
		return x * cos(x);

	case 3:
		return 2 * x + 1;

	case 4:
		if ((x > -0.5) && (x < 0.5)) {
			return 1 - 2 * abs(x);
		}

		else {
			return 0;
		}

	}
}

void swap(std::vector<double> &a, std::vector<double> &b, const std::vector<double> &c) {
	
	a = b;
	b = c;
}

class waveequation {

private:

	double T; //длина отрезка по времени
	double L; //длина отрезка по пространству
	double R;
	int key_init1; //ключ на начальное условие на функцию
	int key_init2; //ключ на начальное условие на производную по времени
	int key_bound1; //ключ на левое граничное условие
	int key_bound2; //ключ на правое граничное условие
	double a;

public:

	waveequation(double T0, double L0, double R0, int key_init1tmp, int key_init2tmp, int key_bound1tmp, int key_bound2tmp, double aa = 1) {

		T = T0; //длина отрезка по времени
		L = L0; //левый конец отрезка
		R = R0; //правый конец отрезка
		key_init1 = key_init1tmp; //ключ на начальное условие на функцию
		key_init2 = key_init2tmp; //ключ на начальное условие на производную по времени
		key_bound1 = key_bound1tmp; //ключ на левое граничное условие
		key_bound2 = key_bound2tmp; //ключ на правое граничное условие
		a = aa;

	}
	
	void cross_scheme(int N, int M, bool flag = true) { //N - по времени, М - по оси икс.

		double h = (R - L) / M;
		double tau = T / N;
		std::ofstream out;
		std::vector<double> yprev(M + 1, 0);
		std::vector<double> ynext(M + 1, 0);
		std::vector<double> y(M + 1, 0);
		double x = L;
		double t = 0;
		double gamma = (tau * a / h);
		double tmp;

		yprev[0] = phi(0, key_bound1);
		yprev[M] = psi(0, key_bound2);
		y[0] = phi(tau, key_bound1);
		y[M] = psi(tau, key_bound2);

		for (int i = 1; i < M; i++) {

			x = i * h;
			yprev[i] = f(x, key_init1);

			if (flag) {
				tmp = ddf(x, key_init1);
			}

			else {
				tmp = (f(x + h, key_init1) - 2 * f(x, key_init1) + f(x - h, key_init1)) / (h * h);
			}

			y[i] = yprev[i] + tau * g(x, key_init2) + a * a * tau * tau * tmp / 2;

		}
		
		out.open("D:\\vichy_umf_output_laba3.txt");

		if (out.is_open()) {

			std::cout << gamma << "  " << flag << "\n";
			out << tau << "\n";
			out << h << "\n";
			out << L << "\n";

			//std::cout << t << std::endl;

			out << t << " ";
			for (int i = 0; i < M + 1; i++) {
				out << yprev[i] << " ";
			}
			out << "\n";

			t += tau; 

			out << t << " ";
			for (int i = 0; i < M + 1; i++) {
				out << y[i] << " ";
			}
			out << "\n";

			for (int i = 2; i < N + 1; i++) {

				t = i * tau;

				ynext[0] = phi(t, key_bound1);
				ynext[M] = psi(t, key_bound2);

				for (int j = 1; j < M; j++) {

					ynext[j] = gamma * gamma * y[j + 1] + 2 * (1 - gamma * gamma) * y[j] + gamma * gamma * y[j - 1] - yprev[j];

				}

				/*yprev = y;
				y = ynext;*/
				swap(yprev, y, ynext);

				out << t << " ";
				for (int k = 0; k < M + 1; k++) {
					out << ynext[k] << " ";
				}
				out << "\n";
			}

		}

		out.close();
	}

	double value_cross_scheme(int N, int M) { //N - по времени, М - по оси икс.

		double h = (R - L) / M;
		double tau = T / N;
		std::ofstream out;
		std::vector<double> yprev(M + 1, 0);
		std::vector<double> ynext(M + 1, 0);
		std::vector<double> y(M + 1, 0);
		double x = L;
		double t = 0;
		double gamma = (tau * a / h);
		double tmp;
		double tmpnorm;
		double norm = 0;
		double dif;

		yprev[0] = phi(0, key_bound1);
		yprev[M] = psi(0, key_bound2);
		y[0] = phi(tau, key_bound1);
		y[M] = psi(tau, key_bound2);

		for (int i = 1; i < M; i++) {

			x = i * h;
			yprev[i] = f(x, key_init1);
			tmp = (f(x + h, key_init1) - 2 * f(x, key_init1) + f(x - h, key_init1)) / (h * h);
			y[i] = yprev[i] + tau * g(x, key_init2) + a * a * tau * tau * tmp / 2;

		}

		//std::cout << gamma << "  " << "\n";

		//std::cout << t << std::endl;

		t += tau;

		for (int i = 2; i < N + 1; i++) {

			t = i * tau;

			ynext[0] = phi(t, key_bound1);
			ynext[M] = psi(t, key_bound2);

			for (int j = 1; j < M; j++) {

				ynext[j] = gamma * gamma * y[j + 1] + 2 * (1 - gamma * gamma) * y[j] + gamma * gamma * y[j - 1] - yprev[j];

			}

			/*yprev = y;
			y = ynext;*/
			swap(yprev, y, ynext);

			if (t <= 0.3) {

				tmpnorm = 0;
				for (int i = 0; i < M + 1; i++) {
					x = i * h;
					dif = abs(ynext[i] - sin(M_PI * x) * cos(M_PI * t));
					if (dif > tmpnorm) {
						tmpnorm = dif;
					}
				}
				if (tmpnorm > norm) {
					norm = tmpnorm;
				}

			}
		}
		
		return norm;
	}
}; 

//tests
//test 1: keyinit1 - 1, keyinit2 - 1, keybound1 - 1, keybound2 - 1. 
//test 2: keyinit1 - 2, keyinit2 - 1, keybound1 - 1, keybound2 - 1. 
//test 3: keyinit1 - 3, keyinit2 - 2, keybound1 - 2, keybound2 - 2. 
//test 4: keyinit1 - 4, keyinit2 - 3, keybound1 - 3, keybound2 - 1.
//test 5: keyinit1 - 5, keyinit2 - 1, keybound1 - 1, keybound2 - 1.
//test 6: keyinit1 - 6, keyinit2 - 4, keybound1 - 1, keybound2 - 1.
//test 7: keyinit1 - 6, keyinit2 - 1, keybound1 - 4, keybound2 - 1.
//test 8: keyinit1 - 1, keyinit2 - 1, keybound1 - 1, keybound2 - 1.  для порядка
void test1() {

	//первый аргумент - по времени, второй аргумент - по оси икс.

	waveequation umfequation(10, 0, 10, 1, 1, 1, 1);
	umfequation.cross_scheme(100, 100, true);

}

void test2() {
	
	//первый аргумент - по времени, второй аргумент - по оси икс.

	waveequation umfequation(1, 0, 1, 2, 1, 1, 1);
	umfequation.cross_scheme(100, 100, true);
}

void test3() {

	//первый аргумент - по времени, второй аргумент - по оси икс.

	waveequation umfequation(1, 0, 1, 4, 3, 3, 1);
	umfequation.cross_scheme(100, 50);

}

void test4() {

	//первый аргумент - по времени, второй аргумент - по оси икс.

	waveequation umfequation(1, 0, 1, 4, 3, 3, 1);
	umfequation.cross_scheme(100, 50);

}

void test5() {

	//первый аргумент - по времени, второй аргумент - по оси икс.

	waveequation umfequation(100, -2, 2, 5, 1, 1, 1);
	umfequation.cross_scheme(10000, 200, false);

}

void test6() {

	waveequation umfequation(2, -1, 1, 5, 1, 1, 1);
	umfequation.cross_scheme(100, 100, false);

}

void test7() {

	waveequation umfequation(10, 0, 4 * M_PI, 6, 1, 4, 1);
	umfequation.cross_scheme(1000, 1200, true);

}

void test8() {
	
	waveequation umfequation(1, 0, 1, 1, 1, 1, 1);
	double d1;
	double d2 = 1;
	double n = 10;
	double m = 5;

	for (int i = 0; i < 6; i++) {
		
		d2 = umfequation.value_cross_scheme(n, m);
		if (i == 0) {
			std::cout << d2 << " ----_ -----" << "\n";
		}
		else {
			std::cout << d2 << " " << d1 / d2 << " " << log(d1 / d2) / log(2) << "\n";
		}
		d1 = d2;
		n *= 2;
		m *= 2;
	}
}
int main()
{
	//test1();
	//test2();
	//test3();
	//test4();
	test5();
	//test6();
	//test7();
	//test8();
}