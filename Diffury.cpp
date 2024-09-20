#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

const double omega = 2.0;
const double sigma = 10;
const double r = 28;
const double b = 8 / 3;
vector<double> operator * (const vector<vector<double>> A, const vector<double> v) {
	vector<double> tmp;
	tmp.resize(v.size());
	for (int i = 0; i < v.size(); i++) {
		tmp[i] = 0;
		for (int j = 0; j < v.size(); j++) {
			tmp[i] += A[i][j] * v[j];
		}
	}
	return tmp;
}
vector<double> operator * (const double c, const vector<double> v) {
	vector<double> tmp;
	tmp.resize(v.size());
	for (int i = 0; i < v.size(); i++) {
		tmp[i] = c * v[i];
	}
	return tmp;
}
vector<double> operator + (const vector<double> v1, const vector<double> v2) {
	vector<double> tmp;
	tmp.resize(v1.size());
	for (int i = 0; i < v1.size(); i++) {
		tmp[i] = v1[i] + v2[i];
	}
	return tmp;
}
vector<double> operator - (const vector<double> v1, const vector<double> v2) {
	vector<double> tmp;
	tmp.resize(v1.size());
	for (int i = 0; i < v1.size(); i++) {
		tmp[i] = v1[i] - v2[i];
	}
	return tmp;
}
vector<vector<double>> InvJacobian(int key, const double t, const vector<double> x, const double tau) {
	double det;
	switch (key)
	{
	case 1:
		det = 1 + omega * omega;
		return { {-1 / det, -tau / det}, {tau * omega * omega / det, -1 / det} };
	case 2:
		det = 1 - 2 * tau + 2 * x[1] * tau - 16 * x[1] * tau * tau;
		return { {(-1 - 2 * x[1] * tau) / det, -2 * x[1] * tau / det }, {-6 * tau / det, (-1 + 2 * tau) / det } };
	case 3:
		det = 1 + 2 * tau * x[0] + 4 * x[1] * tau * tau;
		return{ {-1 / det, 2 * tau * x[1] / det}, {-2 * tau / det, (-1 - 2 * x[0] * tau) / det} };
	case 4:
		det = tau / (-b * sigma + b * r * sigma - x[0] * x[0] * sigma - x[0] * x[1] * sigma);
		return { {-1 + det * (b + x[0] * x[0]), det * b * sigma, det * x[0] * sigma}, {det * (b * r - x[0] * x[1]), -1 + det * b * sigma, det * x[0] * sigma},
			{det * (r * x[0] + x[1]), det * sigma * (x[0] + x[1]), -1 + det * sigma * (1 - r)} };
	case 5:
		det = 1 - 1.6 * tau + 2 * tau * tau * x[0] + 1.6 * tau * x[0] * x[0] + 3.2 * tau * tau * x[0] * x[1];
		return{ {(-1 + 1.6 * tau * (1 - x[0] * x[0])) / det, tau / det}, {(2 * tau * x[0] + 3.2 * x[0] * x[1]) / det, -1 / det} };
	}
}
vector<vector<double>> InvJacobian_symcirc(int key, const double t, const vector<double> x, const double tau) {
	double det;
	switch (key)
	{
	case 1:
		det = (4 + omega * omega);
		return { {-2 / det, -tau / det}, {tau * omega * omega / det, -2 / det} };
	case 2:
		det = 4 - 4 * tau + 4 * x[1] * tau - 16 * x[1] * tau * tau;
		return { {(-2 - 2 * x[1] * tau) / det, -2 * x[1] * tau / det }, {-6 * tau / det, (-2 + 2 * tau) / det } };
	case 3:
		det = 4 + 4 * x[0] * tau + 4 * x[1] * tau * tau;
		return{ {-2 / det, 2 * tau * x[1] / det}, {-2 * tau / det, (-2 - 2 * x[0] * tau) / det} };
	case 4:
		det = (-1 - b * tau) * (1 + tau + sigma * tau + sigma * tau * tau - r * sigma * tau * tau) + x[0];
		return { {-2 + det * (b + x[0] * x[0]), det * b * sigma, det * x[0] * sigma}, {det * (b * r - x[0] * x[1]), -2 + det * b * sigma, det * x[0] * sigma},
			{det * (r * x[0] + x[1]), det * sigma * (x[0] + x[1]), -2 + det * sigma * (1 - r)} };
	case 5:
		det = 4 - 3.2 * tau + 2 * tau * tau * x[0] + 1.6 * tau * x[0] * x[0] + 3.2 * tau * tau * x[0] * x[1];
		return{ {(-2 + 1.6 * tau * (1 - x[0] * x[0]))/ det, tau / det}, {(2 * tau * x[0] + 3.2 * x[0] * x[1]) / det, -2 / det} };
	}
}
double sqr(const vector<double> v1, const vector<double> v2) {
	double tmp = 0;
	for (int i = 0; i < v1.size(); i++) {
		tmp += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	return sqrt(tmp);
}
double maxnorm(const vector<double> v1, const vector<double> v2) {
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
vector<double> f(const double t, vector<double> x, int key) {
	switch (key)
	{
	case 1:
		return { x[1], -omega * omega * x[0] };
	case 2:
		//cout << 2 * x[0] + x[1] * x[1] - 1 << " " << 6 * x[0] - x[1] * x[1] + 1 << endl;
		return { 2 * x[0] + x[1] * x[1] - 1, 6 * x[0] - x[1] * x[1] + 1 };
	case 3:
		return { 1 - x[0] * x[0] - x[1] * x[1], 2 * x[0] };
	case 4:
		return { sigma * (x[1] - x[0]), x[0] * (r - x[2]) - x[1], x[0] * x[1] - b * x[2] };
	case 5:
		return { x[1], 5 * sin(0.01 * t) + 1.6 * x[1] * (1 - x[0] * x[0]) - x[0]*x[0]};
	}
}
vector<double> Gaussf(const double t, const vector<double> x, int key) {
	switch (key)
	{
	case 1:
		return { x[1], -omega * omega * x[0] };
	case 2:
		return { 1, 1, 1 };
	}
}
vector<double> newton_system(int key, const double t, const vector<double> xstart, const double eps, const double tau) {
	//double h = 1e-7;
	int m = xstart.size();
	std::vector<double> xtmp(m, 0), xcur(m, 0);
	xcur = xstart;
	for (int i = 0; i < m; i++) {
		xtmp[i] = xstart[i] - 0.1;
	}
	while (sqr(xcur, xtmp) > eps) {
		xtmp = xcur;
		xcur = xtmp - (InvJacobian(key, t, xtmp, tau) * (tau * f(t, xtmp, key) - xtmp + xstart));
	}
	/*for (int i = 0; i < xcur.size(); i++) {
		cout << xcur[i] << endl;
	}*/
	return xcur;
}
vector<double> newton_system_symcirc(int key, const double t, const vector<double> xstart, const double eps, const double tau, const vector<double> lastf) {
	//double h = 1e-7;
	int m = xstart.size();
	std::vector<double> xtmp(m, 0), xcur(m, 0);
	xcur = xstart;
	for (int i = 0; i < m; i++) {
		xtmp[i] = xstart[i] - 0.1;
	}
	while (sqr(xcur, xtmp) > eps) {
		xtmp = xcur;
		xcur = xtmp - (InvJacobian_symcirc(key, t, xtmp, tau) * (tau * f(t, xtmp, key) - 2 * xtmp + 2 * xstart + tau * lastf));
		/*cout << endl;
		for (int i = 0; i < xcur.size(); i++) {
			cout << xcur[i] << endl;
		}*/
	}
	return xcur;
}
class Cauchy_problem
{
	int n; //порядок ДУ.
	int key; //параметр, отвечающий за выбор правой части
	vector<double> init_cond;
public:
	void initialization(const vector<double> tmp, const int m, int tkey) {
		n = m;
		key = tkey;
		init_cond.resize(n);
		if (tmp.size() != n) {
			cout << "Неверно введены начальные данные" << endl;
		}
		for (int i = 0; i < n; i++) {
			init_cond[i] = tmp[i];
		}
	}
	void print() {
		cout << "порядок уравнения " << n << endl;
		cout << "ключик на правую часть " << key << endl;
		cout << endl;
		cout << "начальные условия" << endl;
		for (int i = 0; i < n; i++) {
			cout << init_cond[i] << endl;
		}
	}
	double explicit_Gaussian_method(int M, const double t0, const double T, int div) { //тау - параметр, M - количество шагов.
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		//cout << key << endl;
		double t = t0;
		y[0] = init_cond;
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		for (int i = 1; i < M + 1; i++) {
			t = t + tau;
			y[i] = y[i - 1] + tau * f(t, y[i - 1], key);
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				//cout << normt << endl;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		/*if (key != 1) {
		}*/
		/*cout << "Gauss_exp" << endl;
		cout << norm << endl;*/
		return norm;
	}
	vector<vector<double>> explicit_Gaussian_method(const int M, const double t0, const double T, const bool flag = true) { //тау - параметр, M - количество шагов.
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		//cout << key << endl;
		double t = t0;
		y[0] = init_cond;
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_expGauss.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			for (int i = 1; i < M + 1; i++) {
				t = t + tau;
				y[i] = y[i - 1] + tau * f(t, y[i - 1], key);
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
			}
		}
		out.close();
		if (key != 1) {
			t = t0 + tau;
			for (int i = 1; i < M; i++) {
				normt = sqr((1 / (2*tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				t = t + tau;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		if (flag) {
			cout << "Gauss_exp" << endl;
			cout << norm << endl;
		}
		return y;
	}
	double nonexplicit_Gaussian_method(int M, const double t0, const double T, int div) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		for (int i = 1; i < M + 1; i++) {
			t = t + tau;
			y[i] = newton_system(key, t, y[i - 1], eps, tau);
			/*for (int k = 0; k < n; k++) {
				cout << y[i][k] << endl;
			}*/
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		return norm;
	}
	vector<vector<double>> nonexplicit_Gaussian_method(const int M, const double t0, const double T, const bool flag = true) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_nonexpGauss.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			for (int i = 1; i < M + 1; i++) {
				t = t + tau;
				y[i] = newton_system(key, t, y[i - 1], eps, tau);
				/*for (int k = 0; k < n; k++) {
					cout << y[i][k] << endl;
				}*/
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
			}
		}
		out.close();
		if (key != 1) {
			t = t0 + tau;
			for (int i = 1; i < M; i++) {
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				t = t + tau;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		if (flag) {
			cout << "Gauss_nonexp" << endl;
			cout << norm << endl;
		}
		return y;
	}
	double symmetrical_circuit(int M, const double t0, const double T, int div) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		vector<double> prevf;
		//ofstream out;       // поток для записи
		for (int i = 1; i < M + 1; i++) {
			prevf = f(t, y[i - 1], key);
			t = t + tau;
			y[i] = newton_system_symcirc(key, t, y[i - 1], eps, tau, prevf);
			/*for (int k = 0; k < n; k++) {
				cout << y[i][k] << endl;
			}*/
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		//out.close();
		return norm;
	}
	vector<vector<double>> symmetrical_circuit(const int M, const double t0, const double T, const bool flag = true) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		vector<double> prevf;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_symcirc.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			for (int i = 1; i < M + 1; i++) {
				prevf = f(t, y[i - 1], key);
				t = t + tau;
				y[i] = newton_system_symcirc(key, t, y[i - 1], eps, tau, prevf);
				/*for (int k = 0; k < n; k++) {
					cout << y[i][k] << endl;
				}*/
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
			}
		}
		out.close();
		if (key != 1) {
			t = t0 + tau;
			for (int i = 1; i < M; i++) {
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				//cout << normt;
				t = t + tau;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		if (flag) {
			cout << "sym_circ" << endl;
			cout << norm << endl;
		}
		return y;
	}
	double RungeKutta2(int M, const double t0, const double T, int div) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		vector<double> k1, k2;
		ofstream out;       // поток для записи
		/*out.open("D:\\vichy_umf_output_RK2.txt");
		if (out.is_open()) {*/
			/*out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}*/
			//out << endl;
		for (int i = 1; i < M + 1; i++) {
			k1 = f(t, y[i - 1], key);
			k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
			t = t + tau;
			y[i] = y[i - 1] + tau * k2;
			/*for (int k = 0; k < n; k++) {
				cout << y[i][k] << endl;
			}*/
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
			}
			/*out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[i][j] << " ";
			}
			out << endl;*/
		}
		/*out.close();*/
		return norm;
	}
	vector<vector<double>> RungeKutta2(const int M, const double t0, const double T, const bool flag = true) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		vector<double> k1, k2;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_RK2.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			for (int i = 1; i < M + 1; i++) {
				k1 = f(t, y[i - 1], key);
				k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
				t = t + tau;
				y[i] = y[i - 1] + tau * k2;
				/*for (int k = 0; k < n; k++) {
					cout << y[i][k] << endl;
				}*/
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
			}
		}
		out.close();
		if (key != 1) {
			t = t0 + tau;
			for (int i = 1; i < M; i++) {
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				t = t + tau;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		if (flag) {
			cout << "RK2" << endl;
			cout << norm << endl;
		}
		return y;
	}
	double RungeKutta4(int M, const double t0, const double T, int div) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		vector<double> k1, k2, k3, k4;
		for (int i = 1; i < M + 1; i++) {
			k1 = f(t, y[i - 1], key);
			k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
			k3 = f(t + tau / 2, y[i - 1] + (tau / 2) * k2, key);
			k4 = f(t + tau, y[i - 1] + tau * k3, key);
			t = t + tau;
			y[i] = y[i - 1] + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
			}
			/*for (int j = 0; j < n; j++) {
				cout << y[i][j] << endl;
			}
			cout << endl;*/
			 
		}
		return norm;
	}
	vector<vector<double>> RungeKutta4(const int M, const double t0, const double T, const bool flag = true) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		vector<double> k1, k2, k3, k4;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_RK4.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			for (int i = 1; i < M + 1; i++) {
				k1 = f(t, y[i - 1], key);
				k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
				k3 = f(t + tau / 2, y[i - 1] + (tau / 2) * k2, key);
				k4 = f(t + tau, y[i - 1] + tau * k3, key);
				t = t + tau;
				y[i] = y[i - 1] + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
				/*for (int j = 0; j < n; j++) {
					cout << y[i][j] << endl;
				}
				cout << endl;*/
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
			}
		}
		out.close();
		if (key != 1) {
			t = t0 + tau;
			for (int i = 1; i < M; i++) {
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				t = t + tau;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		if (flag) {
			cout << "RK4" << endl;
			cout << norm << endl;
		}
		return y;
	}
	double Adams_method(const int M, const double t0, const double T, int div) {
		//сначала получаем три первых значения функции с помощью метода Рунге-Кутты 4 порядка
		//потом используем метод Адамса с прогнозом-коррекции.
		vector<vector<double>> ff;
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		ff.reserve(M);
		for (int i = 0; i < M; i++) {
			ff.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		vector<double> k1, k2, k3, k4;
		double norm = 0;
		double normt;
		ff[0] = f(t, y[0], key);
		for (int i = 1; i < 4; i++) {
			k1 = f(t, y[i - 1], key);
			k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
			k3 = f(t + tau / 2, y[i - 1] + (tau / 2) * k2, key);
			k4 = f(t + tau, y[i - 1] + tau * k3, key);
			y[i] = y[i - 1] + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
			t = t + tau;
			ff[i] = f(t, y[i], key);
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		/*for (int j = 0; j < 3; j++) {
			cout << "AAA" << endl;
			for (int i = 0; i < ff[j].size(); i++) {
				cout << ff[j][i] << endl;
			}
		}*/
		for (int i = 4; i < M + 1; i++) {
			t = t + tau;
			y[i] = y[i - 1] + (tau / 24) * (55 * ff[3] - 59 * ff[2] + 37 * ff[1] - 9 * ff[0]);
			ff[0] = ff[1];
			ff[1] = ff[2];
			ff[2] = ff[3];
			ff[3] = f(t, y[i], key);
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		return norm;
	}
	vector<vector<double>> Adams_method(const int M, const double t0, const double T, const bool flag = true) {
		//сначала получаем три первых значения функции с помощью метода Рунге-Кутты 4 порядка
		//потом используем метод Адамса с прогнозом-коррекции.
		vector<vector<double>> ff;
		vector<vector<double>> y; 
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		ff.reserve(M);
		for (int i = 0; i < M; i++) {
			ff.emplace_back(vector<double>(n, 0));
		}
		/*for (int j = 0; j < 4; j++) {
			cout << "AAA" << endl;
			for (int i = 0; i < ff[j].size(); i++) {
				cout << ff[j][i] << endl;
			}
		}*/
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double norm = 0;
		double normt;
		vector<double> k1, k2, k3, k4;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_Adams.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			ff[0] = f(t, y[0], key);
			for (int i = 1; i < 4; i++) {
				k1 = f(t, y[i - 1], key);
				k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
				k3 = f(t + tau / 2, y[i - 1] + (tau / 2) * k2, key);
				k4 = f(t + tau, y[i - 1] + tau * k3, key);
				y[i] = y[i - 1] + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
				t = t + tau;
				ff[i] = f(t, y[i], key);
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
			}
			/*for (int j = 0; j < 3; j++) {
				cout << "AAA" << endl;
				for (int i = 0; i < ff[j].size(); i++) {
					cout << ff[j][i] << endl;
				}
			}*/
			for (int i = 4; i < M + 1; i++) {
				y[i] = y[i - 1] + (tau / 24) * (55 * ff[3] - 59 * ff[2] + 37 * ff[1] - 9 * ff[0]);
				ff[0] = ff[1];
				ff[1] = ff[2];
				ff[2] = ff[3];
				ff[3] = f(t, y[i], key);
				t = t + tau;
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
			}
		}
		out.close();
		if (key != 1) {
			//t = t0 + tau;
			for (int i = 1; i < M; i++) {
				t = t0 + tau * (i);
				//cout << "y " << y[i][0] << " " << y[i][1] << endl;
				//cout << "f " << f(t, y[i], key)[0] << " " << f(t, y[i], key)[1] << endl;
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				//cout << normt << endl;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		if (flag) {
			cout << "Adams" << endl;
			cout << norm << endl;
		}
		return y;
	}
	vector<vector<double>> Adams_method_modification(const int M, const double t0, const double T, const bool flag = true) {
		//сначала получаем три первых значения функции с помощью метода Рунге-Кутты 4 порядка
		//потом используем метод Адамса с прогнозом-коррекции.
		vector<vector<double>> ff;
		vector<vector<double>> y;
		vector<vector<double>> ytmp;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		int p = 3;
		ytmp.reserve(p+1);
		for (int i = 0; i < p+1; i++) {
			ytmp.emplace_back(vector<double>(n, 0));
		}
		ff.reserve(M);
		for (int i = 0; i < M; i++) {
			ff.emplace_back(vector<double>(n, 0));
		}
		/*for (int j = 0; j < 4; j++) {
			cout << "AAA" << endl;
			for (int i = 0; i < ff[j].size(); i++) {
				cout << ff[j][i] << endl;
			}
		}*/
		double t = t0;
		y[0] = init_cond;
		ff[0] = f(t, y[0], key);
		ytmp[0] = y[0];
		double tau = (T - t0) / M;
		double tauhelp = pow(tau, 5);
		double norm = 0;
		double normt;
		vector<double> k1, k2, k3, k4;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_Adams.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			for (int i = 1; i < 4; i++) {
				ytmp[i + 1] = ytmp[i - 1] + tauhelp * f(t, y[i - 1], key);
				t = tau + tauhelp;

				//cout << "AAA " << i << endl;
				/*k1 = f(t, ytmp[i - 1], key);
				k2 = f(t + tauhelp / 2, ytmp[i - 1] + (tauhelp / 2) * k1, key);
				k3 = f(t + tauhelp / 2, y[i - 1] + (tauhelp / 2) * k2, key);
				k4 = f(t + tauhelp, ytmp[i - 1] + tauhelp * k3, key);
				ytmp[i] = ytmp[i - 1] + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
				t = t + tauhelp;
				ff[i] = f(t, y[i], key);
				//cout << "CCC" << endl;
				if (i % (p/3) == 0) {
					int k = i / (p/3);
					y[k] = ytmp[i];
					ff[k] = f(t, y[k], key);
					out << t << " ";
					for (int j = 0; j < n; j++) {
						out << y[k][j] << " ";
					}
					out << endl;
				}*/
				//cout << "DDD" << endl;
				/*for (int j = 0; j < 3; j++) {
					cout << "AAA" << endl;
					for (int i = 0; i < ff[j].size(); i++) {
					cout << ff[j][i] << endl;
				}*/
			}
			for (int i = 0; i < 4; i++) {
				cout << "y " << y[i][0] << " " << y[i][1] << endl;
			}
			//cout << "BBB" << endl;
			for (int i = 4; i < M + 1; i++) {
				y[i] = y[i - 1] + (tau / 24) * (55 * ff[3] - 59 * ff[2] + 37 * ff[1] - 9 * ff[0]);
				ff[0] = ff[1];
				ff[1] = ff[2];
				ff[2] = ff[3];
				ff[3] = f(t, y[i], key);
				t = t + tau;
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
			}
		}
		out.close();
		if (key != 1) {
			//t = t0 + tau;
			for (int i = 1; i < M; i++) {
				t = t0 + tau * (i);
				//cout << "y " << y[i][0] << " " << y[i][1] << endl;
				//cout << "f " << f(t, y[i], key)[0] << " " << f(t, y[i], key)[1] << endl;
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				cout << t << " " << normt << endl;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		if (flag) {
			cout << "Adams" << endl;
			cout << norm << endl;
		}
		return y;
	}
	double correction_forecast(const int M, const double t0, const double T, int div) {
		vector<vector<double>> ff;
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		ff.reserve(M);
		for (int i = 0; i < M; i++) {
			ff.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double normt;
		double norm = 0;
		vector <double> ftmp;
		vector<double> k1, k2, k3, k4;
		ff[0] = f(t, y[0], key);
		for (int i = 1; i < 4; i++) {
			k1 = f(t, y[i - 1], key);
			k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
			k3 = f(t + tau / 2, y[i - 1] + (tau / 2) * k2, key);
			k4 = f(t + tau, y[i - 1] + tau * k3, key);
			y[i] = y[i - 1] + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
			t = t + tau;
			ff[i] = f(t, y[i], key);
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		for (int i = 4; i < M + 1; i++) {
			y[i] = y[i - 1] + (tau / 24) * (55 * ff[3] - 59 * ff[2] + 37 * ff[1] - 9 * ff[0]);
			ff[0] = ff[1];
			ff[1] = ff[2];
			ff[2] = ff[3];
			ftmp = f(t, y[i], key);
			y[i] = y[i - 1] + (tau / 24) * (9 * ftmp + 19 * ff[2] - 5 * ff[1] + ff[0]);
			t = t + tau;
			ff[3] = f(t, y[i], key);
			if (key == 1) {
				normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		return norm;
	}
	vector<vector<double>> correction_forecast(const int M, const double t0, const double T, const bool flag = true) {
		vector<vector<double>> ff;
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		ff.reserve(M);
		for (int i = 0; i < M; i++) {
			ff.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double eps = pow(10, -8);
		double tau = (T - t0) / M;
		double normt;
		double norm = 0;
		vector<double> ftmp;
		vector<double> y0;
		vector<double> k1, k2, k3, k4;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_corfor.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			ff[0] = f(t, y[0], key);
			for (int i = 1; i < 4; i++) {
				k1 = f(t, y[i - 1], key);
				k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
				k3 = f(t + tau / 2, y[i - 1] + (tau / 2) * k2, key);
				k4 = f(t + tau, y[i - 1] + tau * k3, key);
				y[i] = y[i - 1] + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
				t = t + tau;
				ff[i] = f(t, y[i], key);
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
			}
			for (int i = 4; i < M + 1; i++) {
				y0 = y[i - 1] + (tau / 24) * (55 * ff[3] - 59 * ff[2] + 37 * ff[1] - 9 * ff[0]);
				ff[0] = ff[1];
				ff[1] = ff[2];
				ff[2] = ff[3];
				t = t + tau;
				ftmp = f(t, y0, key);
				y[i] = y[i - 1] + (tau / 24) * (9 * ftmp + 19 * ff[2] - 5 * ff[1] + ff[0]);
				ff[3] = f(t, y[i], key);
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
				if (key == 1) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					if (normt > norm) {
						norm = normt;
					}
				}
			}
		}
		out.close();
		if (key != 1) {
			t = t0 + tau;
			for (int i = 1; i < M; i++) {
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				t = t + tau;
				if (normt > norm) {
					norm = normt;
				}
			}
		}
		if (flag) {
			cout << "correction_forecast" << endl;
			cout << norm << endl;
		}
		return y;
	}
	void RungeKutta2_automatic_step_hypotesic(const int M, const double t0, const double T) {
		vector<vector<double>> y;
		y.reserve(M);
		for (int i = 0; i < M + 1; i++) {
			y.emplace_back(vector<double>(n, 0));
		}
		double t = t0;
		y[0] = init_cond;
		double tau = (T - t0) / M;
		double eps = tau * tau;
		double norm = 0;
		double normt;
		double veps;
		vector<double> k1, k2, kk3;
		vector<double> y_corrector;
		ofstream out;       // поток для записи
		out.open("D:\\vichy_umf_output_RK2.txt");
		if (out.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y[0][j] << " ";
			}
			out << endl;
			for (int i = 1; i < M + 1; i++) {
				k1 = f(t, y[i - 1], key);
				k2 = f(t + tau / 2, y[i - 1] + (tau / 2) * k1, key);
				kk3 = f(t + tau, y[i - 1] - k1 + 2*k2, key);
				y[i] = y[i - 1] + tau * k2;
				y_corrector = (tau / 6) * (k1 + 4 * k2 + kk3);
				veps = sqr(y[i], y_corrector);
				if (veps >= eps) {
					tau = tau * pow(eps / veps, 1/3);
				}
				t = t + tau;
				out << t << " ";
				for (int j = 0; j < n; j++) {
					out << y[i][j] << " ";
				}
				out << endl;
			}
		}
		out.close();
		cout << "RK2" << endl;
		cout << norm << endl;
	}
	void RungeKutta2_automatic_step_Runge_rule(const int M, const double t0, const double T) {
		vector<double> y0, y1, y2, ytmp, ys;
		double t = t0;
		bool flag;
		y0 = init_cond;
		double tau = (T - t0) / (2*M);
		double eps = tau * tau;
		double norm = 0;
		double normt, veps;
		vector<double> k1, k2, kk1, kk2;
		vector<vector<double>> step;
		step.push_back({ t, tau });
		ofstream out, file1, file2;       // поток для записи
		out.open("D:\\vichy_umf_output_RK2_Runge_rule_graph.txt");
		file1.open("D:\\vichy_umf_output_RK2_Runge_rule_graph_error");
		if (out.is_open() && file1.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y0[j] << " ";
			}
			out << 0;
			out << endl;
			while(t < T){
				flag = true;
				while (flag) {
					k1 = f(t, y0, key);
					k2 = f(t + tau / 2, y0 + (tau / 2) * k1, key);
					y1 = y0 + tau * k2;
					k1 = f(t + tau, y1, key);
					k2 = f(t + 1.5 * tau, y1 + (tau / 2) * k1, key);
					y2 = y1 + tau * k2;
					kk1 = f(t, y0, key);
					kk2 = f(t + 2 * tau, y0 + tau * kk1, key);
					ytmp = y0 + 2*tau * kk2;
					veps = sqr(ytmp, y2);
					if (veps > eps) {
						tau = tau / 2;
						//cout << tau << endl;
					}
					if (veps <= eps) {
						flag = false;
						//cout << tau << endl;
						tau = tau * 2;
						step.push_back({ t, tau });
						ys = y2;
					}
				}
				y0 = ys;
				t = t + tau;
				out << t << " ";
				//file1 << t << " ";
				for (int j = 0; j < n; j++) {
					out << ys[j] << " ";
				}
				normt = sqr(ys, { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
				out << normt;
				out << endl;
			}
		}
		out.close();
		file1.close();
		/*file1.open("D:\\vichy_umf_output_RK2_Runge_rule_graph_error.txt");
		if (file1.is_open()) {
			if (key == 1) {
				t = t0;
				for (int i = 0; i < y.size(); i++) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					file1 << t << " " << normt << endl;
					t = t + tau;
				}
			}
		}
		file1.close();*/
		file2.open("D:\\vichy_umf_output_RK2_Runge_rule_graph_step.txt");
		if (file2.is_open()) {
			for (int i = 0; i < step.size(); i++) {
				file2 << step[i][0] << " " << step[i][1] << endl;
			}
		}
		file2.close();
		/*if (key != 1) {
			t = t0 + tau;
			for (int i = 1; i < M; i++) {
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				t = t + tau;
				if (normt > norm) {
					norm = normt;
				}
			}
		}*/
		cout << "RK2_Runge_rule" << endl;
		cout << norm << endl;
	}
	void RungeKutta4_automatic_step_Runge_rule(const int M, const double t0, const double T) {
		vector<double> y0, y1, y2, ytmp, ys;
		double t = t0;
		y0 = init_cond;
		double tau = (T - t0) / (2 * M);
		double eps = tau * tau * tau * tau;
		//cout << eps << endl;
		double norm = 0;
		double normt, veps;
		bool flag;
		vector<vector<double>> step;
		vector<double> k1, k2, k3, k4, kk1, kk2, kk3, kk4;
		ofstream out, file1, file2;       // поток для записи
		out.open("D:\\vichy_umf_output_RK4_Runge_rule_graph.txt");
		file1.open("D:\\vichy_umf_output_RK4_Runge_rule_graph_error.txt");
		if (out.is_open() && file1.is_open()) {
			out << t << " ";
			for (int j = 0; j < n; j++) {
				out << y0[j] << " ";
			}
			out << 0;
			out << endl;
			while(t < T) {
				flag = true;
				while (flag) {
					k1 = f(t, y0, key);
					k2 = f(t + tau / 2, y0 + (tau / 2) * k1, key);
					k3 = f(t + tau / 2, y0 + (tau / 2) * k2, key);
					k4 = f(t + tau, y0 + tau * k3, key);
					y1 = y0 + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
					k1 = f(t, y0, key);
					k2 = f(t + 1.5 * tau, y1 + (tau / 2) * k1, key);
					k3 = f(t + 1.5 * tau, y1 + (tau / 2) * k2, key);
					k4 = f(t + 2 * tau, y1 + tau * k3, key);
					y2 = y1 + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
					kk1 = f(t, y0, key);
					kk2 = f(t + tau, y0 + tau * k1, key);
					kk3 = f(t + tau, y0 + tau * k2, key);
					kk4 = f(t + 2 * tau, y0 + 2 * tau * k3, key);
					ytmp = y0 + (tau / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
					veps = sqr(ytmp, y2);
					if (veps > eps) {
						tau = tau / 2;
						//cout << tau << endl;
					}
					if (veps <= eps) {
						flag = false;
						//cout << tau << endl;
						tau = tau * 2;
						step.push_back({ t, tau });
						ys = y2;
					}
				}
				y0 = ys;
				t = t + tau;
				out << t << " ";
				//file1 << t << " ";
				for (int j = 0; j < n; j++) {
					out << ys[j] << " ";
				}
				normt = sqr(ys, { cos(2 * t), -2 * sin(2 * t) });
				if (normt > norm) {
					norm = normt;
				}
				out << normt;
				out << endl;
			}
		}
		/*k1 = f(t, y[2*M - 1], key);
		k2 = f(t + tau / 2, y[2*M - 1] + (tau / 2) * k1, key);
		y[i] = y[i - 1] + tau * k2;*/
		out.close();
		file1.close();
		/*file1.open("D:\\vichy_umf_output_RK2_Runge_rule_graph_error.txt");
		if (file1.is_open()) {
			if (key == 1) {
				t = t0;
				for (int i = 0; i < y.size(); i++) {
					normt = sqr(y[i], { cos(2 * t), -2 * sin(2 * t) });
					file1 << t << " " << normt << endl;
					t = t + tau;
				}
			}
		}
		file1.close();*/
		file2.open("D:\\vichy_umf_output_RK4_Runge_rule_graph_step.txt");
		if (file2.is_open()) {
			for (int i = 0; i < step.size(); i++) {
				file2 << step[i][0] << " " << step[i][1] << endl;
			}
		}
		file2.close();
		/*if (key != 1) {
			t = t0 + tau;
			for (int i = 1; i < M; i++) {
				normt = sqr((1 / (2 * tau)) * (y[i + 1] - y[i - 1]), f(t, y[i], key));
				t = t + tau;
				if (normt > norm) {
					norm = normt;
				}
			}
		}*/
		cout << "RK4_Runge_rule" << endl;
		cout << norm << endl;
	}
	void RungeKutta4_faze_portrait(int M, const double t0, const double T) {
		double a = -2;
		double b = -2;
		double c = -2;
		vector<vector<double>> solve;
		double tau = (T - t0) / M;
		double t = t0;
		cout << "A" << endl;
		fstream out("D:\\vichy_umf_output_faze_portrait.txt");
		//if (out.is_open()) {
			for (int i = 0; i <= 4; i++) {
				for (int j = 0; j <= 4; j++) {
					for (int h = 0; h <= 4; h++) {
						t = t0;
						this->init_cond = { -a + i, -b + j, -c + h };
						solve = this->RungeKutta4(M, t0, T, false);
						cout << "AAA" << endl;
						for (int k = 0; k < solve.size(); k++) {
							out << t << " ";
							for (int p = 0; p < n; p++) {
								out << solve[k][p] << " ";
							}
							out << endl;
							t = t + tau;
						}
					}
				}
			}
		//}
		out.close();
	}
	int approximation_estimate_abserror(int M, const double t0, const double T, const int tkey) {
		int div = 2;
		double norm1, norm2;
		switch (tkey) {
		case 1:
			norm1 = this->explicit_Gaussian_method(M, t0, T, div);
			cout << "Explicit_Gauss_methond" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << "---" << endl;
			for (int i = 0; i < 5; i++) {
				M = M * div;
				norm2 = this->explicit_Gaussian_method(M, t0, T, div);
				//cout << norm2 << endl;
				cout << (T - t0) / M  << " " << norm2 << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(2) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				norm1 = norm2;
			};
			return 0;
		case 2:
			norm1 = this->nonexplicit_Gaussian_method(M, t0, T, div);
			cout << "Nonexplicit_Gauss_methond" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << "---" << endl;
			for (int i = 0; i < 5; i++) {
				M = M * div;
				norm2 = this->nonexplicit_Gaussian_method(M, t0, T, div);
				//cout << norm2 << endl;
				cout << (T - t0) / M << " " << norm2 << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(2) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				norm1 = norm2;
			};
			return 0;
		case 3:
			norm1 = this->symmetrical_circuit(M, t0, T, div);
			//cout << norm1 << endl;norm2;
			cout << "Symmetrical_circuit" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << "---" << endl;
			for (int i = 0; i < 5; i++) {
				M = M * div;
				norm2 = this->symmetrical_circuit(M, t0, T, div);
				//cout << norm2 << endl;
				cout << (T - t0) / M << " " << norm2 << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(2) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				norm1 = norm2;
			};
			return 0;
		case 4:
			norm1 = this->RungeKutta2(M, t0, T, div);
			//cout << norm1 << endl;
			double norm2;
			cout << "RungeKutta2" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << "---" << endl;
			for (int i = 0; i < 5; i++) {
				M = M * div;
				norm2 = this->RungeKutta2(M, t0, T, div);
				//cout << norm2 << endl;
				cout << (T - t0) / M << " " << norm2 << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(2) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				norm1 = norm2;
			};
			return 0;
		case 5:
			norm1 = this->RungeKutta4(M, t0, T, div);
			//cout << norm1 << endl;
			cout << "RungeKutta4" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << "---" << endl;
			for (int i = 0; i < 5; i++) {
				M = M * div;
				norm2 = this->RungeKutta4(M, t0, T, div);
				//cout << norm2 << endl;
				cout << (T - t0) / M << " " << norm2 << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(2) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				norm1 = norm2;
			};
			return 0;
		case 6:
			norm1 = this->Adams_method(M, t0, T, div);
			//cout << norm1 << endl;
			cout << "Adams_method" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << "---" << endl;
			for (int i = 0; i < 5; i++) {
				M = M * div;
				norm2 = this->Adams_method(M, t0, T, div);
				//cout << norm2 << endl;
				cout << (T - t0) / M << " " << norm2 << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(2) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				norm1 = norm2;
			};
			return 0;
		case 7:
			norm1 = this->correction_forecast(M, t0, T, div);
			//cout << norm1 << endl;
			cout << "correction_forecast" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << "---" << endl;
			for (int i = 0; i < 5; i++) {
				M = M * div;
				norm2 = this->correction_forecast(M, t0, T, div);
				//cout << norm2 << endl;
				cout << (T - t0) / M << " " << norm2 << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(2) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				norm1 = norm2;
			};
			return 0;
		}
	}
	int approximation_estimate_Eitkenrule(int M, const double t0, const double T, const int ttkey) {
		int div = 2;
		vector <double> vec1, vec2, vec3;
		double norm1, norm2;
		switch (ttkey) {
		case 1:
			vec1 = this->explicit_Gaussian_method(M, t0, T, false)[M - 1];
			M = M * div;
			vec2 = this->explicit_Gaussian_method(M, t0, T, false)[M - 1];
			norm1 = sqr(vec1, vec2);
			cout << "Explicit_Gauss_methond" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << endl;
			for (int i = 0; i < 6; i++) {
				M = M * div;
				vec1 = vec2;
				vec3 = this->explicit_Gaussian_method(M, t0, T, false)[M - 1];
				norm2 = sqr(vec2, vec3);
				cout << (T - t0) / M << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(div) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				vec2 = vec3;
				norm1 = norm2;
			}
			return 0;
		case 2:
			vec1 = this->nonexplicit_Gaussian_method(M, t0, T, false)[M - 1];
			M = M * div;
			vec2 = this->nonexplicit_Gaussian_method(M, t0, T, false)[M - 1];
			norm1 = sqr(vec1, vec2);
			cout << "Nonexplicit_Gauss_methond" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << endl;
			for (int i = 0; i < 6; i++) {
				vec1 = vec2;
				M = M * div;
				vec3 = this->nonexplicit_Gaussian_method(M, t0, T, false)[M - 1];
				norm2 = sqr(vec2, vec3);
				cout << (T - t0) / M << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(div) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				vec2 = vec3;
				norm1 = norm2;
			}
			return 0;
		case 3:
			vec1 = this->symmetrical_circuit(M, t0, T, false)[M];
			M = M * div;
			vec2 = this->symmetrical_circuit(M, t0, T, false)[M];
			norm1 = sqr(vec1, vec2);
			cout << "Symmetrical circuit" << endl;
			cout << (T - t0) / M << " " << norm1 << " --- " << endl;
			for (int i = 0; i < 6; i++) {
				vec1 = vec2;
				M = M * div;
				vec3 = this->symmetrical_circuit(M, t0, T, false)[M];
				norm2 = sqr(vec3, vec2);
				cout << (T - t0) / M << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(div) << endl;
				//cout << "i+1" << norm1 / norm2 << endl;
				vec2 = vec3;
				norm1 = norm2;
			}
			return 0;
			case 4:
				vec1 = this->RungeKutta2(M, t0, T, false)[M];
				M = M * div;
				vec2 = this->RungeKutta2(M, t0, T, false)[M];
				norm1 = sqr(vec1, vec2);
				cout << "Runge Kutta 2" << endl;
				cout << (T - t0) / M << " " << norm1 << " --- " << endl;
				for (int i = 0; i < 6; i++) {
					vec1 = vec2;
					M = M * div;
					vec3 = this->RungeKutta2(M, t0, T, false)[M];
					norm2 = sqr(vec3, vec2);
					cout << (T - t0) / M << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(div) << endl;
					//cout << "i+1" << norm1 / norm2 << endl;
					vec2 = vec3;
					norm1 = norm2;
				}
				return 0;
			case 5:
				vec1 = this->RungeKutta4(M, t0, T, false)[M];
				M = M * div;
				vec2 = this->RungeKutta4(M, t0, T, false)[M];
				norm1 = sqr(vec1, vec2);
				cout << "Runge Kutta 4" << endl;
				cout << (T - t0) / M << " " << norm1 << " --- " << endl;
				for (int i = 0; i < 6; i++) {
					vec1 = vec2;
					M = M * div;
					vec3 = this->RungeKutta4(M, t0, T, false)[M];
					norm2 = sqr(vec3, vec2);
					cout << (T - t0) / M << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(div) << endl;
					//cout << "i+1" << norm1 / norm2 << endl;
					vec2 = vec3;
					norm1 = norm2;
				}
				return 0;
			case 6:
				vec1 = this->Adams_method(M, t0, T, false)[M];
				M = M * div;
				vec2 = this->Adams_method(M, t0, T, false)[M];
				norm1 = sqr(vec1, vec2);
				cout << "Adams method" << endl;
				cout << (T - t0) / M << " " << norm1 << " --- " << endl;
				for (int i = 0; i < 6; i++) {
					vec1 = vec2;
					M = M * div;
					vec3 = this->Adams_method(M, t0, T, false)[M];
					norm2 = sqr(vec3, vec2);
					cout << (T - t0) / M << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(div) << endl;
					//cout << "i+1" << norm1 / norm2 << endl;
					vec2 = vec3;
					norm1 = norm2;
				}
				return 0;
			case 7:
				vec1 = this->correction_forecast(M, t0, T, false)[M];
				M = M * div;
				vec2 = this->correction_forecast(M, t0, T, false)[M];
				norm1 = sqr(vec1, vec2);
				cout << "Correction forecast" << endl;
				cout << (T - t0) / M << " " << norm1 << " --- " << endl;
				for (int i = 0; i < 6; i++) {
					vec1 = vec2;
					M = M * div;
					vec3 = this->correction_forecast(M, t0, T, false)[M];
					norm2 = sqr(vec3, vec2);
					cout << (T - t0) / M << " " << norm1 / norm2 << " " << log(norm1 / norm2) / log(div) << endl;
					//cout << "i+1" << norm1 / norm2 << endl;
					vec2 = vec3;
					norm1 = norm2;
				}
				return 0;
		}

	}
};
void approximation_estimate1(int M, const double t0, const double T, Cauchy_problem zad) {
	int div = 2;
	double norm1 = zad.explicit_Gaussian_method(M, t0, T, div);
	//cout << norm1 << endl;
	double norm2;
	cout << "Explicit_Gauss_methond" << endl;
	for (int i = 0; i < 5; i++) {
		M = M * div;
		norm2 = zad.explicit_Gaussian_method(M, t0, T, div);
		//cout << norm2 << endl;
		cout << log(norm1 / norm2) / log(2) << endl;
		//cout << "i+1" << norm1 / norm2 << endl;
		norm1 = norm2;
	}
}
void approximation_estimate2(int M, const double t0, const double T, Cauchy_problem zad) {
	int div = 2;
	double norm1 = zad.nonexplicit_Gaussian_method(M, t0, T, div);
	//cout << norm1 << endl;
	double norm2;
	cout << "Nonexplicit_Gauss_methond" << endl;
	for (int i = 0; i < 5; i++) {
		M = M * div;
		norm2 = zad.nonexplicit_Gaussian_method(M, t0, T, div);
		//cout << norm2 << endl;
		cout << log(norm1 / norm2) / log(2) << endl;
		//cout << "i+1" << norm1 / norm2 << endl;
		norm1 = norm2;
	}
}
void approximation_estimate3(int M, const double t0, const double T, Cauchy_problem zad) {
	int div = 2;
	double norm1 = zad.symmetrical_circuit(M, t0, T, div);
	//cout << norm1 << endl;
	double norm2;
	cout << "Symmetrical_circuit" << endl;
	for (int i = 0; i < 5; i++) {
		M = M * div;
		norm2 = zad.symmetrical_circuit(M, t0, T, div);
		//cout << norm2 << endl;
		cout << log(norm1 / norm2) / log(2) << endl;
		//cout << "i+1" << norm1 / norm2 << endl;
		norm1 = norm2;
	}
}
void approximation_estimate4(int M, const double t0, const double T, Cauchy_problem zad) {
	int div = 2;
	double norm1 = zad.RungeKutta2(M, t0, T, div);
	//cout << norm1 << endl;
	double norm2;
	cout << "RungeKutta2" << endl;
	for (int i = 0; i < 5; i++) {
		M = M * div;
		norm2 = zad.RungeKutta2(M, t0, T, div);
		//cout << norm2 << endl;
		cout << log(norm1 / norm2) / log(2) << endl;
		//cout << "i+1" << norm1 / norm2 << endl;
		norm1 = norm2;
	}
}
void approximation_estimate5(int M, const double t0, const double T, Cauchy_problem zad) {
	int div = 2;
	double norm1 = zad.RungeKutta4(M, t0, T, div);
	//cout << norm1 << endl;
	double norm2;
	cout << "RungeKutta4" << endl;
	for (int i = 0; i < 5; i++) {
		M = M * div;
		norm2 = zad.RungeKutta4(M, t0, T, div);
		//cout << norm2 << endl;
		cout << log(norm1 / norm2) / log(2) << endl;
		//cout << "i+1" << norm1 / norm2 << endl;
		norm1 = norm2;
	}
}
void approximation_estimate6(int M, const double t0, const double T, Cauchy_problem zad) {
	int div = 2;
	double norm1 = zad.Adams_method(M, t0, T, div);
	//cout << norm1 << endl;
	double norm2;
	cout << "Adams_method" << endl;
	for (int i = 0; i < 5; i++) {
		M = M * div;
		norm2 = zad.Adams_method(M, t0, T, div);
		//cout << norm2 << endl;
		cout << log(norm1 / norm2) / log(2) << endl;
		//cout << "i+1" << norm1 / norm2 << endl;
		norm1 = norm2;
	}
}
void approximation_estimate7(int M, const double t0, const double T, Cauchy_problem zad) {
	int div = 2;
	double norm1 = zad.correction_forecast(M, t0, T, div);
	//cout << norm1 << endl;
	double norm2;
	cout << "correction_forecast" << endl;
	for (int i = 0; i < 5; i++) {
		M = M * div;
		norm2 = zad.correction_forecast(M, t0, T, div);
		//cout << norm2 << endl;
		cout << log(norm1 / norm2) / log(2) << endl;
		//cout << "i+1" << norm1 / norm2 << endl;
		norm1 = norm2;
	}
} 
int main()
{
	setlocale(LC_ALL, "Russian");
	Cauchy_problem dequation;
	//key = 2 (0, 1) - седло (0, -1) - фокус.
	//key = 3 (0, 1) - центр (0, -1) - седло.
	//const vector<double> start = { 1, 0 };
	//const int n = start.size();
	//int key = 1; // задаем задачу Коши
	//dequation.initialization(start, n, key);
	//20, 0, 10, 10, 0, 10
	/*dequation.RungeKutta2(200, 0, 10);
	dequation.RungeKutta4(100, 0, 10);
	dequation.RungeKutta2_automatic_step_Runge_rule(60, 0, 10);
	dequation.RungeKutta4_automatic_step_Runge_rule(30, 0, 10);*/
	//dequation.RungeKutta2(10, 0, 10);
	//dequation.RungeKutta2_automatic_step_Runge_rule(5, 0, 10);//_graphics(10, 0, 2);
	//cout << "Absolute Error" << endl;
	//for (int i = 1; i < 8; i++) {
	//	dequation.approximation_estimate_abserror(10, 0, 1, i);
	//};
	//cout << endl;
	//const vector<double> start1 = { 0 , 0 };
	//const int n1 = start1.size();
	//int key1 = 2; // задаем задачу Коши
	//Cauchy_problem dequation1;
	//dequation1.initialization(start1, n1, key1);
	//cout << "Eitken Rule" << endl;
	//for (int i = 1; i < 8; i++) {
	//	dequation1.approximation_estimate_Eitkenrule(5, 0, 0.5, i);
	//};
	//cout << endl;
	//approximation_estimate1(8, 0, 1, dequation);
	//approximation_estimate2(8, 0, 1, dequation);
	//approximation_estimate3(8, 0, 1, dequation);
	//approximation_estimate4(8, 0, 1, dequation);
	//approximation_estimate5(8, 0, 1, dequation);
	//approximation_estimate6(8, 0, 1, dequation);
	//approximation_estimate7(8, 0, 1, dequation);
	//dequation.explicit_Gaussian_method(10000, 0, 1); //0.01 - 400
	//dequation.nonexplicit_Gaussian_method(10000, 0, 1); //0.01 - 400
	//dequation.symmetrical_circuit(12, 0, 1); //0.01 - 12, 
	//dequation.RungeKutta2(12, 0, 1);//0.01 - 12,
	//dequation.RungeKutta4(3, 0, 1);//0.01 - 3 
	//const vector<double> start = { 0, 0 };
	//const int n = start.size();
	//int key = 2; // задаем задачу Коши
	//dequation.initialization(start, n, key);
	//dequation.explicit_Gaussian_method(100, 0, 0.5);
	//dequation.nonexplicit_Gaussian_method(100, 0, 0.5);
	//dequation.symmetrical_circuit(100, 0, 0.5);
	//dequation.RungeKutta2(100, 0, 0.5);
	//dequation.RungeKutta4(100, 0, 0.5);
	//dequation.Adams_method(100, 0, 0.5);
	//dequation.correction_forecast(100, 0, 0.5);
	//dequation.tmp();
	//const vector<double> start = { 1, 1 };
	//const int n = start.size();
	//const int key = 3; // задаем задачу Коши
	//dequation.initialization(start, n, key);
	////dequation.print();
	//dequation.explicit_Gaussian_method(100, 0, 0.9);
	//dequation.nonexplicit_Gaussian_method(100, 0, 0.9);
	//dequation.symmetrical_circuit(100, 0, 0.9);
	//dequation.RungeKutta2(100, 0, 0.9);
	//dequation.RungeKutta4(100, 0, 0.5);
	//dequation.Adams_method(100, 0, 0.9);
	//dequation.correction_forecast(100, 0, 0.9);
	//const vector<double> start = { 1, 1 };
	//const int n = start.size();
	//int key = 3; // задаем задачу Коши
	//dequation.initialization(start, n, key);
	//dequation.RungeKutta2_automatic_step_Runge_rule(100, 0, 1);
	//dequation.RungeKutta4_automatic_step_Runge_rule(100, 0, 1);
	//dequation.print();
	/*dequation.explicit_Gaussian_method(100, 0, 1);
	dequation.nonexplicit_Gaussian_method(100, 0, 1);
	dequation.symmetrical_circuit(100, 0, 1);
	dequation.RungeKutta2(100, 0, 1);
	dequation.RungeKutta4(100, 0, 1);
	dequation.Adams_method(100, 0, 1);
	dequation.correction_forecast(100, 0, 1);*/
	////dequation.tmp();
	const vector<double> start = { 1, 1, 1 };
	const int n = start.size();
	int key = 4; // задаем задачу Коши
	dequation.initialization(start, n, key);
	dequation.RungeKutta4_faze_portrait(360, 0, 12);
	////dequation.print();
	//dequation.explicit_Gaussian_method(1000, 0, 1);
	//dequation.nonexplicit_Gaussian_method(1000, 0, 1);
	//dequation.symmetrical_circuit(1000, 0, 1);
	//dequation.RungeKutta2(1000, 0, 1);
	//dequation.RungeKutta4(1000, 0, 1);
	//dequation.Adams_method(1000, 0, 1);
	//dequation.correction_forecast(1000, 0, 1);
	//dequation.tmp();
	//const vector<double> start = { -2, 3 };
	//const int n = start.size();
	//int key = 5; // задаем задачу Коши
	//dequation.initialization(start, n, key);
	////dequation.print();
	////int B = 100000;
	////cout << 100.0 / B << endl;
	//dequation.explicit_Gaussian_method(10000, 0, 100);
	//dequation.nonexplicit_Gaussian_method(10000, 0, 100);
	//dequation.symmetrical_circuit(10000, 0, 100);
	//dequation.RungeKutta2(10000, 0, 100);
	//dequation.RungeKutta4(10000, 0, 100);
	//dequation.Adams_method(10000, 0, 100);
	//dequation.correction_forecast(10000, 0, 100);
	//const vector<double> start = { 0,0 };
	//const int n = start.size();
	//int key = 5; // задаем задачу Коши
	//dequation.initialization(start, n, key);
	//dequation.RungeKutta4_faze_portrait(1000, 0, 5);
}
