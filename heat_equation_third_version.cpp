// heat_equation.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include <iostream>
#include <vector>
#include <cmath>
#include<vector>
#include <fstream>

const double PI = 3.14159265358979323846;
const double Q = 10;
const double c = 0.5;
const double rho = 4;
const double t0 = 0.5;
const double alpha = 0.1;
const double betha = 1;
const double gamma = 3;
const double u0 = 0.3;
const double uL = 0.8;
const double L = 1;
const double x1 = 1 / 3;
const double x2 = 3 / 4;
const double k1 = 0.5;
const double k2 = 2.5;

std::vector<double> operator * (const std::vector<std::vector<double>>& A, const std::vector<double>& v) {

	std::vector<double> tmp;
	tmp.resize(v.size());

	for (int i = 0; i < v.size(); i++) {

		tmp[i] = 0;

		for (int j = 0; j < v.size(); j++) {
			tmp[i] += A[i][j] * v[j];
		}

	}

	return tmp;
}

std::vector<double> operator * (const double c, const std::vector<double>& v) {

	std::vector<double> tmp;

	tmp.resize(v.size());

	for (int i = 0; i < v.size(); i++) {
		tmp[i] = c * v[i];
	}

	return tmp;
}

std::vector<double> operator + (const std::vector<double>& v1, const std::vector<double>& v2) {

	std::vector<double> tmp;

	tmp.resize(v1.size());

	for (int i = 0; i < v1.size(); i++) {
		tmp[i] = v1[i] + v2[i];
	}

	return tmp;
}

std::vector<double> operator - (const std::vector<double>& v1, const std::vector<double>& v2) {

	std::vector<double> tmp;
	
	tmp.resize(v1.size());

	for (int i = 0; i < v1.size(); i++) {
		tmp[i] = v1[i] - v2[i];
	}

	return tmp;
}

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

double P(double t, const int key) {

	switch (key) {

	case 1:

		if (t < t0) {
			return Q;
		}
		else {
			return 0;
		}

	case 2:

		if (t < t0) {
			return 2 * Q * t;
		}

		else {
			return 0;
		}

	case 3:

		if (t < t0) {
			return 2 * Q * (t0 - t);
		}

		else {
			return 0;
		}

	case 4:
		if (t <= 0.5 * t0) {
			return 2 * Q * t;
		}

		else if ((t > 0.5 * t0) && (t < t0)) {
			return 2 * Q * (t0 - t);
		}

		else {
			return 0;
		}

	case 5:

		return 0;

	}
}

double u_init(double x, const int key, double t = 0) {

	switch (key) {

	case 1:

		return u0 + x * (L - x);

	case 2:

		return 0;

	case 3:

		return u0;

	case 4: //(uL - u0) * x + u0;// -(u0 / 2) * (x * x - L * x)

		if (x < 0) {
			return 1;
		}

		else if ((x >= 0) && (x <= 0.5)) { 
			return 1 - 2 * x / 3;
		}

		else if ((x >= 0.5) && (x < 1)) {
			return 4 * (1 - x) / 3;
		}

		else {
			return 0;
		}

	case 5:

		return (uL - u0) * x + u0;

	case 6:

		return exp(x) * sin(t);

	}

}
double exact_u(double x, double t) {
	//(uL - u0) * x + u0;// -(u0 / 2) * (x * x - L * x)
	/*if (x < 0) {
		return 1;
	}
	else if ((x >= 0) && (x <= 0.5)) {
		return 1 - 2 * x / 3;
	}
	else if ((x >= 0.5) && (x < 1)) {
		return 4 * (1 - x) / 3;
	}
	else {
		return 0;
	}*/
	return sin(PI * x / L) * exp(-PI * PI * t / (L*L));
	//return (uL - u0) * x + u0;

}
double K(double x, const int okey = 1, double y0 = 0) {

	switch (okey) {

	case 1:

		if (x <= x1) {
			return k1;
		}

		else if ((x > x1) && (x < x2)) {
			return k1 * (x - x2) / (x1 - x2) + k2 * (x - x1) / (x2 - x1);
		}

		else {
			return k2;
		}

	case 2:

		return 0.5 * y0 * y0;

	case 3:

		if (x <= 0.5) {
			return 2;
		}

		else if(x >= 0.5) {
			return 1;
		}

	case 4:

		return 2;

	}
}
std::vector<double> progonka(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {

	int n = b.size();
	std::vector<double> alph(n, 0);
	std::vector<double> beth(n, 0);
	std::vector<double> solve(n, 0);
	double tmp;

	//vectinput();
	alph[0] = 0;
	beth[0] = 0;
	tmp = b[0];
	alph[1] = c[0] / tmp;
	beth[1] = d[0] / tmp;

	//std::cout << alph[1] << " " << beth[1] << "\n";
	for (int i = 1; i < n - 1; i++) {
		tmp = b[i] - a[i] * alph[i];
		alph[i + 1] = c[i] / tmp;//(b[i] - a[i] * alph[i]);
		beth[i + 1] = (d[i] + a[i] * beth[i]) / tmp; // (b[i] - a[i] * alph[i]);
	}

	solve[n - 1] = (d[n - 1] + a[n - 1] * beth[n - 1]) / (b[n - 1] - a[n - 1] * alph[n - 1]);

	for (int i = n - 2; i >= 0; i--) {
		solve[i] = alph[i + 1] * solve[i + 1] + beth[i + 1];
	}

	return solve;
}
class heatequation {
private:

	double T; //отрезок по времени
	int key_init_condition;//ключ на начальное условие
	int key_bound_condition;//ключ на граничные условия
	int P_key1; //ключ на функцию P;
	int P_key2;
public:

	heatequation(const double T1, int keyinit, int keybound, int keyp1 = 1, int keyp2 = 1) { //инициализация

		T = T1;
		key_init_condition = keyinit;
		key_bound_condition = keybound;
		P_key1 = keyp1;
		P_key2 = keyp2;

	};
	
	void integrate_interpolation_method(const double N, const double M, const double sigma, const int keyp = 1) {//N - по времени, М - по оси икс.

		double h = L / M;
		double tau = T / N;
		std::ofstream out;
		std::ofstream outenergy;
		std::vector<double> yprev(M + 1, 0);
		std::vector<double> tmpvec(M + 1, 0);
		std::vector<double> y(M + 1, 0);
		std::vector<double> A(M + 1, 0);
		std::vector<double> B(M + 1, 0);
		std::vector<double> C(M + 1, 0);
		std::vector<double> F(M + 1, 0);
		double t, x, w0, a0, w1, w2, wn, kappa, mu, am, a1, a2, tmp, sum1, sum2;
		double uu0 = 10;

		for (int i = 0; i < M + 1; i++) {
			x = i * h;
			yprev[i] = u_init(x, key_init_condition);
			//std::cout << x << " " << yprev[i] << "\n";
		}

		out.open("D:\\vichy_umf_output_laba2.txt");
		outenergy.open("D:\\vichy_umf_output_laba2_energy.txt");

		if (out.is_open() && outenergy.is_open()) {

			t = 0;

			out << h << "\n";
			out << t << " ";
			for (int i = 0; i < M + 1; i++) {
				out << yprev[i] << " ";
			}
			out << std::endl;

			for (int i = 1; i < N + 1; i++) {

				t = tau * i;

				out << t << " ";
				outenergy << t << " ";

				if (sigma == 0) {

					if (key_bound_condition == 1) {
						y[0] = u0;
						y[M] = u0;
					}

					else if (key_bound_condition == 2) {
						y[0] = u0;
						
						am = 2 * K(L, keyp) * K(L - h, keyp) / (K(L, keyp) + K(L - h, keyp));
						wn = am * (yprev[M] - yprev[M - 1]) / h;
						mu = (c * rho * yprev[M] * h * 0.5 / tau + (P(t - tau, P_key2) - wn)) / (c * rho * h * 0.5 / tau);
						y[M] = -mu;
					}

					else if (key_bound_condition == 3) {
						a0 = 2 * K(h, keyp) * K(0, keyp) / (K(h, keyp) + K(0, keyp));
						w0 = a0 * (yprev[1] - yprev[0]) / h;
						mu = (c * rho * yprev[0] * h * 0.5 / tau + (P(t - tau, P_key1) + w0)) / (c * rho * h * 0.5 / tau);
						y[0] = -mu;

						y[M] = u0;
					}

					else if (key_bound_condition == 4) {
						a0 = 2 * K(h, keyp) * K(0, keyp) / (K(h, keyp) + K(0, keyp));
						w0 = a0 * (yprev[1] - yprev[0]) / h;
						mu = (c * rho * yprev[0] * h * 0.5 / tau + (P(t - tau, P_key1) + w0)) / (c * rho * h * 0.5 / tau);
						y[0] = -mu;

						am = 2 * K(L, keyp) * K(L - h, keyp) / (K(L, keyp) + K(L - h, keyp));
						wn = am * (yprev[M] - yprev[M - 1]) / h;
						mu = (c * rho * yprev[M] * h * 0.5 / tau + (P(t - tau, P_key2) - wn)) / (c * rho * h * 0.5 / tau);
						y[M] = -mu;
					}
					
					else if (key_bound_condition == 5) {
						y[0] = 1;
						y[M] = 0;
					}

					for (int j = 1; j < M; j++) {
						x = j * h;
						a1 = 2 * K(x, keyp) * K(x - h, keyp) / (K(x, keyp) + K(x - h, keyp));
						a2 = 2 * K(x + h, keyp) * K(x, keyp) / (K(x + h, keyp) + K(x, keyp));
						w1 = a1 * (yprev[j] - yprev[j - 1]) / h;
						w2 = a2 * (yprev[j + 1] - yprev[j]) / h;
						y[j] = (c * rho * h * yprev[j] / tau + (w2 - w1)) / (c * rho * h / tau);
						/*x = j * h;
						a1 = 2 * K(x) * K(x - h) / (K(x) + K(x - h));
						a2 = 2 * K(x + h) * K(x) / (K(x + h) + K(x));
						w1 = a1 * (yprev[j] - yprev[j - 1]) / h;
						w2 = a2 * (yprev[j + 1] - yprev[j]) / h;
						A[j] = sigma * a1 / h;
						B[j] = sigma * a2 / h;
						C[j] = A[j] + B[j] + c * rho * h / tau;
						F[j] = c * rho * h * yprev[j] / tau + (1 - sigma) * (w2 - w1);*/
					}
				}
				else {
					if (key_bound_condition == 1) {

						A[0] = 0;
						B[0] = 0;
						C[0] = -1;
						F[0] = -u0;

						A[M] = 0;
						B[M] = 0;
						C[M] = -1;
						F[M] = -u0;
					}

					else if (key_bound_condition == 2) {

						A[0] = 0;
						B[0] = 0;
						C[0] = -1;
						F[0] = -u0;

						am = 2 * K(L, keyp) * K(L - h, keyp) / (K(L, keyp) + K(L - h, keyp));
						wn = am * (yprev[M] - yprev[M - 1]) / h;
						kappa = sigma * am / h / (c * rho * h * 0.5 / tau + sigma * am / h);
						mu = (c * rho * yprev[M] * h * 0.5 / tau + sigma * P(t, P_key2) + (1 - sigma) * (P(t - tau, P_key2) - wn)) / (c * rho * h * 0.5 / tau + sigma * am / h);

						A[M] = -kappa;
						B[M] = 0;
						C[M] = -1;
						F[M] = -mu;
					}

					else if (key_bound_condition == 3) {

						a0 = 2 * K(h, keyp) * K(0, keyp) / (K(h, keyp) + K(0, keyp));
						w0 = a0 * (yprev[1] - yprev[0]) / h;
						kappa = sigma * a0 / h / (c * rho * h * 0.5 / tau + sigma * a0 / h);
						mu = (c * rho * yprev[0] * h * 0.5 / tau + sigma * P(t, P_key1) + (1 - sigma) * (P(t - tau, P_key1) + w0)) / (c * rho * h * 0.5 / tau + sigma * a0 / h);
						A[0] = 0;
						B[0] = -kappa;
						C[0] = -1;
						F[0] = -mu;

						A[M] = 0;
						B[M] = 0;
						C[M] = -1;
						F[M] = -u0;
					}
					else if (key_bound_condition == 4) {

						a0 = 2 * K(h, keyp) * K(0, keyp) / (K(h, keyp) + K(0, keyp));
						w0 = a0 * (yprev[1] - yprev[0]) / h;
						kappa = sigma * a0 / h / (c * rho * h * 0.5 / tau + sigma * a0 / h);
						mu = (c * rho * yprev[0] * h * 0.5 / tau + sigma * P(t, P_key1) + (1 - sigma) * (P(t - tau, P_key1) + w0)) / (c * rho * h * 0.5 / tau + sigma * a0 / h);

						A[0] = 0;
						B[0] = -kappa;
						C[0] = -1;
						F[0] = -mu;

						am = 2 * K(L, keyp) * K(L - h, keyp) / (K(L, keyp) + K(L - h, keyp));
						wn = am * (yprev[M] - yprev[M - 1]) / h;
						kappa = sigma * am / h / (c * rho * h * 0.5 / tau + sigma * am / h);
						mu = (c * rho * yprev[M] * h * 0.5 / tau + sigma * P(t, P_key2) + (1 - sigma) * (P(t - tau, P_key2) - wn)) / (c * rho * h * 0.5 / tau + sigma * am / h);

						A[M] = -kappa;
						B[M] = 0;
						C[M] = -1;
						F[M] = -mu;

					}
					else if (key_bound_condition == 5) {

						A[0] = 0;
						B[0] = 0;
						C[0] = -1;
						F[0] = -1;

						A[M] = 0;
						B[M] = 0;
						C[M] = -1;
						F[M] = 0;

					}

					for (int j = 1; j < M; j++) {
						 
						x = j * h;
						a1 = 2 * K(x, keyp) * K(x - h, keyp) / (K(x, keyp) + K(x - h, keyp));
						a2 = 2 * K(x + h, keyp) * K(x, keyp) / (K(x + h, keyp) + K(x, keyp));
						//std::cout << a1 << " " << a2 << "\n";
						w1 = a1 * (yprev[j] - yprev[j - 1]) / h;
						w2 = a2 * (yprev[j + 1] - yprev[j]) / h;

						A[j] = sigma * a1 / h;
						B[j] = sigma * a2 / h;
						C[j] = A[j] + B[j] + c * rho * h / tau;
						F[j] = c * rho * h * yprev[j] / tau + (1 - sigma) * (w2 - w1);

					}

					y = progonka(A, C, B, F);
				}
			
				sum1 = 0;
				sum2 = 0;

				for (int k = 0; k < M + 1; k++) {
					sum1 += y[k];
					sum2 += yprev[k];
					out << y[k] << " ";
				}

				out << "\n";
				outenergy << h * (abs(sum1 - 0.5 * (y[M] + y[0])) - abs(sum2 - 0.5 * (yprev[M] + yprev[0]))) << "\n";

				yprev = y;
			}
		}

		out.close();
		outenergy.close();
	}

	void iteration_method(const double N, const double M) { //N - по времени, М - по оси икс.

		double h = L / M;
		double tau = T / N;
		std::ofstream out;
		std::vector<double> yprev(M + 1, 0);
		std::vector<double> tmpvec(M + 1, 0);
		std::vector<double> y(M + 1, 0);
		std::vector<double> A(M + 1, 0);
		std::vector<double> B(M + 1, 0);
		std::vector<double> C(M + 1, 0);
		std::vector<double> F(M + 1, 0);
		double t, x, w1, w2, wn, kappa, mu, am, a1, a2;
		double uu0 = 10;

		for (int i = 0; i < M + 1; i++) {
			t = i * tau;
			yprev[i] = 0;//u_init(x, key_init_condition);
		}

		out.open("D:\\vichy_umf_output_laba2_itermethod.txt");

		if (out.is_open()) {

			t = 0;

			out << h << "\n";
			out << t << " ";
			for (int k = 0; k < M + 1; k++) {
				out << yprev[k] << " ";
			}
			out << std::endl;

			for (int i = 1; i < N + 1; i++) {

				t = i * tau;
				out << t << " ";
				//out << y[0] << " ";

				x = 0;
				A[0] = 0;
				B[0] = 0;
				C[0] = -1;
				F[0] = -uu0 * sqrt(t);

				A[M] = 0;
				C[M] = -1;
				B[M] = 0;
				F[M] = 0;

				for (int s = 0; s < 2; s++) {

					for (int j = 1; j < M; j++) {

						x = j * h;
						//std::cout << j+1 << " " << x << std::endl;
						a1 = 0.5 * (K(x, 2, yprev[j]) + K(x, 2, yprev[j - 1]));
						a2 = 0.5 * (K(x, 2, yprev[j + 1]) + K(x, 2, yprev[j]));

						A[j] = a1 / h;
						B[j] = a2 / h;
						C[j] = A[j] + B[j] + c * rho * h / tau;
						F[j] = c * rho * h * yprev[j] / tau;

					}

					y = progonka(A, C, B, F);

					yprev = y;

				}

				for (int j = 0; j < y.size(); j++) {
					//y= tmpvec[j];
					out << y[j] << " ";
				}
				out << std::endl;

				yprev = y;

			}
		}

		out.close();
	}

	void modify_iteration_method(const double N, const double M) { // для нахождения необходимо для сходимости количества итераций, лучше взять эпсилон поменьше

		double h = L / M;
		double tau = T / N;
		std::ofstream out;
		std::ofstream iterout;
		double eps = 10E-9;
		std::vector<double> yprev(M + 1, 0);
		std::vector<double> tmpvec(M + 1, 0);
		std::vector<double> y(M + 1, 0);
		std::vector<double> A(M + 1, 0);
		std::vector<double> B(M + 1, 0);
		std::vector<double> C(M + 1, 0);
		std::vector<double> F(M + 1, 0);
		double t, x, w1, w2, wn, kappa, mu, am, a1, a2;
		int iter;
		double uu0 = 10;

		for (int i = 0; i < M + 1; i++) {
			t = i * tau;
			yprev[i] = 0;//u_init(x, key_init_condition);
		}

		out.open("D:\\vichy_umf_output_laba2_itermethod.txt");
		iterout.open("D:\\vichy_umf_output_laba2_itermethod_iteration.txt");

		if ((out.is_open()) && (iterout.is_open())) {

			t = 0;

			out << h << "\n";
			out << t << " ";
			for (int k = 0; k < M + 1; k++) {
				out << yprev[k] << " ";
			}
			out << std::endl;

			for (int i = 1; i < N + 1; i++) {
				t += tau;
				iter = 0;
				y[0] = u0;

				out << t << " ";
				iterout << t << " ";
				//out << y[0] << " ";

				x = 0;
				A[0] = 0;
				B[0] = 0;
				C[0] = -1;
				F[0] = -uu0 * sqrt(t);

				A[M] = 0;
				C[M] = -1;
				B[M] = 0;
				F[M] = 0;

				tmpvec = yprev;

				while(sqr(y, tmpvec) > eps) {

					iter += 1;
					tmpvec = yprev;

					//std::cout << sqr(y, yprev) << " " << iter << "\n";
					for (int j = 1; j < M; j++) {

						x = j * h;
						//std::cout << j+1 << " " << x << std::endl;
						a1 = 0.5 * (K(x, 2, yprev[j]) + K(x, 2, yprev[j - 1]));
						a2 = 0.5 * (K(x, 2, yprev[j + 1]) + K(x, 2, yprev[j]));

						A[j] = a1 / h;
						B[j] = a2 / h;
						C[j] = A[j] + B[j] + c * rho * h / tau;
						F[j] = c * rho * h * yprev[j] / tau;

					}

					y = progonka(A, C, B, F);
					yprev = y;

				}

				for (int j = 0; j < y.size(); j++) {
					//y= tmpvec[j];
					out << y[j] << " ";
				}
				iterout << iter << "\n";
				out << std::endl;

				yprev = y;

			}
		}

		out.close();
		iterout.close();
	}

	double modify_integrate_interpolation_method(const double N, const double M, const double sigma, const int keyp) {//N - по времени, М - по оси икс.

		double h = L / M;
		double tau = T / N;
		std::ofstream out;
		std::ofstream outenergy;
		std::vector<double> arrnorm;
		std::vector<double> yprev(M + 1, 0);
		std::vector<double> tmpvec(M + 1, 0);
		std::vector<double> y(M + 1, 0);
		std::vector<double> A(M + 1, 0);
		std::vector<double> B(M + 1, 0);
		std::vector<double> C(M + 1, 0);
		std::vector<double> F(M + 1, 0);
		double t, x, w0, a0, w1, w2, wn, kappa, mu, am, a1, a2, tmp, norm, tmpp, normh;

		normh = 0;
		norm = 0;
		double uu0 = 10;

		for (int i = 0; i < M + 1; i++) {
			x = i * h;
			yprev[i] = sin(PI * x / L);
		}

		out.open("D:\\vichy_umf_output_laba2.txt");

		if (out.is_open()) {

			t = 0;
			
			for (int i = 1; i < N + 1; i++) {

				t = tau * i;
				
				if (sigma == 0) {

					y[0] = 0;
					y[M] = 0;

					for (int j = 1; j < M; j++) {

						x = j * h;
						a1 = 0.5 * (K(x, keyp) + K(x - h, keyp));
						a2 = 0.5 * (K(x + h, keyp) + K(x, keyp));
						w1 = a1 * (yprev[j] - yprev[j - 1]) / h;
						w2 = a2 * (yprev[j + 1] - yprev[j]) / h;

						y[j] = (c * rho * h * yprev[j] / tau + (w2 - w1)) / (c * rho * h / tau);
						/*x = j * h;
						a1 = 2 * K(x) * K(x - h) / (K(x) + K(x - h));
						a2 = 2 * K(x + h) * K(x) / (K(x + h) + K(x));
						w1 = a1 * (yprev[j] - yprev[j - 1]) / h;
						w2 = a2 * (yprev[j + 1] - yprev[j]) / h;
						A[j] = sigma * a1 / h;
						B[j] = sigma * a2 / h;
						C[j] = A[j] + B[j] + c * rho * h / tau;
						F[j] = c * rho * h * yprev[j] / tau + (1 - sigma) * (w2 - w1);*/

					}
				}
				else {

					A[0] = 0;
					B[0] = 0;
					C[0] = -1;
					F[0] = 0;

					A[M] = 0;
					B[M] = 0;
					C[M] = -1;
					F[M] = 0;

					for (int j = 1; j < M; j++) {

						x = j * h;
						a1 = 0.5 * (K(x, keyp) + K(x - h, keyp));
						a2 = 0.5 * (K(x + h, keyp) + K(x, keyp));
						w1 = a1 * (yprev[j] - yprev[j - 1]) / h;
						w2 = a2 * (yprev[j + 1] - yprev[j]) / h;

						A[j] = sigma * a1 / h;
						B[j] = sigma * a2 / h;
						C[j] = A[j] + B[j] + c * rho * h / tau;
						F[j] = c * rho * h * yprev[j] / tau + (1 - sigma) * (w2 - w1);

					}

					y = progonka(A, C, B, F);
				}

				if (t <= 0.1) {

					norm = 0;

					for (int i = 0; i < y.size(); i++) {
						tmpp = fabs(y[i] - exact_u(i*h, t));
						if (norm < tmpp) {
							norm = tmpp;
						}
					}

					arrnorm.push_back(norm);
				}

				yprev = y;
			}

			normh = 0;

			for (int i = 0; i < arrnorm.size(); i++) {
				if (normh < arrnorm[i]) {
					normh = arrnorm[i];
				}
			}

		}

		out.close();

		return normh;
	}
};



void order(double sigma) {

	std::cout << sigma << "\n";
	heatequation umfequation(1, 5, 1);
	int m, n;
	int div = 2;
	double prev, ths;
	prev = 1;

	if (not sigma) {
		m = 4;
		n = 50;
	}

	else if(sigma == 0.5){
		//sigma = 0.5 +/- ok m = 10, n= 50;
		m = 10;
		n = 25;
	}

	else if (sigma == 1) {
		m = 10;
		n = 20;
	}

	/*if (flag) {
		m = 4;
		n = 1000;
	}*/

	std::cout << "sigma " << sigma << "\n";

	for (int p = 0; p < 7; p++) {

		std::cout << n << " " << m << "\n";

		ths = umfequation.modify_integrate_interpolation_method(n, m, sigma, 4);
		
		if (sigma == 0.5) {
			n *= 2;
			m *= 2;
		}

		else {
			n *= 4;
			m *= 2;
		}

		if (not p) {
			std::cout << "step " << p + 1 << " " << ths << " --- " << "---" << "\n";
		}

		else {
			std::cout << "step " << p + 1 << " " << ths << " " << prev/ths << " " << log(prev / ths) / log(2) << "\n";
		}

		prev = ths;
		
	}

	std::cout << "\n";
}

//при инциализации: 1 аргумент - длина отрезка по времени, 2 аргумент - ключ на начальное условие, 3 аргумент - на граничные условие, 4 и 5 аргументы - необязательные, обозначают
//поток слева и справа соответственно

void test0() {
	heatequation umfequation(1, 1, 1);
	umfequation.integrate_interpolation_method(2000, 100, 0);
}

void test05() {
	heatequation umfequation(1, 1, 1);
	umfequation.integrate_interpolation_method(100, 100, 0.5);
}

void test1() {
	heatequation umfequation(1, 1, 1);
	umfequation.integrate_interpolation_method(100, 100, 1);
}

void test11() {
	heatequation umfequation(1, 1, 2, 1, 4);
	umfequation.integrate_interpolation_method(100, 100, 0.5);
}

void test2() {
	heatequation umfequation(1, 1, 3, 4, 5);
	umfequation.integrate_interpolation_method(100, 100, 1);
}

void test3() {
	heatequation umfequation(1, 1, 4, 4, 4);
	umfequation.integrate_interpolation_method(100, 200, 1);
}

void test4() {
	heatequation umfequation(1, 1, 4, 5, 5);
	umfequation.integrate_interpolation_method(100, 100, 0.5);
}

void test5() {
	heatequation umfequation(1, 2, 1, 5, 5);
	umfequation.iteration_method(2 * 1000, 50);
}

void test55() {
	heatequation umfequation(1, 2, 1, 5, 5);
	umfequation.modify_iteration_method(100, 10);
}

void test6() {
	heatequation umfequation(1, 1, 4, 5, 5);
	umfequation.integrate_interpolation_method(100, 100, 1);
}

void test8() {
	heatequation umfequation(1, 1, 2, 4, 3);
	umfequation.integrate_interpolation_method(100, 100, 0.5);
}

void test7() {
	heatequation umfequation(50, 1, 1, 4, 5);
	umfequation.integrate_interpolation_method(100, 100, 0.5);
}

void test9() {
	heatequation umfequagtion(1, 4, 1);
	umfequagtion.integrate_interpolation_method(10, 10, 1, 3);
}

void test_mytest() {
	heatequation umfequation(10, 3, 3, 4, 4);
	umfequation.integrate_interpolation_method(100, 100, 0.5);
}

void test10() { //немонотонный расчет
	heatequation umfequation(1, 1, 2, 5, 4);
	umfequation.integrate_interpolation_method(10, 100, 0.5);
}

void test_order() {
	order(0);
	order(0.5);
	order(1);
}

void anotheronetset() {
	heatequation umfequation(1, 4, 5);
	umfequation.integrate_interpolation_method(10, 100, 1, 3);
}

int main()
{
	//test9();
	test10();
	//test5();
	//test_order();
	//anotheronetset();
	//test_mytest();
}