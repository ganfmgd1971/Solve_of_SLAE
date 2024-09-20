#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
using namespace std;
//const double PI = 3.141592653589793;
double ddg(double x) {
    if (x < 0) {
        return -8 / ((-2 + x) * (-2 + x) * (-2 + x));
    } else if (x == 0) {
        return  1;
    }
    else {
        return exp(x);
    }
}
double f(double x, int k) {
    switch (k)
    {
    case 1:
        return x * x; //-1, 1
    case 2:
        return 1 / (1 + x * x); //-1, 1
    case 3:
        return 1 / atan(1 + 10 * x * x); // -3, 3
    case 4:
        return pow(4 * x * x * x + 2 * x * x - 4 * x + 2, sqrt(2)) + asin(1 / (5 + x - x * x)) - 5; //-1, 1
    case 5:
        return exp(x);
    case 6:
        return tan((2 * x * x * x * x - 5 * x + 6) / 8) + atan((7 * x * x - 11 * x + 1 - sqrt(2)) / (-7 * x * x + 11 * x + sqrt(2))); //-1 0
    case 7:
        return 1;
    case 8:
        return pow(4, x);
    case 9:
        if (x < 0) {
            return (1 + x / 2) / (1 - x / 2);
        }
        else {
            return exp(x);
        }

    case 10:
        double lla = -1;
        double rrb = 1;
        double p = ddg(lla) * (x - rrb) * (x - rrb) * (x - rrb) / (6 * (rrb - lla)) - ddg(rrb) * (x - lla) * (x - lla) * (x - lla) / (6 * (rrb - lla));
        if (x < 0) {
            return (1 + x / 2) / (1 - x / 2) + p;
        }
        else {
            return exp(x) + p;
        }
    }
    /*R1: = Tan((2 * Degree(x, 4) - 5 * x + 6) / 8);
    R2: = ArcTan((7 * Degree(x, 2) - 11 * x + 1 - Sqrt(2)) /
        (-7 * Degree(x, 2) + 11 * x + Sqrt(2)));
VarF: = R1 + R2;*/
}
void printvector(vector<double> &y, int N) {
    for (int i = 0; i < N; i++) {
        cout << y[i] << " ";
    }
    cout << endl;
}
class interp {
    int N;
    vector<double> X;
    vector<double> F;
    vector<double> A; //коэффициенты для сплайн интерполяции
    vector<double> B;
    vector<double> C;
    vector<double> D;
public:
    void construct(int P) {
        N = P;
        X.resize(N + 1);
        F.resize(N + 1);
    }
    void input(string str) {
        ifstream in(str);
        if (in.is_open()) {
            in >> N;
            X.resize(N);
            F.resize(N);
            for (int i = 0; i < N; i++) {
                in >> X[i];
            }
        }
        printvector(X, N);
    }
    void make_even_grid(double a, double b, int key) {
        double h = (b - a) / N;
        //N = M;
        for (int i = 0; i < N + 1; i++) {
            X[i] = a + i * h;
            F[i] = f(X[i], key);
        }
        //printvector(X, N);
        /*printvector(F, N); */
    }
    void make_chebyshev(double a, double b, int key) {
        //N = M;
        double PI = acos(-1.0);
        for (int i = 0; i < N + 1; i++) {
            X[i] = (a + b) / 2 + ((b - a) / 2) * cos((2 * i + 1) * PI / (2 * (N + 1)));
            F[i] = f(X[i], key);
        }
        /*printvector(X, N);
        printvector(F, N);*/
    }
    void output(string str) {
        ofstream out;
        out.open(str);      // открываем файл для записи
        if (out.is_open())
        {
            for (int i = 0; i < N + 1; i++) {
                out << F[i] << " ";
            }
            out << endl;
            for (int i = 0; i < N + 1; i++) {
                out << X[i] << " ";
            }
            //out.close();
        }
        out.close();
    }
    double Lagrange(double x) {
        //vector<double> Lagrange_polynom(N+1, 0);
        double res = 0;
        vector<double> ck(N + 1, 0);
        for (int i = 0; i <= N; i++) {
            double num = 1;
            double denom = 1;
            for (int j = 0; j <= N; j++) {
                if (i != j) {
                    num *= x - X[j];
                    denom *= X[i] - X[j];
                }
            }
            ck[i] = num / denom;
            res += ck[i] * F[i];
        }
        return res;
    }
    /*vector<double> Lagrange(double x) {
        vector<double> Lagrange_polynom(N+1, 0);
        vector<double> ck(N + 1, 0);
        for (int i = 0; i <= N; i++) {
            double num = 1;
            double denom = 1;
            for (int j = 0; j <= N; j++) {
                if (i != j) {
                    num *= x - X[j];
                    denom *= X[i] - X[j];
                }
            }
            ck[i] = num / denom;
            res += ck[i] * F[i];
        }
        return Lagrange;
    }*/
    vector<double> progonka(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d) {
        vector<double> z(10, 0);
        int n = b.size();
        vector<double> alpha(n, 0);
        vector<double> betha(n, 0);
        vector<double> y(n, 0);
        //vectinput();
        double tmp = b[0];
        alpha[1] = c[0] / tmp;
        betha[1] = d[0] / tmp;
        for (int i = 1; i < n - 1; i++) {
            tmp = b[i] - a[i] * alpha[i];
            alpha[i + 1] = c[i] / tmp;
            betha[i + 1] = (d[i] + a[i] * betha[i]) / tmp;
        }
        //std::cout << "AAA" << std::endl;
        y[n - 1] = (d[n - 1] + a[n - 1] * betha[n - 1]) / (b[n - 1] - a[n - 1] * alpha[n - 1]);
        //std::cout << "BBB" << std::endl;
        for (int i = n - 2; i >= 0; --i) {
            y[i] = alpha[i + 1] * y[i + 1] + betha[i + 1];
        }
        //std::cout << "CCC" << std::endl;
        //printvector(y, n);
        return y;
    }
    void coef_spline_interp() {
        //int n = N;
        int m = N;
        vector<double> g(m);
        vector<double> h(m);
        A.resize(m);
        B.resize(m);
        C.resize(m - 1);
        D.resize(m);
        vector<double> tmpvec1(m + 1, 0);
        vector<double> tmpvec2(m - 1, 0);
        vector<double> atmp(m - 1, 0);
        vector<double> btmp(m - 1, 0);
        vector<double> ctmp(m - 1, 0);
        vector<double> dtmp(m - 1, 0);
        for (int i = 0; i < m; i++) {
            h[i] = X[i + 1] - X[i];
            A[i] = F[i];
            B[i] = 0;
            D[i] = 0;
            g[i] = (F[i + 1] - F[i]) / h[i];
        }
        for (int i = 0; i < m - 1; i++) {
            C[i] = 0;
            atmp[i] = h[i];
            btmp[i] = -2 * (h[i + 1] + h[i]);
            ctmp[i] = h[i + 1];
            dtmp[i] = -3 * (g[i + 1] - g[i]);
        }
        atmp[0] = 0;
        ctmp[m - 2] = 0;
        C = progonka(atmp, btmp, ctmp, dtmp);
        C.emplace(C.begin(), 0);
        C.push_back(0);
        //printvector(C, m+1);
        tmpvec1[m] = 0;
        for (int i = 0; i < m; i++) {
            B[i] = g[i] - (C[i + 1] + 2 * C[i]) * h[i] / 3;
            D[i] = (C[i + 1] - C[i]) / (3 * h[i]);
        }
    }
    double calc_spline(double x) {
        int m = N;
        int j = 0;
        for (int i = 0; i < m; i++) {
            if (x > X[i]) {
                j = i;
            }
        }
        return A[j] + B[j] * (x - X[j]) + C[j] * (x - X[j]) * (x - X[j]) + D[j] * (x - X[j]) * (x - X[j]) * (x - X[j]);
    }
    void output_Lagrange(double a, double b, int key, string str) {
        vector<double> xhelp(N, 0);
        vector<double> Lagr(N, 0);
        double h = (b - a) / N;
        double tmp = a;
        double epsilon = 0;
        for (int i = 0; i < N; i++) {
            xhelp[i] = tmp + h / 2;
            tmp += h;
        }
        for (int i = 0; i < N; i++) {
            Lagr[i] = Lagrange(xhelp[i]);
            /*tmp = abs(Lagr[i] - f(xhelp[i], key));
            if (tmp > epsilon) {
                epsilon = tmp;
            }*/
        }
        //std::cout << epsilon << std::endl;
        //printvector(xhelp, N);
        ofstream out;
        out.open(str);      // открываем файл для записи
        if (out.is_open())
        {
            for (int i = 0; i < N; i++) {
                out << Lagr[i] << " ";
            }
            out << endl;
            for (int i = 0; i < N; i++) {
                out << xhelp[i] << " ";
            }
            //out.close();
        }
        out.close();
    }
    void spec_output_Lagrange(double a, double b, int key, string str) {
        vector<double> xhelp(2 * N, 0);
        vector<double> Lagr(2 * N, 0);
        double h = (b - a) / N;
        double tmp = a;
        double epsilon = 0;
        for (int i = 0; i < 2 * N; i++) {
            xhelp[i] = tmp + h / 4;
            tmp += h / 2;
        }
        for (int i = 0; i < 2 * N; i++) {
            Lagr[i] = Lagrange(xhelp[i]);
            epsilon += abs(Lagr[i] - f(xhelp[i], key));
        }
        //std::cout << epsilon << std::endl;
        printvector(xhelp, 2 * N);
        ofstream out;
        out.open(str);      // открываем файл для записи
        if (out.is_open())
        {
            for (int i = 0; i < 2 * N; i++) {
                out << Lagr[i] << " ";
            }
            out << endl;
            for (int i = 0; i < 2 * N; i++) {
                out << xhelp[i] << " ";
            }
            //out.close();
        }
        out.close();
    }
    void output_spline(double a, double b, int key, string str) {
        vector<double> xhelp(N, 0);
        vector<double> spline(N, 0);
        double h = (b - a) / N;
        double tmp = a;
        double epsilon = 0;
        for (int i = 0; i < N; i++) {
            xhelp[i] = tmp + h / 2;
            tmp += h;

        }
        coef_spline_interp();
        for (int i = 0; i < N; i++) {
            spline[i] = calc_spline(xhelp[i]);
            //epsilon += abs(spline[i] - f(xhelp[i], key));
            //std::cout << epsilon << std::endl;
        }
        //printvector(xhelp, N);
        //std::cout << epsilon << std::endl;
        ofstream out;
        out.open(str);      // открываем файл для записи
        if (out.is_open())
        {
            for (int i = 0; i < N; i++) {
                out << spline[i] << " ";
            }
            out << endl;
            for (int i = 0; i < N; i++) {
                out << xhelp[i] << " ";
            }
            //out.close();
        }
        out.close();
    }
    double check_error_Lagrange(double a, double b, int key) {
        vector<double> xhelp(N, 0);
        double h = (b - a) / N;
        double tmp = a;
        double epsilon = 0;
        for (int i = 0; i < N; i++) {
            xhelp[i] = tmp + h / 2;
            tmp += h;
            //std::cout << xhelp[i] << std::endl;
        }
        for (int i = 0; i < N; i++) {
            tmp = abs(Lagrange(xhelp[i]) - f(xhelp[i], key));
            //std::cout << epsilon << std::endl;
            if (tmp > epsilon) {
                epsilon = tmp;
            }
        }
        return epsilon;
    }
    double check_error_spline(double a, double b, int key) {
        vector<double> xhelp(N, 0);
        double h = (b - a) / N;
        double tmp = a;
        double epsilon = 0;
        for (int i = 0; i < N; i++) {
            xhelp[i] = tmp + h / 2;
            tmp += h;
            //std::cout << xhelp[i] << std::endl;
        }
        coef_spline_interp();
        for (int i = 0; i < N; i++) {
            //std::cout << calc_spline(xhelp[i]) << " " << f(xhelp[i], key) << std::endl;
            tmp = abs(calc_spline(xhelp[i]) - f(xhelp[i], key));
            //std::cout << epsilon << std::endl;
            if (tmp > epsilon) {
                epsilon = tmp;
            }
        }
        return epsilon;
    }
    void make_super_grid(int p, double a, double b, int key) {
        N = N * p;
        //double q = 1 / p;
        double h = (b - a) / N;
        /*X.resize(N);
        F.resize(N);*/
        for (int i = 0; i < N; i++) {
            X[i] = a + i * h;
            F[i] = f(X[i], key);
        }
    }
    void error() {
        double t = 0.3;

    }
};
int main()
{
    interp C;
    string str1 = "D:\\output_laba4_Lagrange_even.txt";
    string str2 = "D:\\output_laba4_Lagrange_chebyshev.txt";
    string str3 = "D:\\output_laba4_spline_even.txt";
    string str4 = "D:\\output_laba4_spline_chebyshev.txt";
    int N = 5;
    int M;
    C.construct(N);
    double la = 0;
    double rb = 5;
    int key = 8;
    C.make_even_grid(la, rb, key);
    std::cout << C.check_error_Lagrange(la, rb, key) << std::endl;
    std::cout << C.Lagrange(0.5) - 2 << std::endl;
    C.output_Lagrange(la, rb, key, str1);
    //f(x) = 4^x интерполяция равномерная сетка лагранж в точке 0.5
    /*C.make_even_grid(la, rb, key);
    std::cout << C.Lagrange(0.5) << std::endl;
    std::cout << C.check_error_Lagrange(la, rb, key) << std::endl;
    C.spec_output_Lagrange(la, rb, key, str1);*/
    //C.output_Lagrange(-1, 1, 2, str1);
    ////C.output(str1);
    //std::cout << C.check_error_spline(-1, 1, 2) << std::endl;
    //C.output_spline( - 1, 1, 2, str3);
    ////C.output(str3);
    //C.make_chebyshev(-1, 1, 2);
    //std::cout << C.check_error_Lagrange(-1, 1, 2) << std::endl;
    //C.output_Lagrange(-1, 1, 2, str2);
    //C.output(str2);
    /*C.make_even_grid(la, rb, key);
    std::cout << C.check_error_Lagrange(la, rb, key) << std::endl;
    C.output_Lagrange(la, rb, key, str1);*/
   /* C.make_even_grid(la, rb, key);
    std::cout << C.check_error_Lagrange(la, rb, key) << std::endl;
    C.output_Lagrange(la, rb, key, str1);
    C.make_chebyshev(la, rb, key);
    std::cout << C.check_error_Lagrange(la, rb, key) << std::endl;
    C.output_Lagrange(la, rb, key, str2);
    C.make_even_grid(la, rb, key);
    std::cout << C.check_error_spline(la, rb, key) << std::endl;
    C.output_spline(la, rb, key, str3);*/
    //double p = 1;
    //double tmp = 1;
    //double eps;
    //for (int j = 0; j < 25; j++) {
    //    std::cout << j + 1 << std::endl;
    //    M = (j + 1) * N;
    //    C.construct(M);
    //    C.make_even_grid(la, rb, key);
    //    eps = C.check_error_spline(la, rb, key);
    //    //std::cout << eps;
    //    std::cout << tmp / eps << std::endl;
    //    tmp = eps;
    //    //C.output_spline(la, rb, key, str3);
    //    std::cout << std::endl;
    //}
   /* C.make_even_grid(la, rb, key);
    std::cout << C.check_error_spline(la, rb, key) << std::endl;
    C.output_spline(la, rb, key, str3);*/
    /*double tmp = 1;
    int p;
    for (int i = 1; i < 6; i++) {
        p = std::pow(2, i);
        N = N * p;
        interp D;
        D.construct(N);
        D.even_grid(la, rb, key);
        double delta = D.check_error_spline(la, rb, key);
        std::cout << "error " << delta << " ont " << tmp/delta<<  std::endl;
        tmp = delta;
    }
    C.output_spline(la, rb, key, str3);*/
}