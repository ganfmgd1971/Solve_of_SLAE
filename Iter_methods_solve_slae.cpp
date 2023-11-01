// ConsoleApplication87.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
template <class T> using matrix = std::vector<std::vector<T>>;
template<typename T>
T evaluatemiddle(matrix<T>& A, int n) {
    T res = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res += abs(A[i][j]);
        }
    }
    res = res / (n * n);
    return res;
}
template<typename T>
void printmatrix(matrix<T>& A, int n, T epsilon) {
    //epsilon = epsilon * 100000;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            std::cout << A[i][j] << " ";

        }
        std::cout << std::endl;
    }
}
template<typename T>
T norm(bool flag, matrix<T> A, int n) { //можно передавать имя нормы, в которой мы считаем матрицу
    T tmp, max;
    tmp = 0;
    if (flag) {
        for (int i = 0; i < n; i++) {
            tmp += abs(A[i][0]);
        }
        max = tmp;
        for (int j = 1; j < n; j++) {
            tmp = 0;
            for (int i = 0; i < n; i++) {
                tmp += abs(A[i][j]);
            }
            if (tmp > max) {
                max = tmp;
            }
        }
        return max;
    }
    else {
        //T condition_inf() { //убрать B +
        T tmp, max;
        tmp = 0;
        for (int j = 0; j < n; j++) {
            tmp += abs(A[0][j]);
        }
        max = tmp;
        for (int i = 1; i < n; i++) {
            tmp = 0;
            for (int j = 0; j < n; j++) {
                tmp += abs(A[i][j]);
            }
            if (tmp > max) {
                max = tmp;
            }
        }
        return max;
    }
}
template<typename T>
void printvector(std::vector<T>& y) {
    for (int i = 0; i < y.size(); i++) {
        std::cout << y[i] << std::endl;
    }
    std::cout << std::endl;
}
template<typename T>
void swap(T& a, T& b) {
    T tmp = a;
    a = b;
    b = tmp;
}
template <typename T>
std::vector<T> operator*(const std::vector<T>& lhs, const T& rhs) {
    std::vector<T> result;
    result.resize(lhs.size());
    std::copy(lhs.cbegin(), lhs.cend(), result.begin());
    for (auto& t : result)
        t *= rhs;
    return result;
};
template <typename T>
std::vector<T> operator*(matrix<T>& lhs, std::vector<T>& rhs) {
    std::vector<T> result;
    result.resize(rhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        result[i] = 0;
    }
    for (int i = 0; i < rhs.size(); i++) {
        for (int j = 0; j < rhs.size(); j++) {
            result[i] += lhs[i][j] * rhs[j];
        }
    }
    return result;
}
template<typename T>
T vector_norm(std::vector<T>& a) {
    int n = a.size();
    T result = 0;
    for (int j = 0; j < n; j++) {
        result += a[j] * a[j];
    }
    return sqrt(result);
}
template<typename T>
matrix<T> operator*(matrix<T>& lhs, matrix<T>& rhs) {
    matrix<T> result;
    result.reserve(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        result.emplace_back(std::vector<T>(lhs.size(), 0));
    }
    for (int i = 0; i < lhs.size(); i++) {
        for (int j = 0; j < lhs.size(); j++) {
            result[i][j] = 0;
            for (int k = 0; k < lhs.size(); k++) {
                result[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }
    return result;
};
template<typename T>
matrix<T> operator+(matrix<T>& lhs, matrix<T>& rhs) {
    matrix<T> result;
    result.reserve(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        result.emplace_back(std::vector<T>(lhs.size(), 0));
    }
    for (int i = 0; i < lhs.size(); i++) {
        for (int j = 0; j < lhs.size(); j++) {
            result[i][j] = lhs[i][j] + rhs[i][j];
        }
    }
    return result;
}
template<typename T>
matrix<T> operator-(matrix<T>& lhs, matrix<T>& rhs) {
    matrix<T> result;
    result.reserve(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        result.emplace_back(std::vector<T>(lhs.size(), 0));
    }
    for (int i = 0; i < lhs.size(); i++) {
        for (int j = 0; j < lhs.size(); j++) {
            result[i][j] = lhs[i][j] - rhs[i][j];
        }
    }
    return result;
}
template <typename T>
std::vector<T> operator-(std::vector<T>& lhs, std::vector<T>& rhs) {
    std::vector<T> result(lhs.size(), 0);
    for (int i = 0; i < lhs.size(); i++) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
};
template<typename T>
std::vector<T> operator+(std::vector<T>& lhs, std::vector<T>& rhs) {
    std::vector<T> result;
    result.resize(lhs.size());
    for (int i = 0; i < lhs.size(); i++) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
};
template<typename T>
T column_norm(matrix<T>& a, int column) {
    T column_norm = 0;
    int n = a.size();
    for (int j = 0; j < n; j++)
        column_norm += a[j][column] * a[j][column];
    return sqrt(column_norm);
}
template <typename T>
T vector_projection(matrix<T>& a, int column_a, matrix<T>& q, int column_q) {
    T projection = 0;
    T numerator = 0;
    T norm = 0;
    for (int i = 0; i < a.size(); i++) {
        numerator += a[i][column_a] * q[i][column_q];
        norm += q[i][column_q] * q[i][column_q];
    }
    projection = numerator / sqrt(norm);
    return projection;
}
template <typename T>
T evresidual(std::vector<T>& a, std::vector<T>& b) {
    T tmp = 0;
    for (int i = 0; i < a.size(); i++) {
        tmp += (b[i] - a[i]) * (b[i] - a[i]);
    }
    return sqrt(tmp);
}

template<typename T>
void transpose(matrix<T>& A, matrix<T>& B, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i][j] = A[j][i];
        }
    }
}
template<typename T>
matrix<T> operator*(T& lhs, matrix<T>& rhs) {
    int m = rhs.size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            rhs[i][j] = lhs * rhs[i][j];
        }
    }
    return rhs;
}
template <typename T>
std::vector<T> operator*(T& lhs, std::vector<T>& rhs) {
    int m = rhs.size();
    std::vector <T> result(m, 0);
    for (int i = 0; i < m; i++) {
        result[i] = lhs * rhs[i];
    }
    return result;
}
template <typename T>
matrix<T> Inverse_diag(matrix<T>& D) {
    matrix<T> result;
    int M = D.size();
    result.reserve(M);
    for (int i = 0; i < M; i++) {
        result.emplace_back(std::vector<T>(M, 0));
    }
    for (int i = 0; i < M; i++) {
        D[i][i] = 1 / D[i][i];
    }
}
template <typename T>
class Hsole { //Heterogeneous system of linear equations
    matrix<T> A;
    int N;
    std::vector<T> B;
    std::vector<T> x;
public:
    /*void simple_it(T epsilon, T tau) {
        matrix<T> E;
        E.reserve(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                E.emplace_back(std::vector<T>(N, 0));
            }
            E[i][i] = 1;
        }
        matrix <T> C = -1 * (A - E);
        if ((norm(true, C) <= 1) && (norm(false C) <= 1)) {
            std::cout << "can get solve with help of method of simple iterations " << std::endl;
        }
        else {
            std::cout << "we can't get solve with help of method of simple iterations " << std::endl;
        }

    }*/
    void input(std::string str) { //аргумент строка
        std::ifstream in;
        T tmp;
        in.open(str); //input_laba1.txt"
        if (in.is_open()) {
            in >> N;
            A.reserve(N);
            for (int i = 0; i < N; i++) {
                A.emplace_back(std::vector<T>(N, 0));
            }
            B.resize(N);
            if (in.is_open()) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        in >> tmp;
                        A[i][j] = tmp;
                    }
                    in >> tmp;
                    B[i] = tmp;
                }
            }
            in.close();
        }
    }
    void special_input(int M) { //аргумент строка
        N = M;
        T tmp;
        std::vector<T> z;
        // in.open(str); //input_laba1.txt"
        A.reserve(N);
        for (int i = 0; i < N; i++) {
            A.emplace_back(std::vector<T>(N, 0));
        }
        B.resize(N);
        z.resize(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = 1.0 / (i + 1 + j + 1 - 1);
            }
            z[i] = 1;
        }
        B = A * z;
    }
    void output(std::string str) {
        std::ofstream out;
        out.open(str);
        if (out.is_open()) {
            out << N << std::endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    out << A[i][j] << " ";
                }
                out << B[i] << " " << std::endl;
                //out << x[i] << std::endl;
            }
        }
        out.close();
    }
    std::vector<T> getB() {
        return B;
    }
    std::vector<T> back_gauss(matrix<T>& V, std::vector<T>& F) {
        std::vector<T> y;
        y.resize(N);
        for (int i = 0; i < N; i++) {
            y[i] = 0;
        }
        y[N - 1] = F[N - 1] / V[N - 1][N - 1];
        for (int k = N - 2; k >= 0; k--) {// обратный ход
            y[k] = F[k];
            for (int l = N - 1; l > k; l--) {
                y[k] = y[k] - V[k][l] * y[l];
            }
            y[k] = y[k] / V[k][k];
        }
        return y;
    }
    void Gauss(T epsilon, bool flag) { //параметр эпсилон
        std::vector<T> y;
        T max_element;
        int max_number;
        T tmp, d;
        for (int k = 0; k < N; k++) { //можно добавить флажок на полный частичный выбор
            max_element = abs(A[k][k]);
            max_number = k;
            for (int i = k + 1; i < N; i++) {
                if (abs(A[i][k]) > max_element) {
                    max_element = abs(A[i][k]);
                    max_number = i;
                }
            }
            if (max_element < epsilon * 10000) {
                std::cout << "BAD MATRIX" << std::endl; //добавить про строки
                std::cout << "row " << max_number << "cow " << max_number << std::endl;
                exit(1);
            }
            if (k != max_number) {
                for (int i = 0; i < N; i++) {
                    swap(A[k][i], A[max_number][i]);
                }
                swap(B[k], B[max_number]);
            }
            for (int i = k + 1; i < N; i++) {
                d = A[i][k] / A[k][k];
                for (int j = k; j < N; j++) {
                    A[i][j] -= d * A[k][j];
                }
                B[i] -= d * B[k];
            }
        }
        y = back_gauss(A, B);
        x = y;
        if (flag) {
            std::cout << "solve of Gauss " << std::endl;
            printvector(x);
        }
    }
    void qr_decomposition(T epsilon, matrix<T>& q, matrix<T>& r, int n, bool flag) {
        epsilon = epsilon * 10000;
        matrix<T> t;
        t.reserve(n);
        for (int i = 0; i < n; i++) {
            t.emplace_back(std::vector<T>(n, 0));
        }
        for (int i = 0; i < n; i++) {
            t[i][i] = 1;
        }
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                T sqr = sqrt(r[i][i] * r[i][i] + r[j][i] * r[j][i]);
                if (sqr < epsilon)
                    continue;
                T c = r[i][i] / sqr;
                T s = r[j][i] / sqr;
                for (int k = 0; k < n; k++) {
                    T buf = t[i][k];
                    t[i][k] = c * buf + s * t[j][k];
                    t[j][k] = -1 * s * buf + c * t[j][k];
                }

                for (int k = i; k < n; k++) {
                    T buf = r[i][k];
                    r[i][k] = c * buf + s * r[j][k];
                    if (abs(r[i][k]) < epsilon)
                        r[i][k] = 0;
                    r[j][k] = -1 * s * buf + c * r[j][k];
                    if (abs(r[j][k]) < epsilon)
                        r[j][k] = 0;
                }
            }
        }
        if (abs(r[n - 1][n - 1]) < epsilon) {
            std::cout << "BAD MATRIX" << std::endl;
            exit(1);
        }
        transpose(t, q, n);
        if (flag) {
            std::cout << "Q" << std::endl;
            printmatrix(q, N, epsilon);
            std::cout << std::endl;
            std::cout << "R" << std::endl;
            printmatrix(r, N, epsilon / 100);
            std::cout << std::endl;
        }
    }
    void QR_solve(T epsilon, bool flag) {
        std::vector<T> y;
        matrix<T> r, q, transpq;
        std::vector<T> btmp(N, 0);
        matrix<T> tmp100;
        tmp100.reserve(N);
        r.reserve(N);
        q.reserve(N);
        transpq.reserve(N);
        for (int i = 0; i < N; i++) {
            r.emplace_back(std::vector<T>(N, 0));
            q.emplace_back(std::vector<T>(N, 0));
            transpq.emplace_back(std::vector<T>(N, 0));
            tmp100.emplace_back(std::vector<T>(N, 0));
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                r[i][j] = A[i][j];
            }
        }
        qr_decomposition(epsilon, q, r, N, flag);
        transpose(q, transpq, N);
        btmp = transpq * B;
        y = back_gauss(r, btmp);
        x = y;
        if (flag) {
            std::cout << "solve of QR " << std::endl;
            printvector(x);
        }
        tmp100 = q * r;
        if (flag) {
            std::cout << "relative error of QR" << std::endl;
            std::cout << norm(true, tmp100 - A, N) / norm(true, A, N) << std::endl;
        }
    }
    void Solve_slau(T epsilon, bool flag) {
        if (flag) {
            Gauss(epsilon, true);
        }
        else {
            QR_solve(epsilon, true);
        }
    }
    matrix<T> Inverse(T epsilon, std::string str) {
        matrix<T> Ainv;
        std::vector<T> y;
        Ainv.reserve(N);
        std::vector<T> e;
        std::vector<T> etmp;
        matrix<T> r, tmpr, q, transpq;
        r.reserve(N);
        q.reserve(N);
        transpq.reserve(N);
        for (int i = 0; i < N; i++) {
            Ainv.emplace_back(std::vector<T>(N, 0));
            r.emplace_back(std::vector<T>(N, 0));
            q.emplace_back(std::vector<T>(N, 0));
            transpq.emplace_back(std::vector<T>(N, 0));
            tmpr.emplace_back(std::vector<T>(N, 0));
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                r[i][j] = A[i][j];
            }
        }
        qr_decomposition(epsilon, q, r, N, false);
        transpose(q, transpq, N);
        std::cout << std::endl;
        /*Ainv = R^-1*QT */
        for (int k = 0; k < N; k++) {
            for (int ii = 0; ii < N; ii++) {
                for (int jj = 0; jj < N; jj++) {
                    tmpr[ii][jj] = r[ii][jj];
                }
            }
            etmp.resize(N);
            for (int i = 0; i < N; i++) {
                etmp[i] = 0;
            }
            etmp[k] = 1;
            /*etmp = transpq*e;*/
            y = back_gauss(tmpr, etmp);
            /*printvector(y);*/
            /*transpose(q, transpq, N);
            btmp = transpq * B;
            y = back_gauss(r, btmp);
            x = y;*/
            for (int i = 0; i < N; i++) {
                Ainv[i][k] = y[i];
            }
            //this->input(str);
        }
        Ainv = Ainv * transpq;
        return Ainv;
    }
    void checkresidual(T epsilon) {
        Gauss(epsilon, false);
        std::vector<T> y = this->x;
        std::vector<T> h = A * y;
        T tmp = evresidual(h, B);
        std::cout << std::endl;
        std::cout << "Residual " << tmp << std::endl;
    }
    T conditionality(bool flag) {
        T tmp1 = norm(flag, A, N);
        T tmp2 = norm(flag, Inverse(A), N);
        return tmp1 * tmp2;
    }
    void difB() {
        for (int j = 0; j < N; j++) {
            B[j] = 0.01 + B[j];
        }
    }
    void changeB(std::vector<T>& vec) {
        for (int i = 0; i < N; i++) {
            B[i] = vec[i];
        }
    }
    T specresidual(std::vector<T>& z) {
        std::vector<T> tmpvec1 = z;
        std::vector<T> tmpvec2 = A * tmpvec1;
        tmpvec1 = tmpvec2 - B;
        return vector_norm(tmpvec1);
    }
    void check_conditionality(T epsilon, std::string str) {
        std::vector<T> y;
        T t1, t2, y1, y2;
        t1 = vector_norm(this->B);
        //printvector(B);
        //std::cout << t1 << std::endl;
        Gauss(epsilon, false);
        y1 = vector_norm(this->x);
        //printvector(x);
        //std::cout << y1 <<std::endl;
        /*printvector(y);
        printvector(B);*/
        this->input(str);
        difB();
        t2 = vector_norm(this->B) - t1;
        //printvector(B);
        //std::cout << t2 << std::endl;
        Gauss(epsilon, false);
        y2 = vector_norm(this->x) - y1;
        //printvector(x);
        //std::cout << y2 << std::endl;
        std::cout << "assessment of conditionality: " << (y2 / y1) / (t2 / t1) << std::endl;
    };
    void check_multiply(T epsilon, std::string str) {
        matrix<T> Ainverse = Inverse(epsilon, str);
        matrix<T> R;
        R.reserve(N);
        for (int i = 0; i < N; i++) {
            R.emplace_back(std::vector<T>(N, 0));
        }
        R = Ainverse * A;
        epsilon = pow(10, -6);
        printmatrix(R, N, epsilon);
    }
    void evaluate_conditionality(bool flag, T epsilon, std::string str) {
        matrix<T> Ainverse = Inverse(epsilon, str);
        std::cout << "cond: " << norm(flag, A, N) * norm(flag, Ainverse, N) << std::endl;
    }
    void check_solve() {
        std::vector<T> u;
        u.resize(N);
        for (int i = 0; i < N; i++) {
            u[i] = 1 - x[i];
        }
        std::cout << std::endl;
        std::cout << "error of solve" << std::endl;
        std::cout << vector_norm(u) << std::endl;
    }
    void simple_iteration(T eps, T tau, std::vector<T>& x0, bool flag) {
        int kest;
        std::vector<T> xcur(N, 0);
        std::vector<T> tmpv;
        std::vector<T> y;
        std::vector<T> difx;
        std::vector<T> bufv(N, 0);
        matrix<T> BUF;
        matrix<T> E;
        matrix<T> C;
        std::vector<T> x1;
        T eps0 = eps / 10;
        T normxp;
        E.reserve(N);
        C.reserve(N);
        difx.reserve(N);
        BUF.reserve(N);
        for (int i = 0; i < N; i++) {
            E.emplace_back(std::vector<T>(N, 0));
            C.emplace_back(std::vector<T>(N, 0));
            BUF.emplace_back(std::vector<T>(N, 0));
        }
        T normA = norm(flag, A, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                E[i][j] = 0;
                BUF[i][j] = A[i][j];
            }
            E[i][i] = 1;
        }
        y = B;
        for (int i = 0; i < N; i++) {
            if (BUF[i][i] < 0) {
                for (int j = 0; j < N; j++) {
                    BUF[i][j] = -BUF[i][j];
                }
                y[i] = -y[i];
            }
        }
        if (flag) {
            double maxA = fabs(BUF[0][0]);
            for (int i = 1; i < N; i++) {
                if (maxA < fabs(BUF[i][i])) {
                    maxA = fabs(BUF[i][i]);
                }
            }
            //std::cout << maxA << std::endl;
            tau = 1.0 / maxA;
        }
        else {
            tau = BUF[0][0];
            for (int i = 0; i < N; i++) {
                double sum = BUF[i][i];
                for (int j = 0; j < N; j++) {
                    if (i != j)
                        sum += fabs(BUF[i][j]);
                }
                if (1 / sum < tau)
                    tau = 1. / sum;
            }
        }//tau = 0.0329702;
        BUF = tau * BUF;
        T r0;
        /*else if (way == 1) {
            tau = BUF[0][0];
            for (int i = 0; i < N; i++) {
                double sum = BUF[i][i];
                for (int j = 0; j < N; j++) {
                    if (i != j)
                        sum += fabs(C[i][j]);
                }
                if (1 / sum < tau)
                    tau = 1. / sum;
            }
        }
        else {
            tau *= 0.9;
            tau = 1. / norm(flag, C, N);
        }*/
        C = E - BUF;
        tmpv = B;
        T normC = norm(flag, C, N);
        y = tau * y;
        x1 = x0;
        //T rk = specresidual();
        T normv = 10;
        normxp = (1 - normC) / normC;
        //std::cout << normC << std::endl;
        int iter_count = 0;
        //std::cout << tau << std::endl;
        std::cout << std::endl;
        std::cout << "first" << std::endl;
        std::cout << tau << std::endl;
        while (normv > (1 - normC) * eps / normC) {
            iter_count = iter_count + 1;
            tmpv = x1;
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x1[j];
                }
            }
            x1 = bufv + y;
            difx = x1 - tmpv;
            normv = vector_norm(difx);
            if (iter_count == 1) {
                kest = round(log((1 - normC) * eps / normv) / log(normC));
            }
        }
        std::cout << "kest " << kest << std::endl;
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of simple iterations" << std::endl;
        printvector(x1);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
        std::cout << "second" << std::endl;
        x1 = x0;
        normv = 10;
        normxp = 1;
        iter_count = 0;
        std::cout << tau << std::endl;
        while (normv > eps) {
            iter_count = iter_count + 1;
            tmpv = x1;
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x1[j];
                }
            }
            x1 = bufv + y;
            difx = x1 - tmpv;
            normv = vector_norm(difx);
        }
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of simple iterations" << std::endl;
        printvector(x1);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
        /* x = x0;
         normv = 10;
         normxp = 1;
         iter_count = 0;
         while (normv > normxp * eps + eps0) {
             iter_count = iter_count + 1;
             tmpv = x;
             normxp = vector_norm(tmpv);
             for (int i = 0; i < N; i++) {
                 bufv[i] = 0;
                 for (int j = 0; j < N; j++) {
                     bufv[i] += C[i][j] * x[j];
                 }
             }
             x = bufv + y;
             difx = x - tmpv;
             normv = vector_norm(difx);
         }
         std::cout << "Norm of C " << normC << std::endl;
         std::cout << "Iterations " << iter_count << std::endl;
         std::cout << "Solve of simple iterations" << std::endl;
         printvector(x);
         std::cout << "residual " << specresidual() << std::endl;
         std::cout << std::endl;*/
        std::cout << "third" << std::endl;
        x1 = x0;
        normv = 10;
        normxp = 10;
        iter_count = 0;
        while (normv > eps) {
            iter_count = iter_count + 1;
            tmpv = x1;
            normxp = vector_norm(tmpv);
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x1[j];
                }
            }
            x1 = bufv + y;
            difx = x1 - tmpv;
            normv = specresidual(x1);
            //normv = vector_norm(difx);
        }
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of simple iterations" << std::endl;
        printvector(x1);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
    }
    void Jacoby(T eps, T tau, std::vector<T>& x0, bool flag) {
        int kest;
        std::vector<T> tmpv(N, 0);
        std::vector<T> y(N, 0);
        std::vector<T> difx(N, 0);
        std::vector<T> bufv(N, 0);
        T eps0 = 0.0001;
        matrix<T> C;
        std::vector<T> x1;
        T t;
        C.reserve(N);
        for (int i = 0; i < N; i++) {
            C.emplace_back(std::vector<T>(N, 0));
        }
        for (int i = 0; i < N; i++) {
            t = A[i][i];
            for (int j = 0; j < N; j++) {
                C[i][j] = -A[i][j] / t;
            }
            C[i][i] = 0;
            y[i] = B[i] / t;
        }
        x1 = x0;
        T normv = 10;
        T normxp = 1;
        T normC = norm(flag, C, N) / 2;
        int iter_count = 0;
        //T rk = specresidual();
        std::cout << "first" << std::endl;
        while (normv > (1 - normC) * eps / normC) {
            iter_count = iter_count + 1;
            tmpv = x1;
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x1[j];
                }
            }
            x1 = bufv + y;
            difx = x1 - tmpv;
            normv = vector_norm(difx);
            if (iter_count == 1) {
                kest = round(log((1 - normC) * eps / normv) / log(normC));
            }
        }
        std::cout << "kest " << kest << std::endl;
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of Jacoby" << std::endl;
        printvector(x1);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
        std::cout << "second" << std::endl;
        x1 = x0;
        normv = 10;
        normxp = 1;
        iter_count = 0;
        while (normv > eps) {
            iter_count = iter_count + 1;
            tmpv = x1;
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x1[j];
                }
            }
            x1 = bufv + y;
            difx = x1 - tmpv;
            normv = vector_norm(difx);
        }
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of Jacoby" << std::endl;
        printvector(x1);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
        /*std::cout << "second" << std::endl;
        normv = 10;
        normxp = 1;
        iter_count = 0;
        x = x0;
        while (normv > normxp * eps + eps0) {
            iter_count = iter_count + 1;
            tmpv = x;
            normxp = vector_norm(tmpv);
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x[j];
                }
            }
            x = bufv + y;
            difx = x - tmpv;
            normv = vector_norm(difx);
        }
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of Jacoby" << std::endl;
        printvector(x);
        std::cout << "residual " << specresidual() << std::endl;
        std::cout << std::endl;*/
        std::cout << "third" << std::endl;
        normv = 10;
        normxp = 1;
        iter_count = 0;
        x1 = x0;
        while (normv > eps) {
            iter_count = iter_count + 1;
            tmpv = x1;
            normxp = vector_norm(tmpv);
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x1[j];
                }
            }
            x1 = bufv + y;
            difx = x1 - tmpv;
            normv = specresidual(x1);
        }
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of Jacoby" << std::endl;
        printvector(x);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
    }
    void Zeydel_relax(T omega, T eps, std::vector<T>& x0, bool flag) {
        std::vector<T> tmpv(N, 0);
        std::vector<T> y(N, 0);
        std::vector<T> difx(N, 0);
        std::vector<T> bufx(N, 0);
        int iter_count = 0;
        T sum;
        T sum1;
        T sum2;
        T buf;
        T nC1 = 0;
        T nC2 = 0;
        T nC = 0;
        T q;
        T r0;
        std::vector<T> x1;
        int kest;
        for (int i = 0; i < N; i++) {
            sum1 = 0;
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j];
            }
            sum2 = 0;
            for (int j = i + 1; j < N; j++) {
                sum2 += A[i][j];
            }
            buf = A[i][i];
            sum1 = sum1 / buf;
            sum2 = sum2 / buf;
            sum = sum1 + sum2;
            if (sum1 > nC1)
                nC1 = sum1;

            if (sum2 > nC2)
                nC2 = sum2;

            if (sum > nC)
                nC = sum;
        }
        nC1 *= omega;
        q = nC2 / (1 - nC1);
        nC2 *= omega;
        nC2 += abs(1 - omega);
        T tmpcof;
        if (flag) {
            tmpcof = nC1 / (1 - nC2);
        }
        else {
            tmpcof = (1 - nC) / nC;
        }
        nC *= omega;
        nC += abs(1 - omega);
        x1 = x0;
        T normv = 10;
        //T normC = norm(flag, C, N);
        /*std::cout << std::endl;
        std::cout << "first" << std::endl;
        while (normv > (1-normC)*eps/normC) {
            iter_count = iter_count + 1;
            tmpv = x;
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x[j];
                }
            }
            x = bufv + y;
            difx = x - tmpv;
            normv = vector_norm(difx);
        }
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of simple iterations" << std::endl;
        printvector(x);
        std::cout << "residual " << specresidual() << std::endl;
        std::cout << std::endl;
        std::cout << "second" << std::endl;
        x = x0;
        normv = 10;
        normxp = 1;
        iter_count = 0;
        while (normv > normxp*eps+eps0) {
            iter_count = iter_count + 1;
            tmpv = x;
            normxp = vector_norm(tmpv);
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x[j];
                }
            }
            x = bufv + y;
            difx = x - tmpv;
            normv = vector_norm(difx);
        }
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of simple iterations" << std::endl;
        printvector(x);
        std::cout << "residual " << specresidual() << std::endl;
        std::cout << std::endl;
        std::cout << "third" << std::endl;
        x = x0;
        normv = 10;
        normxp = 1;
        iter_count = 0;
        while (normv/(normxp + eps0) > eps) {
            iter_count = iter_count + 1;
            tmpv = x;
            normxp = vector_norm(tmpv);
            for (int i = 0; i < N; i++) {
                bufv[i] = 0;
                for (int j = 0; j < N; j++) {
                    bufv[i] += C[i][j] * x[j];
                }
            }
            x = bufv + y;
            difx = x - tmpv;
            normv = vector_norm(difx);
        }
        std::cout << "Norm of C " << normC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of simple iterations" << std::endl;
        printvector(x);
        std::cout << "residual " << specresidual() << std::endl;
        std::cout << std::endl;*/
        T eps0 = 0.000001;
        std::cout << "first" << std::endl;
        while (normv > tmpcof * eps) {
            iter_count = iter_count + 1;
            sum = 0;
            for (int k = 1; k < N; k++) {
                sum = sum + A[0][k] * x1[k] / A[0][0];
            }
            sum = omega * sum;
            bufx[0] = (1 - omega) * x1[0] - sum + omega * B[0] / A[0][0];
            for (int i = 1; i < N - 1; i++) {
                sum1 = 0;
                sum2 = 0;
                for (int j = 0; j <= i - 1; j++) {
                    sum1 = sum1 + A[i][j] * bufx[j] / A[i][i];
                }
                sum1 = -1 * omega * sum1;
                for (int j = i + 1; j < N; j++) {
                    sum2 = sum2 + A[i][j] * x1[j] / A[i][i];
                }
                sum2 = -1 * omega * sum2;
                bufx[i] = sum1 + sum2 + (1 - omega) * x1[i] + omega * B[i] / A[i][i];
            }
            sum = 0;
            for (int k = 0; k < N - 1; k++) {
                sum = sum + A[N - 1][k] * bufx[k] / A[N - 1][N - 1];
            }
            sum = omega * sum;
            bufx[N - 1] = (1 - omega) * x1[N - 1] - sum + omega * B[N - 1] / A[N - 1][N - 1];
            difx = bufx - x1;
            normv = vector_norm(difx);
            x1 = bufx;
            if (iter_count == 1) {
                kest = log((1 - q) * eps / normv) / log(q);
            }
        }
        std::cout << "kest " << kest << std::endl;
        std::cout << "Norm of C " << nC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of Relaxation" << " with omega " << omega << std::endl;
        printvector(x1);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
        x1 = x0;
        normv = 10;
        T normxp = 1;
        iter_count = 0;
        std::cout << "second" << std::endl;
        while (normv > eps) {
            iter_count = iter_count + 1;
            sum = 0;
            for (int k = 1; k < N; k++) {
                sum = sum + A[0][k] * x1[k] / A[0][0];
            }
            sum = omega * sum;
            bufx[0] = (1 - omega) * x1[0] - sum + omega * B[0] / A[0][0];
            for (int i = 1; i < N - 1; i++) {
                sum1 = 0;
                sum2 = 0;
                for (int j = 0; j <= i - 1; j++) {
                    sum1 = sum1 + A[i][j] * bufx[j] / A[i][i];
                }
                sum1 = -1 * omega * sum1;
                for (int j = i + 1; j < N; j++) {
                    sum2 = sum2 + A[i][j] * x1[j] / A[i][i];
                }
                sum2 = -1 * omega * sum2;
                bufx[i] = sum1 + sum2 + (1 - omega) * x1[i] + omega * B[i] / A[i][i];
            }
            sum = 0;
            for (int k = 0; k < N - 1; k++) {
                sum = sum + A[N - 1][k] * bufx[k] / A[N - 1][N - 1];
            }
            sum = omega * sum;
            bufx[N - 1] = (1 - omega) * x1[N - 1] - sum + omega * B[N - 1] / A[N - 1][N - 1];
            difx = bufx - x1;
            normv = vector_norm(difx);
            x1 = bufx;
        }
        std::cout << "Norm of C " << nC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of Relaxation" << " with omega " << omega << std::endl;
        printvector(x);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
        /* std::cout << "second" << std::endl;
         x = x0;
         normv = 10;
         normxp = 1;
         iter_count = 0;
         while (normv > normxp * eps + eps0) {
             iter_count = iter_count + 1;
             sum = 0;
             for (int k = 1; k < N; k++) {
                 sum = sum + A[0][k] * x[k] / A[0][0];
             }
             sum = omega * sum;
             bufx[0] = (1 - omega) * x[0] - sum + omega * B[0] / A[0][0];
             for (int i = 1; i < N - 1; i++) {
                 sum1 = 0;
                 sum2 = 0;
                 for (int j = 0; j <= i - 1; j++) {
                     sum1 = sum1 + A[i][j] * bufx[j] / A[i][i];
                 }
                 sum1 = -1 * omega * sum1;
                 for (int j = i + 1; j < N; j++) {
                     sum2 = sum2 + A[i][j] * x[j] / A[i][i];
                 }
                 sum2 = -1 * omega * sum2;
                 bufx[i] = sum1 + sum2 + (1 - omega) * x[i] + omega * B[i] / A[i][i];
             }
             sum = 0;
             for (int k = 0; k < N - 1; k++) {
                 sum = sum + A[N - 1][k] * bufx[k] / A[N - 1][N - 1];
             }
             sum = omega * sum;
             bufx[N - 1] = (1 - omega) * x[N - 1] - sum + omega * B[N - 1] / A[N - 1][N - 1];
             normxp = vector_norm(bufx);
             difx = bufx - x;
             normv = vector_norm(difx);
             x = bufx;
         }
         std::cout << "Norm of C " << nC << std::endl;
         std::cout << "Iterations " << iter_count << std::endl;
         std::cout << "Solve of Relaxation" << " with omega " << omega << std::endl;
         printvector(x);
         std::cout << "residual " << specresidual() << std::endl;
         std::cout << std::endl;*/
        x1 = x0;
        normv = 10;
        normxp = 1;
        iter_count = 0;
        std::cout << "third" << std::endl;
        while (normv > eps) {
            iter_count = iter_count + 1;
            sum = 0;
            for (int k = 1; k < N; k++) {
                sum = sum + A[0][k] * x1[k] / A[0][0];
            }
            sum = omega * sum;
            bufx[0] = (1 - omega) * x1[0] - sum + omega * B[0] / A[0][0];
            for (int i = 1; i < N - 1; i++) {
                sum1 = 0;
                sum2 = 0;
                for (int j = 0; j <= i - 1; j++) {
                    sum1 = sum1 + A[i][j] * bufx[j] / A[i][i];
                }
                sum1 = -1 * omega * sum1;
                for (int j = i + 1; j < N; j++) {
                    sum2 = sum2 + A[i][j] * x1[j] / A[i][i];
                }
                sum2 = -1 * omega * sum2;
                bufx[i] = sum1 + sum2 + (1 - omega) * x1[i] + omega * B[i] / A[i][i];
            }
            sum = 0;
            for (int k = 0; k < N - 1; k++) {
                sum = sum + A[N - 1][k] * bufx[k] / A[N - 1][N - 1];
            }
            sum = omega * sum;
            bufx[N - 1] = (1 - omega) * x1[N - 1] - sum + omega * B[N - 1] / A[N - 1][N - 1];
            normxp = vector_norm(bufx);
            difx = bufx - x1;
            //vector_norm(difx);
            x1 = bufx;
            normv = specresidual(x1);
        }
        std::cout << "Norm of C " << nC << std::endl;
        std::cout << "Iterations " << iter_count << std::endl;
        std::cout << "Solve of Relaxation" << " with omega " << omega << std::endl;
        printvector(x1);
        std::cout << "residual " << specresidual(x1) << std::endl;
        std::cout << std::endl;
    }
    void vectinput(std::vector<T>& a, std::vector<T>& b, std::vector<T>& c, std::vector<T>& d, int M) {
        d[0] = 6;
        for (int i = 0; i < M; i++) {
            a[i] = 1;
            c[i] = 1;
            b[i] = 4;
            d[i] = 10 - 2 * (i % 2);
        }
        d[M - 1] = 9 - 3 * (M % 2);
    }
    std::vector<T> matrixmult(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& vec, int M) {
        std::vector<T> result;
        result.resize(M);
        result[0] = b[0] * vec[0] + c[0] * vec[1];
        for (int i = 1; i < M - 1; ++i) {
            result[i] = a[i] * vec[i - 1] + b[i] * vec[i] + c[i] * vec[i + 1];
        }
        result[M - 1] = a[M - 1] * vec[M - 2] + b[M - 1] * vec[M - 1];
        return result;
    }
    T checkres(std::vector<T>& a, std::vector<T>& b, std::vector<T>& c, std::vector<T>& d, std::vector<T>& r, int M) {
        T res;
        std::vector<T> resvec(M, 0);
        resvec = matrixmult(a, b, c, r, M);
        resvec = resvec - d;
        res = vector_norm(resvec);
        return res;
    }
    /* void Threediag_relax(T eps, std::vector<T>& a, std::vector<T>& b, std::vector<T>& c, std::vector<T>& d, bool flag) {
         std::vector<T> alpha(N, 0);
         std::vector<T> betha(N, 0);
         vectinput();
         T tmp = b[0];
         alpha[1] = c[0] / tmp;
         betha[1] = d[0] / tmp;
         for (int i = 1; i < N - 1; i++) {
             tmp = b[i] - a[i] * alpha[i];
             alpha[i + 1] = c[i] / tmp;
             betha[i + 1] = (d[i] + a[i] * betha[i]) / tmp;
         }
         x[N - 1] = (d[N - 1] + a[N - 1] * betha[N - 1]) / (b[N - 1] - a[N - 1] * alpha[N - 1]);
         for (int i = N - 2; i >= 0; --i) {
             x[i] = alpha[i + 1] * x[i + 1] + betha[i + 1];
         }
         printvector(x);
     }*/
    void Threediag_relaxation(T eps, T omega, int M) {
        int iter = 0;
        std::vector<T> a(M, 0);
        std::vector<T> b(M, 0);
        std::vector<T> c(M, 0);
        std::vector<T> d(M, 0);
        std::vector<T> alpha(M, 0);
        std::vector<T> betha(M, 0);
        std::vector<T> xcur(M, 0);
        std::vector<T> xprev(M, 1);
        std::vector<T> xdif(M, 0);
        T p, q, r;
        vectinput(a, b, c, d, M);
        xdif = xcur - xprev;
        while (vector_norm(xdif) >= eps) {
            iter = iter + 1;
            xprev = xcur;
            xcur[0] = -c[0] / b[0] * omega * xprev[1] + d[0] * omega / b[0] - xprev[0] * (omega - 1);
            for (int i = 1; i < M - 1; i++) {
                p = -a[i] / b[i] * omega * xcur[i - 1];
                q = -c[i] / b[i] * omega * xprev[i + 1];
                r = d[i] * omega / b[i];
                xcur[i] = p + q + r - xprev[i] * (omega - 1);
            }
            xcur[M - 1] = -(a[M - 1] / b[M - 1]) * omega * xcur[M - 2] + d[M - 1] * omega / b[M - 1] - xprev[M - 1] * (omega - 1);
            xdif = xcur - xprev;
        }
        std::cout << "Iterations " << iter << " " << "residual with omega " << omega << " " << checkres(a, b, c, d, xcur, M) << std::endl;
    }
};
int main() {
    //"D:\\input4_laba1.txt" D:\\input3_laba1.txt"
    Hsole<double> C;
    double epsilon = pow(10, -4);
    std::string str_input = "D:\\input2_laba2.txt";
    std::string str_output = "D:\\output_laba2.txt";
    std::vector<double> x0 = { 0, 0, 0, 0 };
    std::vector<double> x11 = { 5, 4, 2, 3 };
    std::vector<double> x2 = { 5, 1, -20, 1 };
    C.input(str_input);
    C.simple_iteration(epsilon, 0.000166013, x11, true);
    //C.input(str_input);
    std::cout << "AAA" << std::endl;
    C.simple_iteration(epsilon, 0.000166013, x11, false);
    //C.input(str_input);
    C.Jacoby(epsilon, 0.5, x11, false); //0.5
    //C.input(str_input);
    C.Zeydel_relax(1, epsilon, x2, true);
    //C.input(str_input);
    C.Zeydel_relax(1.5, epsilon, x0, true);
    //C.input(str_input);
    C.Zeydel_relax(0.7, epsilon, x0, true);
    C.Threediag_relaxation(epsilon, 1, 211);
    C.Threediag_relaxation(epsilon, 0.5, 211);
    C.Threediag_relaxation(epsilon, 1.5, 211);
}
