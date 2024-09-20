// ConsoleApplication87.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
template <class T> using matrix = std::vector<std::vector<T>>;
//template<typename T>
//void printmatrix(matrix<T> &A, int n) {
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            std::cout << A[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//}
template<typename T>
void printmatrix(matrix<T>& A, int n, T epsilon) {
    epsilon = epsilon * 10;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (abs(A[i][j]) < epsilon) {
                std::cout << 0 << " ";
            }
            else {
                std::cout << A[i][j] << " ";
            }
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
template <typename T>
matrix<T> operator*(const T& lhs, const matrix<T>& A) {
    matrix<T> B;
    int m = A.size();
    B.reserve(m);
    for (int i = 0; i < m; i++) {
        B.emplace_back(std::vector<T>(m, 0));
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            B[i][j] = lhs * A[i][j];
        }
    }
    return B;
};
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
    result.reserve(rhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        result.emplace_back(std::vector<T>(rhs.size(), 0));
    }
    for (int i = 0; i < rhs.size(); i++) {
        for (int j = 0; j < rhs.size(); j++) {
            result[i][j] = 0;
            for (int k = 0; k < rhs.size(); k++) {
                result[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }
    return result;
};
template<typename T>
T operator*(std::vector<T>& lhs, std::vector <T>& rhs) {
    T result = 0;
    for (int i = 0; i < lhs.size(); i++) {
        result += lhs[i] * rhs[i];
    }
    return result;
};
template<typename T>
std::vector<T> operator- (std::vector<T>& lhs, std::vector<T> rhs) {
    int m = lhs.size();
    std::vector<T> result(m, 0);
    for (int i = 0; i < m; i++) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}
template<typename T>
std::vector<T> operator+ (std::vector<T>& lhs, std::vector<T> rhs) {
    int m = lhs.size();
    std::vector<T> result(m, 0);
    for (int i = 0; i < m; i++) {
        result[i] = lhs[i] + rhs[i];
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
            result[i][j] = 0;
            result[i][j] = lhs[i][j] - rhs[i][j];
        }
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
T residual(std::vector<T>& a, std::vector<T>& b) {
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
void specmult(matrix<T>& A, matrix<T>& B) {
    matrix<T> C;
    int N = A.size();
    C.reserve(N);
    for (int i = 0; i < N; i++) {
        C.emplace_back(std::vector<T>(N, 0));
    }
    A = B * A;
    transpose(B, C, N);
    A = A * C;
}
template <typename T>
bool my_eq(const std::vector<T>& x, const std::vector<T>& y, const T eps)
{
    bool t = true;
    int n = x.size();
    if (x[0] * y[0] < 0) {
        for (int i = 0; i < n; i++) {
            t *= abs(x[i] + y[i]) < eps;
        }
    }
    else {
        for (size_t i = 0; i < n; ++i) {
            t *= abs(x[i] - y[i]) < eps;
        }
    }
    return t;
};
template <typename T>
T eqiual(std::vector<T>& a, std::vector<T>& b) {
    int n = a.size();
    T tmp;
    T max = 0;
    for (int i = 0; i < n; i++) {
        tmp = abs(abs(a[i]) - abs(b[i]));
        if (tmp > max) {
            max = tmp;
        }
    }
    return max;
}
template<typename T>
void transposemultiply(matrix<T>& A, int k, int l, int M, T a, T b, int &buf) {
    T tmp;
    for (int i = k; i < M; i++) {
        tmp = A[k+1][i];
        A[k+1][i] = a * A[k+1][i] + b * A[l][i];
        A[l][i] = -b * tmp + a * A[l][i];
        buf += 4;
    }
    for (int i = k; i < M; i++) {
        tmp = A[i][k+1];
        A[i][k+1] = a * A[i][k+1] + b * A[i][l];
        A[i][l] = -b * tmp + a * A[i][l];
        buf += 4;
    }
}
template<typename T>
void transposemultiply(matrix<T>& A, int k, int l, int M, T a, T b) {
    T tmp;
    for (int i = k; i < M; i++) {
        tmp = A[k + 1][i];
        A[k + 1][i] = a * A[k + 1][i] + b * A[l][i];
        A[l][i] = -b * tmp + a * A[l][i];
    }
    for (int i = k; i < M; i++) {
        tmp = A[i][k + 1];
        A[i][k + 1] = a * A[i][k + 1] + b * A[i][l];
        A[i][l] = -b * tmp + a * A[i][l];
    }
}
template <typename T>
class Eigenvecandval { //Heterogeneous system of linear equations
    matrix<T> A;
    int N;
    std::vector<T> b;
    std::vector<T> x;
public:
    void input(std::string str) { //аргумент строка
        std::ifstream in;
        T tmp;
        in.open(str); //input_laba1.txt"
        in >> N;
        if (in.is_open()) {
            A.reserve(N);
            for (int i = 0; i < N; i++) {
                A.emplace_back(std::vector<T>(N, 0));
            }
            //B.resize(N);
            if (in.is_open()) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        in >> tmp;
                        A[i][j] = tmp;
                    }
                }
            }
            in.close();
        }
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
                //out << B[i] << " ";
                //out << x[i] << std::endl;
            }
        }
        out.close();
    }
    T determinant(matrix<T>& A, T lambda) {
        matrix<T> B = A;
        T p = 1/pow(10, 10);
        B = p * A;
        lambda = lambda * p;
        T determ = 1.0;
        T d;
        for (int i = 0; i < N; i++) {
            B[i][i] -= lambda;
        }
        for (int k = 0; k < N; k++) {
            for (int i = k+1; i < N; i++) {
                d = B[i][k] / B[k][k];
                for (int j = k; j < N; j++) {
                    B[i][j] -= d * B[k][j];
                }
                //v[i] -= d * v[k];
            }
        }
        for (int k = 0; k < N; k++) {
            determ *= B[k][k];
        }
        return determ;
    }
    void spec_input(int M) {
        A.reserve(M);
        for (int i = 0; i < M; i++) {
            A.emplace_back(std::vector<T>(M, 0));
        }
        T d = pow(10, 10);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                A[i][j] = d / (i + j + 1);
            }
        }
        N = M;
        printmatrix(A, M, 0.0000001);
    }
    std::vector<T> Gauss(matrix<T>& F, std::vector<T>& v, T epsilon, bool flag) { //параметр эпсилон
        std::vector<T> y;
        epsilon = epsilon;
        T max_element;
        T d;
        int max_number;
        for (int k = 0; k < N; k++) { //можно добавить флажок на полный частичный выбор
            max_element = abs(F[k][k]);
            max_number = k;
            for (int i = k + 1; i < N; i++) {
                if (abs(F[i][k]) > max_element) {
                    max_element = abs(F[i][k]);
                    max_number = i;
                }
            }
            //if (max_element < epsilon) {
            //    std::cout << "BAD MATRIX" << std::endl; //добавить про строки
            //    std::cout << "row " << max_number << "cow " << max_number << std::endl;
            //    exit(1);
            //}
            if (k != max_number) {
                for (int i = 0; i < N; i++) {
                    swap(F[k][i], F[max_number][i]);
                    /*tmp = A[k][i];
                    A[k][i] = A[max_number][i];
                    A[max_number][i] = tmp;*/
                }
                swap(v[k], v[max_number]);
                /*tmp = B[k];
                    B[k] = B[max_number];
                    B[max_number] = tmp;*/
            }
            for (int i = k + 1; i < N; i++) {
                d = F[i][k] / F[k][k];
                for (int j = k; j < N; j++) {
                    F[i][j] -= d * F[k][j];
                }
                v[i] -= d * v[k];
            }
        }
        y = back_gauss(F, v);
        if (flag) {
            std::cout << "solve of Gauss " << std::endl;
            printvector(y);
        }
        return y;
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
        /*printvector(y);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    std::cout << A[i][j] << " ";
                }
                std::cout << B[i] << std::endl;
            }*/
    }
    matrix<T> qr_decomposition(T epsilon, matrix<T>& r, int n) {
        //epsilon = epsilon * 10000;
        matrix<T> t;
        matrix<T> q;
        t.reserve(n);
        q.reserve(n);
        for (int i = 0; i < n; i++) {
            t.emplace_back(std::vector<T>(n, 0));
            q.emplace_back(std::vector<T>(n, 0));
        }
        //std::cout << "qr " << std::endl;
        for (int i = 0; i < n; i++) {
            t[i][i] = 1;
        }
        if (n == 1) {
            return r;
        }
        for (int i = 0; i < n-1; i++) {
            for (int j = i + 1; j < n; j++) {
                //std::cout << n << " qr " << i << " " << j << std::endl;
                T sqr = sqrt(r[i][i] * r[i][i] + r[j][i] * r[j][i]);
                if (sqr < epsilon) {
                    continue;
                };
                T c = r[i][i] / sqr;
                T s = r[j][i] / sqr;
                //std::cout << "AAA" << std::endl;
                for (int k = 0; k < n; k++) {
                    T buf = t[i][k];
                    t[i][k] = c * buf + s * t[j][k];
                    t[j][k] = -1 * s * buf + c * t[j][k];
                };
                //std::cout << "BBB" << std::endl;
                for (int k = i; k < n; k++) {
                    T buf = r[i][k];
                    r[i][k] = c * buf + s * r[j][k];
                    if (abs(r[i][k]) < epsilon)
                        r[i][k] = 0;
                    r[j][k] = -1 * s * buf + c * r[j][k];
                    if (abs(r[j][k]) < epsilon)
                        r[j][k] = 0;
                };
                //std::cout << "CCC" << std::endl;
            }
        }
        /*if (abs(r[n - 1][n - 1]) < epsilon) {
            std::cout << "BAD MATRIX" << std::endl;
            exit(1);
        }*/
        transpose(t, q, n);
        //std::cout << "DDD" << std::endl;
        r = r * q;
        //std::cout << "EEE" << std::endl;
        /*if (flag) {
            std::cout << "Q" << std::endl;
            printmatrix(q, N, epsilon);
            std::cout << std::endl;
            std::cout << "R" << std::endl;
            printmatrix(r, N, epsilon / 100);
            std::cout << std::endl;
        }*/
        return r;
        //std::vector<std::vector<T>> test_A;
        //for (int i = 0; i < N; i++) {
        //    for (int j = 0; j < N; j++) {
        //        test_A = 
        //    }
        //}
        ////std::vector<std::vector<double>> test_A = { {0, 0}, {0, 0} };
        //for (int i = 0; i < a.size(); i++) {
        //    for (int j = 0; j < a.size(); j++) {
        //        for (int inner = 0; inner < a.size(); inner++)
        //            test_A[i][j] += q[i][inner] * r[inner][j];
        //    }
        //}
    }
    matrix<T> qr_decomposition(T epsilon, matrix<T>& r, int n, int &vtmp) {
        //epsilon = epsilon * 10000;
        matrix<T> t;
        matrix<T> q;
        t.reserve(n);
        q.reserve(n);
        for (int i = 0; i < n; i++) {
            t.emplace_back(std::vector<T>(n, 0));
            q.emplace_back(std::vector<T>(n, 0));
        }
        //std::cout << "qr " << std::endl;
        for (int i = 0; i < n; i++) {
            t[i][i] = 1;
        }
        if (n == 1) {
            return r;
        }
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                //std::cout << n << " qr " << i << " " << j << std::endl;
                if ((abs(r[i + 1][i]) < epsilon) && (abs(r[j][i]) < epsilon))
                {
                    T c = 1;
                    T s = 0;
                    continue;
                }
                T sqr = sqrt(r[i][i] * r[i][i] + r[j][i] * r[j][i]);
                vtmp += 2;
                if (sqr < epsilon) {
                    continue;
                };
                T c = r[i][i] / sqr;
                T s = r[j][i] / sqr;
                vtmp += 2;
                //std::cout << "AAA" << std::endl;
                for (int k = 0; k < n; k++) {
                    T buf = t[i][k];
                    t[i][k] = c * buf + s * t[j][k];
                    t[j][k] = -1 * s * buf + c * t[j][k];
                    vtmp += 4;
                };
                //std::cout << "BBB" << std::endl;
                for (int k = i; k < n; k++) {
                    T buf = r[i][k];
                    r[i][k] = c * buf + s * r[j][k];
                    if (abs(r[i][k]) < epsilon) {
                        r[i][k] = 0;
                    }
                    else {
                        vtmp += 2;
                    }
                    r[j][k] = -1 * s * buf + c * r[j][k];
                    if (abs(r[j][k]) < epsilon) {
                        r[j][k] = 0;
                    }
                    else {
                        vtmp += 2;
                    }
                };
                //std::cout << "CCC" << std::endl;
            }
        }
        /*if (abs(r[n - 1][n - 1]) < epsilon) {
            std::cout << "BAD MATRIX" << std::endl;
            exit(1);
        }*/
        transpose(t, q, n);
        //std::cout << "DDD" << std::endl;
        r = r * q;
        //std::cout << "EEE" << std::endl;
        /*if (flag) {
            std::cout << "Q" << std::endl;
            printmatrix(q, N, epsilon);
            std::cout << std::endl;
            std::cout << "R" << std::endl;
            printmatrix(r, N, epsilon / 100);
            std::cout << std::endl;
        }*/
        return r;
        //std::vector<std::vector<T>> test_A;
        //for (int i = 0; i < N; i++) {
        //    for (int j = 0; j < N; j++) {
        //        test_A = 
        //    }
        //}
        ////std::vector<std::vector<double>> test_A = { {0, 0}, {0, 0} };
        //for (int i = 0; i < a.size(); i++) {
        //    for (int j = 0; j < a.size(); j++) {
        //        for (int inner = 0; inner < a.size(); inner++)
        //            test_A[i][j] += q[i][inner] * r[inner][j];
        //    }
        //}
    }
    void qr_decomposition(T epsilon, matrix<T>& q, matrix<T>& r, int n, bool flag) {
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
                if (sqr < epsilon) {
                    continue;
                }
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
        /*if (abs(r[n - 1][n - 1]) < epsilon*epsilon) {
            std::cout << "BAD MATRIX" << std::endl;
            exit(1);
        }*/
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
    matrix<T> Inverse(matrix<T>& F, T epsilon, std::string str) {
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
                r[i][j] = F[i][j];
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
            /*printvector(y);
            transpose(q, transpq, N);
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
    void check_multiply(T epsilon, std::string str) {
        matrix<T> Ainverse = Inverse(epsilon, str);
        matrix<T> R;
        R.reserve(N);
        for (int i = 0; i < N; i++) {
            R.emplace_back(std::vector<T>(N, 0));
        }
        R = Ainverse * A;
        epsilon = pow(10, -6);
        std::cout << "multiply of matrix and inverse matrix" << std::endl;
        printmatrix(R, N, epsilon);
    }
    void QR_method(T eps, bool flag1, bool flag2, bool flag5) { //flag1 - Хессенберг, flag2 - cдвиг, flag3 - критерий останова 
        //flag3 - true - корень, false - сумма модулей.
        matrix<T> B = A;
        //matrix<T> E;
        std::vector<T> eigenvec(N, 0);
        T eps1 = eps;
        int mults = 0;
        /*E.reserve(N);
        for (int i = 0; i < N; i++) {
            E.emplace_back(std::vector<T>(N, 0));
        }*/
        double d;
        double alpha;
        double betha;
        if (flag1) {
            for (int i = 0; i < N - 2; i++) {
                for (int j = i + 2; j < N; j++) {
                    if (abs(B[i + 1][i]) < eps && abs(B[j][i]) < eps)
                    {
                        continue;
                    }
                    //else {
                    mults += 4;
                    d = sqrt(B[i + 1][i] * B[i + 1][i] + B[j][i] * B[j][i]);
                    alpha = B[i + 1][i] / d;
                    betha = B[j][i] / d;
                    transposemultiply(B, i, j, N, alpha, betha, mults);
                    //mults += 8 * (N - i);
                //}
                }
            }
        }
        std::cout << mults << std::endl;
        T buf;
        T max;
        int iter = 0;
        //while (max > eps*eps) {
        //    T sigma = B[N-1][N-1];
        //    for (int i = 0; i < N; i++) {
        //        B[i][i] -= sigma;
        //    }
        //    B = qr_decomposition(eps, B, N);
        //    for (int i = 0; i < N; i++) {
        //        B[i][i] += sigma;
        //    }
        //    max = 0;
        //    for (int k = 0; k < N-1; k++) // Нахождение максимального по модулю элемента текущей строки слева A[i][i]
        //    {
        //        buf = abs(B[N-1][k]);
        //        if (buf > max) {
        //            max = buf;
        //        }
        //    }
        //}
        //for (int i = 0; i < N; i++) {
        //    std::cout << B[i][i] << std::endl;
        //}
        for (int j = N - 1; j >= 0; j--) {
            max = 1;
            //std::cout << j << std::endl;
            while (max > 10*eps) {
                iter += 1;
                T sigma = B[j][j];
                if (flag2) {
                    for (int i = 0; i < j + 1; ++i) {
                        B[i][i] -= sigma;
                    }
                }
                //std::cout << "A " << j << std::endl;
                B = qr_decomposition(eps, B, j + 1, mults);
                //std::cout << "B " << j << std::endl;
                if (flag2) {
                    for (int i = 0; i < j + 1; ++i) {
                        B[i][i] += sigma;
                    }
                }
                //std::cout << "C " << j << std::endl;
                //max = 0;
                buf = 0;
                for (int k = 0; k < j; k++) // Нахождение максимального по модулю элемента текущей строки слева A[i][i]
                {
                    if (flag5) {
                        buf += B[j][k] * B[j][k];
                    } else {
                        buf += abs(B[j][k]);
                    }
                }
                if (flag5) {
                    max = sqrt(buf);
                }
                else {
                    max = buf;
                }
                //max = sqrt(max);
                //std::cout << max << std::endl;
                /*if (iter > 150) {
                    break;
                }*/
                //std::cout << "D " << j << std::endl;
            }
            eigenvec[j] = B[j][j];
        }
        std::sort(begin(eigenvec), end(eigenvec));
        std::vector<T> rqp = { 0.00109315, 0.226675, 21.4744, 1228.97, 47296.9, 1.2875 * pow(10, 6), 2.53089 * pow(10, 7), 
            3.57418 * pow(10, 8), 3.4293 * pow(10, 9), 1.75192 * pow(10, 10) };
        for (int i = 0; i < N; i++) {
            std::cout << "eigen numbers " << eigenvec[i] << std::endl;
            std::cout << "error" << determinant(A, eigenvec[i]) << std::endl;
            //std::cout << "true_difference " << abs(eigenvec[i] - rqp[i]) << std::endl;
        }
        std::cout << iter << std::endl;
        std::cout << mults << std::endl;
        std::cout << std::endl;

    }
    void QR_method_with_error(T eps, bool flag1, bool flag2, bool flag3) { //flag1 - Хессенберг, flag2 - cдвиг, flag3 - критерий останова 
        //flag3 - true - корень, false - сумма модулей.
        matrix<T> B = A;
        //matrix<T> E;
        std::vector<T> eigenvec(N, 0);
        T eps1 = eps;
        int mults = 0;
        /*E.reserve(N);
        for (int i = 0; i < N; i++) {
            E.emplace_back(std::vector<T>(N, 0));
        }*/
        double d;
        double alpha;
        double betha;
        if (flag1) {
            for (int i = 0; i < N - 2; i++) {
                for (int j = i + 2; j < N; j++) {
                    if (abs(B[i + 1][i]) < eps && abs(B[j][i]) < eps)
                    {
                        continue;
                    }
                    //else {
                    mults += 4;
                    d = sqrt(B[i + 1][i] * B[i + 1][i] + B[j][i] * B[j][i]);
                    alpha = B[i + 1][i] / d;
                    betha = B[j][i] / d;
                    transposemultiply(B, i, j, N, alpha, betha, mults);
                    //mults += 8 * (N - i);
                //}
                }
            }
        }
        std::cout << mults << std::endl;
        T buf;
        T max;
        int iter = 0;
        //while (max > eps*eps) {
        //    T sigma = B[N-1][N-1];
        //    for (int i = 0; i < N; i++) {
        //        B[i][i] -= sigma;
        //    }
        //    B = qr_decomposition(eps, B, N);
        //    for (int i = 0; i < N; i++) {
        //        B[i][i] += sigma;
        //    }
        //    max = 0;
        //    for (int k = 0; k < N-1; k++) // Нахождение максимального по модулю элемента текущей строки слева A[i][i]
        //    {
        //        buf = abs(B[N-1][k]);
        //        if (buf > max) {
        //            max = buf;
        //        }
        //    }
        //}
        //for (int i = 0; i < N; i++) {
        //    std::cout << B[i][i] << std::endl;
        //}
        for (int j = N - 1; j >= 0; j--) {
            max = 1;
            //std::cout << j << std::endl;
            while (max > eps) {
                iter += 1;
                T sigma = B[j][j];
                if (flag2) {
                    for (int i = 0; i < j + 1; ++i) {
                        B[i][i] -= sigma;
                    }
                }
                //std::cout << "A " << j << std::endl;
                B = qr_decomposition(eps, B, j + 1, mults);
                //std::cout << "B " << j << std::endl;
                if (flag2) {
                    for (int i = 0; i < j + 1; ++i) {
                        B[i][i] += sigma;
                    }
                }
                //std::cout << "C " << j << std::endl;
                max = 0;
                buf = 0;
                for (int k = 0; k < j; ++k) // Нахождение максимального по модулю элемента текущей строки слева A[i][i]
                {
                    if (flag3) {
                        buf += B[k+1][k] * B[k+1][k];
                    }
                    else {
                        buf += abs(B[k+1][k]);
                    }
                }
                if (!flag3) {
                    max = sqrt(buf);
                }
                else {
                    max = buf;
                }
                //max = sqrt(max);
                //std::cout << "D " << j << std::endl;
            }
            eigenvec[j] = B[j][j];
        }
        std::cout << iter << std::endl;
        std::cout << mults << std::endl;
        //T tmp;
        //T lambda0;
        //std::string str = "AAA";
        //std::vector<T> ek;
        //ek.resize(N);
        //std::vector<T> prevvec;
        //prevvec.resize(N);
        ////matrix<T> E;
        ////E.reserve(N);
        ///*for (int i = 0; i < N; i++) {
        //    E[i][i] = 1;
        //}*/
        //std::vector<T> l;
        //std::vector<T> r;
        //for (int j = 0; j < N; j++) {
        //    /*if (j <= 1) {
        //        eps = eps * eps;
        //    }*/
        //    lambda0 = eigenvec[j];
        //    //std::cout << lambda0 << std::endl;
        //    B = A;
        //    for (int i = 0; i < N; i++) {
        //        prevvec[i] = 0;
        //    }
        //    //printmatrix(B, N, eps);
        //    //std::cout << std::endl;
        //    for (int i = 0; i < N; i++) {
        //        B[i][i] = B[i][i] - lambda0;
        //    }
        //    //printmatrix(B, N, eps);
        //    B = Inverse(B, eps * eps, str);
        //    //printmatrix(B, N, eps);
        //    //std::cout << std::endl;
        //    for (int i = 0; i < N; i++) {
        //        ek[i] = 0;
        //    }
        //    ek[0] = 1;
        //    l = ek - prevvec;
        //    r = ek + prevvec;
        //    while (true) {
        //        prevvec = ek;
        //        ek = B * ek;
        //        tmp = 1 / vector_norm(ek);
        //        ek = ek * tmp;
        //        l = ek - prevvec;
        //        r = ek + prevvec;
        //        /*printvector(ek);
        //        printvector(l);
        //        std::cout << vector_norm(l) << std::endl;
        //        printvector(r);
        //        std::cout << vector_norm(r) << std::endl;*/
        //        if (vector_norm(l) < eps || vector_norm(r) < eps) { //(vector_norm(r) < eps)
        //            break;
        //        }
        //    }
        //    ek = ek * (-1.0);
        //    std::vector<T> yvec = A * ek;
        //    std::vector<T> tmpvecc = ek*eigenvec[j];
        //    yvec = yvec - tmpvecc;
        //    std::cout << vector_norm(yvec) << std::endl;
        //    /*std::cout << "eigen value number " << j + 1 << " " << eigenvec[j] << std::endl;
        //    std::cout << "eigen vector number " << j + 1 << std::endl;
        //    printvector(ek);*/
        //    
        //    //this->input(str);
        //}
        //for (int i = 0; i < N; i++) {
        //    std::cout << "eigen numbers " << eigenvec[i] << std::endl;
        //}
        //std::cout << iter << std::endl;
        //std::cout << mults << std::endl;
       
    }
    void method_backward_iterations(std::vector<T> lambda, T eps, std::string str) {
        //printmatrix(A, N, eps);
        /*T tmp = 1 / vector_norm(e0);
        e0 = e0*tmp;*/
        matrix<T> B;
        T tmp;
        T lambda0;
        std::vector<T> ek;
        ek.resize(N);
        std::vector<T> prevvec;
        prevvec.resize(N);
        matrix<T> E;
        E.reserve(N);
        for (int i = 0; i < N; i++) {
            E.emplace_back(std::vector<T>(N, 0));
        }
        for (int i = 0; i < N; i++) {
            E[i][i] = 1;
        }
        std::vector<T> l;
        std::vector<T> r;
        for (int j = 0; j < N; j++) {
            lambda0 = lambda[j];
            std::cout << lambda0 << std::endl;
            B = A;
            for (int i = 0; i < N; i++) {
                prevvec[i] = 0;
            }
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            for (int i = 0; i < N; i++) {
                B[i][i] = B[i][i] - lambda0;
            }
            //printmatrix(B, N, eps);
            B = Inverse(B, eps, str);
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            for (int i = 0; i < N; i++) {
                ek[i] = 0;
            }
            ek[0] = 1;
            l = ek - prevvec;
            r = ek + prevvec;
            while (true) {
                prevvec = ek;
                ek = B * ek;
                tmp = 1 / vector_norm(ek);
                ek = ek * tmp;
                l = ek - prevvec;
                r = ek + prevvec;
                /*printvector(ek);
                printvector(l);
                std::cout << vector_norm(l) << std::endl;
                printvector(r);
                std::cout << vector_norm(r) << std::endl;*/
                if (vector_norm(l) < eps || vector_norm(r) < eps) { //(vector_norm(r) < eps)
                    break;
                }
            }
            ek = ek * (-1.0);
            std::cout << "eigen vector number " << j + 1 << std::endl;
            printvector(ek);
            //this->input(str);
        }
    }
    void special_method_backward_iterations(std::vector<T> lambda, matrix<T> e0, T eps, std::string str) {
        //printmatrix(A, N, eps);
        /*T tmp = 1 / vector_norm(e0);
        e0 = e0*tmp;*/
        matrix<T> B;
        T tmp;
        T lambda0;
        std::vector<T> ek;
        ek.resize(N);
        std::vector<T> prevvec;
        prevvec.resize(N);
        matrix<T> E;
        E.reserve(N);
        for (int i = 0; i < N; i++) {
            E.emplace_back(std::vector<T>(N, 0));
        }
        for (int i = 0; i < N; i++) {
            E[i][i] = 1;
        }
        std::vector<T> l;
        std::vector<T> r;
        for (int j = 0; j < N; j++) {
            lambda0 = lambda[j];
            std::cout << lambda0 << std::endl;
            B = A;
            for (int i = 0; i < N; i++) {
                prevvec[i] = 0;
            }
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            for (int i = 0; i < N; i++) {
                B[i][i] = B[i][i] - lambda0;
            }
            //printmatrix(B, N, eps);
            B = Inverse(B, eps, str);
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            /*if (j == 0) {*/
            ek = e0[j];
              /*  for (int i = 0; i < N; i++) {
                    ek[i] = 0;
                }
                ek[0] = 1;*/
            //}
            l = ek - prevvec;
            r = ek + prevvec;
            while (true) {
                prevvec = ek;
                ek = B * ek;
                tmp = 1 / vector_norm(ek);
                ek = ek * tmp;
                l = ek - prevvec;
                r = ek + prevvec;
                printvector(ek);
                /*printvector(ek);
                printvector(l);
                std::cout << vector_norm(l) << std::endl;
                printvector(r);
                std::cout << vector_norm(r) << std::endl;*/
                if (vector_norm(l) < eps || vector_norm(r) < eps) { //(vector_norm(r) < eps)
                    break;
                }
            }
            ek = ek * (-1.0);
            std::cout << "eigen vector number " << j + 1 << std::endl;
            printvector(ek);
            //this->input(str);
        }
    }
    void method_backward_iterations_with_Relay(matrix<T> e0, T eps, std::string str) {
        //printmatrix(A, N, eps);
        /*T tmp = 1 / vector_norm(e0);
        e0 = e0*tmp;*/
        matrix<T> B;
        T tmp;
        T lambda0;
        std::vector<T> ek;
        ek.resize(N);
        std::vector<T> prevvec;
        prevvec.resize(N);
        matrix<T> E;
        E.reserve(N);
        for (int i = 0; i < N; i++) {
            E.emplace_back(std::vector<T>(N, 0));
        }
        std::vector<T> l;
        std::vector<T> r;
        std::vector<T> tmpvec;
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                prevvec[i] = 0;
            }
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            //printmatrix(B, N, eps);
            //B = Inverse(B, eps, str);
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            ek = e0[j];
            tmp = 1 / vector_norm(ek);
            ek = ek * tmp;
            l = ek - prevvec;
            r = ek + prevvec;
            while (true) {
                B = A;
                tmpvec = B * ek;
                lambda0 = tmpvec * ek;
                for (int i = 0; i < N; i++) {
                    B[i][i] -= lambda0;
                }
                //l = B * ek;
                //r = ek * lambda0;
                //printvector(l);
                //printvector(r);
                prevvec = ek;
                ek = Gauss(B, ek, eps, false);
                tmp = 1 / vector_norm(ek);
                ek = ek * tmp;
                l = ek - prevvec;
                r = ek + prevvec;
                if (vector_norm(l) < eps || vector_norm(r) < eps) {
                    break;
                }
                //printvector(l);
                //std::cout << vector_norm(l) << std::endl;
                /*printvector(ek);
                printvector(l);
                std::cout << vector_norm(l) << std::endl;
                printvector(r);
                std::cout << vector_norm(r) << std::endl;*/
            }
            ek = ek * (-1.0);
            std::cout << "eigen value number " << j + 1 << " " << lambda0 << std::endl;
            std::cout << "eigen vector number " << j + 1 << std::endl;
            printvector(ek);
            //this->input(str);
        }
    }
    void special_method_backward_iterations_with_Relay(matrix<T> e0, T eps, std::string str) {
        //printmatrix(A, N, eps);
        /*T tmp = 1 / vector_norm(e0);
        e0 = e0*tmp;*/
        matrix<T> B;
        std::vector<T> tmvec;
        T tmp;
        T lambda0;
        std::vector<T> ek;
        ek.resize(N);
        std::vector<T> prevvec;
        prevvec.resize(N);
        matrix<T> E;
        E.reserve(N);
        for (int i = 0; i < N; i++) {
            E.emplace_back(std::vector<T>(N, 0));
        }
        std::vector<T> l;
        std::vector<T> r;
        std::vector<T> tmpvec;
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                prevvec[i] = 0;
            }
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            //printmatrix(B, N, eps);
            //B = Inverse(B, eps, str);
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            if (j == 0) {
                ek = e0[j];
            }
            tmp = 1 / vector_norm(ek);
            ek = ek * tmp;
            l = ek - prevvec;
            r = ek + prevvec;
            while (true) {
                B = A;
                tmpvec = B * ek;
                lambda0 = tmpvec * ek;
                for (int i = 0; i < N; i++) {
                    B[i][i] -= lambda0;
                }
                //l = B * ek;
                //r = ek * lambda0;
                //printvector(l);
                //printvector(r);
                prevvec = ek;
                ek = Gauss(B, ek, eps, false);
                tmp = 1 / vector_norm(ek);
                ek = ek * tmp;
                l = ek - prevvec;
                r = ek + prevvec;
                if (vector_norm(l) < eps || vector_norm(r) < eps) {
                    break;
                }
                //printvector(l);
                //std::cout << vector_norm(l) << std::endl;
                /*printvector(ek);
                printvector(l);
                std::cout << vector_norm(l) << std::endl;
                printvector(r);
                std::cout << vector_norm(r) << std::endl;*/
            }
            ek = ek * (-1.0);
            tmvec = ek;
            std::cout << "eigen value number " << j + 1 << " " << lambda0 << std::endl;
            std::cout << "eigen vector number " << j + 1 << std::endl;
            printvector(ek);
            //this->input(str);
        }
    }
    void Modification_method_eig_val_vec(T eps, std::string str) {
        matrix<T> B = A;
        matrix<T> E;
        std::vector<T> eigenvec(N, 0);
        T eps1 = eps;
        int mults = 0;
        E.reserve(N);
        for (int i = 0; i < N; i++) {
            E.emplace_back(std::vector<T>(N, 0));
        }
        double d;
        double alpha;
        double betha;
        for (int i = 0; i < N - 2; i++) {
            for (int j = i + 2; j < N; j++) {
                //printmatrix(B, N, eps);
                //std::cout << std::endl;
                mults += 4;
                d = sqrt(B[i + 1][i] * B[i + 1][i] + B[j][i] * B[j][i]);
                alpha = B[i + 1][i] / d;
                betha = B[j][i] / d;
                transposemultiply(B, i, j, N, alpha, betha);
                mults += 4 * N;
            }
        }
        T buf;
        T max = 1;
        //while (max > eps*eps) {
        //    T sigma = B[N-1][N-1];
        //    for (int i = 0; i < N; i++) {
        //        B[i][i] -= sigma;
        //    }
        //    B = qr_decomposition(eps, B, N);
        //    for (int i = 0; i < N; i++) {
        //        B[i][i] += sigma;
        //    }
        //    max = 0;
        //    for (int k = 0; k < N-1; k++) // Нахождение максимального по модулю элемента текущей строки слева A[i][i]
        //    {
        //        buf = abs(B[N-1][k]);
        //        if (buf > max) {
        //            max = buf;
        //        }
        //    }
        //}
        //for (int i = 0; i < N; i++) {
        //    std::cout << B[i][i] << std::endl;
        //}
        for (int j = N - 1; j >= 0; j--) {
            max = 1;
            //std::cout << j << std::endl;
            while (max > eps) {
                T sigma = B[j][j];
                for (int i = 0; i < j + 1; ++i) {
                    B[i][i] -= sigma;
                }
                //std::cout << "A " << j << std::endl;
                B = qr_decomposition(eps, B, j + 1);
                //std::cout << "B " << j << std::endl;
                for (int i = 0; i < j + 1; ++i) {
                    B[i][i] += sigma;
                }
                //std::cout << "C " << j << std::endl;
                max = 0;
                for (int k = 0; k < j; ++k) // Нахождение максимального по модулю элемента текущей строки слева A[i][i]
                {
                    buf = abs(B[j][k]);
                    if (buf > max) {
                        max = buf;
                    }
                }
                //std::cout << "D " << j << std::endl;
            }
            eigenvec[j] = B[j][j];
        }
        printvector(eigenvec);
        //matrix<T> B;
        T tmp;
        T lambda0;
        std::vector<T> ek;
        ek.resize(N);
        std::vector<T> prevvec;
        prevvec.resize(N);
        //matrix<T> E;
        //E.reserve(N);
        for (int i = 0; i < N; i++) {
            E[i][i] = 1;
        }
        std::vector<T> l;
        std::vector<T> r;
        for (int j = 0; j < N; j++) {
            /*if (j <= 1) {
                eps = eps * eps;
            }*/
            lambda0 = eigenvec[j];
            //std::cout << lambda0 << std::endl;
            B = A;
            for (int i = 0; i < N; i++) {
                prevvec[i] = 0;
            }
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            for (int i = 0; i < N; i++) {
                B[i][i] = B[i][i] - lambda0;
            }
            //printmatrix(B, N, eps);
            B = Inverse(B, eps * eps, str);
            //printmatrix(B, N, eps);
            //std::cout << std::endl;
            for (int i = 0; i < N; i++) {
                ek[i] = 0;
            }
            ek[0] = 1;
            l = ek - prevvec;
            r = ek + prevvec;
            while (true) {
                prevvec = ek;
                ek = B * ek;
                tmp = 1 / vector_norm(ek);
                ek = ek * tmp;
                l = ek - prevvec;
                r = ek + prevvec;
                /*printvector(ek);
                printvector(l);
                std::cout << vector_norm(l) << std::endl;
                printvector(r);
                std::cout << vector_norm(r) << std::endl;*/
                if (vector_norm(l) < eps || vector_norm(r) < eps) { //(vector_norm(r) < eps)
                    break;
                }
            }
            ek = ek * (-1.0);
            std::cout << "eigen value number " << j + 1 << " " << eigenvec[j] << std::endl;
            std::cout << "eigen vector number " << j + 1 << std::endl;
            printvector(ek);
            //this->input(str);
        }
    }
};
//std::vector<std::vector <T>> qr_decomposition(std::vector<std::vector <T>>& a, std::vector<T> b, int N) {
//    /* std::vector<std::vector<double>> r = { {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0} };
//     std::vector<std::vector<double>> q = { {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0} };
//     std::vector<std::vector<double>> a_ort = { {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0} };*/
//     /*std::vector<std::vector<double>> r = { {0, 0}, {0, 0} };
//     std::vector<std::vector<double>> q = { {0, 0}, {0, 0} };
//     std::vector<std::vector<double>> a_ort = { {0, 0}, {0, 0} };*/
//    T max_element, tmp;
//    int max_number;
//    T epsilon = pow(10, -10);
//    std::vector<std::vector<T>> r, q, a_ort;
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            r[i][j] = 0;
//            q[i][j] = 0;
//            a_ort[i][j] = 0;
//        }
//    }
//    T projections_used = 0;
//    std::vector<T> projections;
//    //std::cout << "first q column: ";
//    for (int i = 0; i < N; i++) {
//        q[i][0] = a[i][0] / column_norm(a, 0);
//        //std::cout << q[i][0] << " ";
//    }
//    //std::cout << std::endl;
//    T projection_count = 0;
//    for (int i = 1; i < N; i++) {
//        for (int j = 0; j < i; j++) {
//            projections.push_back(vector_projection(a, i, q, j));
//            //std::cout << "new projection: " << vector_projection(a, i, q, j)<<std::endl;
//        }
//        for (int j = 0; j < N; j++)
//            a_ort[j][i] = a[j][i];
//        for (int k = 0; k < projections.size(); k++) //(*)
//            for (int j = 0; j < N; j++) {
//                a_ort[j][i] -= projections[k] * q[j][k];
//                //std::cout << "a_ort["<<j<<"]["<<i<<"] = " << a_ort[j][i] << " ";
//            }
//        //std::cout << std::endl;
//        projections.clear();
//        projections_used = projections.size(); //изменяем индекс первой нужной проекции вектора для начала (*)
//        //std::cout << "second q column: ";
//        for (int j = 0; j < N; j++) {
//            q[j][i] = a_ort[j][i] / column_norm(a_ort, i);
//            //std::cout << q[j][i] << " ";
//        }
//        //std::cout << std::endl;
//
//    }
//
//    T rowcount = 0;
//    r[0][0] = column_norm(a, 0);
//    for (int i = 1; i < N; i++)
//        r[i][i] = column_norm(a_ort, i);
//    rowcount++;
//    for (int i = 1; i < N; i++) {
//        for (int j = 0; j < rowcount; j++) {
//            r[j][i] = vector_projection(a, i, q, j);
//        }
//        rowcount++;
//    }
//
//    //std::vector<std::vector<T>> test_A;
//    //for (int i = 0; i < N; i++) {
//    //    for (int j = 0; j < N; j++) {
//    //        test_A = 
//    //    }
//    //}
//    ////std::vector<std::vector<double>> test_A = { {0, 0}, {0, 0} };
//
//    //for (int i = 0; i < a.size(); i++) {
//    //    for (int j = 0; j < a.size(); j++) {
//    //        for (int inner = 0; inner < a.size(); inner++)
//    //            test_A[i][j] += q[i][inner] * r[inner][j];
//    //    }
//    //}
//    std::vector<T> btmp = transpose(q) * b;
//    for (int k = 0; k < N; k++) {
//        max_element = abs(r[k][k]);
//        max_number = k; 
//        for (int i = k + 1; i < N; i++) {
//            if (abs(r[i][k]) > max_element) {
//                max_element = abs(r[i][k]);
//                max_number = i;
//            }
//        }
//        if (max_element < epsilon) {
//            std::cout << "BAD MATRIX" << std::endl;
//            exit(1);
//        }
//        if (k != max_number) {
//            for (int i = 0; i < N; i++) {
//                tmp = r[k][i];
//                r[k][i] = r[max_number][i];
//                r[max_number][i] = tmp;
//            }
//            tmp = btmp[k];
//            btmp[k] = btmp[max_number];
//            btmp[max_number] = tmp;
//        }
//    }
//    std::vector<T> x;
//    x.resize(N);
//    for (int i = 0; i < N - 1; i++) {
//        x[i] = 0;
//    }
//    x[N - 1] = btmp[N - 1] / r[N - 1][N - 1];
//    for (int k = N - 2; k >= 0; k--) {// обратный ход
//        x[k] = btmp[k];
//        for (int l = N - 1; l > k; l--) {
//            x[k] = x[k] - r[k][l] * x[l];
//        }
//        x[k] = x[k] / r[k][k];
//        /*for (int i = 0; i < N; i++) {
//            for (int j = 0; j < N; j++) {
//                std::cout << A[i][j] << " ";
//            }
//            std::cout << B[i] << std::endl;
//        }*/
//
//    }
//    for (int i = 0; i < N; i++) {
//        std::cout << x[i] << std::endl;
//    }
//}
//int main()
//{
//    std::vector<std::vector<double>> A = { {1, 2, 3, 4, 5, 6, 7}, {8, 9, 10, 11, 12, 13, 14}, {15, 16, 17, 18, 19, 20, 21}, {22, 23, 24, 25, 26, 27, 28}, {29, 30, 31, 32, 33, 34, 35}, {36, 37, 38, 39, 40, 41, 42}, {43, 44, 45, 46, 47, 48, 49} };
//
//    //std::vector<std::vector<double>> A = { {1, 3}, {1, -1} };
//
//    std::vector<std::vector<double>> answer = qr_decomposition(A);
//    for (int i = 0; i < A.size(); i++) {
//        for (int j = 0; j < A.size(); j++) {
//            std::cout << answer[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
//}
int main()
{
    std::string str_input = "D:\\input1_laba3.txt";
    std::string str_output = "D:\\output_laba1.txt";
    Eigenvecandval<double> C; // -0.0119071, 0.709759, -0.607555, 0.356337
    matrix<double> e1 = { {0.5, 0.5, 0.5, 0.5}, {-0.5, -0.01, 0.4, 0.7}, {0, -0.5, -0.5, 0.2}, {-0.8, 0, -0.2, -0.4} }; //для файла инпут 1;
    matrix<double> e2 = { {0.5, 0.5, 0.5, 0.5}, {-0.5, -0.01, 0.4, 0.7}, {0, -0.9, -0.01, 0}, {0, 0, -0.9, -0.1} }; //для файла инпут 2;
    matrix<double> e3 = { {0.5, 0.5, 0.5, 0.5}, {-0.3, -1, 0, 0}, {-0.3, 0, -1, 0}, {-0.8, 0, -0.2, -0.4} }; //для файла инпут 3;
    double epsilon = pow(10, -10);
    double epsilon1 = pow(10, -3);
    double epsilon2 = pow(10, -9);
    double epsilon3 = pow(10, -7);
    C.input(str_input);
    //C.spec_input(10);
    std::vector<double> lambda1 = {0.6, 1.6, 2.6, 3.6}; //для файла инпут 1
    std::vector<double> lambda2 = {170.0, -184.0, -160.0, -90.0}; //для файла инпут 2
    std::vector<double> lambda3 = {180, -120, 100, -90}; //для файла инпут 3
    matrix<double> e10 = { { 0.00343211, -0.704259, -0.620789, 0.344426 }, { -0.502548, -0.0157964, 0.430896, 0.749349 }, { -0.0119071, 0.709759, -0.607555, 0.356337 }, { -0.864461, -0.00338916, -0.244595, -0.439169 } };
    //flag1 - Хессенберг, flag2 - cдвиг, flag3 - критерий останова 
    //flag3 - true - корень, false - сумма модулей.
    /*std::cout << "Hessenberg " << "shift " << "square" << std::endl;
    C.QR_method(epsilon3, true, true, true);
    std::cout << "Without Hessenberg " << "shift " << "square" << std::endl;
    C.QR_method(epsilon3, false, true, true);
    std::cout << "Hessenberg " << "without shift " << "square" << std::endl;
    C.QR_method(epsilon3, true, false, true);
    std::cout << "Without Hessenberg " << "without shift " << "square" << std::endl;
    C.QR_method(epsilon3, false, false, true);
    std::cout << "Hessenberg " << "shift " << "module" << std::endl;
    C.QR_method(epsilon2, true, true, false);
    std::cout << "Without Hessenberg " << "shift " << "moudle" << std::endl;
    C.QR_method(epsilon2, false, true, false);
    std::cout << "Hessenberg " << "without shift " << "module" << std::endl;
    C.QR_method(epsilon2, true, false, false);
    std::cout << "Without Hessenberg " << "without shift " << "module" << std::endl;
    C.QR_method(epsilon2, false, false, false);
    std::cout << std::endl;*/
    C.special_method_backward_iterations(lambda1, e10, epsilon, str_input);
    //C.QR_method_with_error(epsilon2, true, true, false);
    /*s td::cout << "Backward iterations" << std::endl;
    C.method_backward_iterations(lambda2, epsilon, str_input);*/
    std::cout << "Backward iterations with Relay" << std::endl;
    C.method_backward_iterations_with_Relay(e1, epsilon1, str_input);
    /*std::cout << "Modification_method_eig_val_vec" << std::endl;
    C.Modification_method_eig_val_vec(epsilon, str_input);*/
    //модификация обратных итераций.
}
