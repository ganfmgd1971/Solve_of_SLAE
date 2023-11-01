// ConsoleApplication87.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
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
            if (A[i][j] < epsilon) {
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
template<typename T>
T vector_norm(std::vector<T>& a) {
    int n = a.size();
    T result = 0;
    for (int j = 1; j < n; j++) {
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
template <typename T>
class Hsole { //Heterogeneous system of linear equations
    matrix<T> A;
    int N;
    std::vector<T> B;
    std::vector<T> x;
public:
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
                out << B[i] << " ";
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
        //printvector(y);
            /*for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    std::cout << A[i][j] << " ";
                }
                std::cout << B[i] << std::endl;
            }*/
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
                    /*tmp = A[k][i];
                    A[k][i] = A[max_number][i];
                    A[max_number][i] = tmp;*/
                }
                swap(B[k], B[max_number]);
                /*tmp = B[k];
                    B[k] = B[max_number];
                    B[max_number] = tmp;*/
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
        T tmp = residual(h, B);
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
        std::cout << "multiply of matrix and inverse matrix" << std::endl;
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
    std::string str_input = "D:\\inputlaba1.txt";
    std::string str_output = "D:\\output_laba1.txt";
    //std::cout << "double, N = 5" << std::endl;
    //Hsole<double> C;
    //double epsilon = pow(10, -16);
    //C.special_input(5);
    //C.Gauss(epsilon, true);
    //C.check_solve();
    //C.special_input(5);
    //C.QR_solve(epsilon, true);
    //C.check_solve();
    //matrix<double> Ainv = C.Inverse(epsilon, str_input);
    //std::cout << "Inverse matrix" << std::endl;
    //printmatrix(Ainv, Ainv.size(),epsilon);
    //C.check_multiply(epsilon, str_input);
    //std::cout << std::endl;
    //std::cout << "double, N = 10" << std::endl;
    //Hsole<double> D;
    //epsilon = pow(10, -25);
    ///*std::string str_input = "D:\\inputlaba1.txt";
    //std::string str_output = "D:\\output_laba1.txt";*/
    //D.special_input(10);
    //D.Gauss(epsilon, true);
    //D.check_solve();
    //D.special_input(10);
    //D.QR_solve(epsilon, true);
    //D.check_solve();
    //matrix<double> Ainv1 =D.Inverse(epsilon, str_input);
    //std::cout << "Inverse matrix" << std::endl;
    //printmatrix(Ainv1, Ainv1.size(), epsilon);
    //D.check_multiply(epsilon, str_input);
    /*epsilon = pow(10, -25);*/
    //D.check_multiply(epsilon, str_input);
    std::cout << std::endl;
    std::cout << "float, N = 5" << std::endl;
    Hsole<float> E;
    float epsilon1 = pow(10, -10);
    //*std::string str_input = "D:\\inputlaba1.txt";
    //std::string str_output = "D:\\output_laba1.txt";*/
    E.special_input(5);
    E.Gauss(epsilon1, true);
    E.check_solve();
    E.special_input(5);
    E.QR_solve(epsilon1, true);
    E.check_solve();
    matrix<float> Ainv2 = E.Inverse(epsilon1, str_input);
    std::cout << "Inverse matrix" << std::endl;
    printmatrix(Ainv2, Ainv2.size(), epsilon1);
    E.check_multiply(epsilon1, str_input);
    //D.check_multiply(epsilon1, str_input);
    std::cout << std::endl;
    std::cout << "float, N = 10" << std::endl;
    Hsole<float> F;
    epsilon1 = pow(10, -15);
    //*std::string str_input = "D:\\inputlaba1.txt";
    //std::string str_output = "D:\\output_laba1.txt";*/
    F.special_input(10);
    F.Gauss(epsilon1, true);
    F.check_solve();
    F.special_input(10);
    F.QR_solve(epsilon1, true);
    F.check_solve();
    matrix<float> Ainv3 = F.Inverse(epsilon1, str_input);
    std::cout << "Inverse matrix" << std::endl;
    printmatrix(Ainv3, Ainv3.size(), epsilon1);
    F.check_multiply(epsilon1, str_input);
    /* C.input(str_input);
     C.Gauss(epsilon, true);
     C.output(str_output);
     C.input(str_input);
     C.QR_solve(epsilon, true);
     C.output(str_output);
     C.input(str_input);
     matrix<double> Ainv = C.Inverse(epsilon, str_input);
     std::cout << "Inverse matrix" << std::endl;
     printmatrix(Ainv, Ainv.size());
     C.input(str_input);
     std::cout << std::endl;
     C.checkresidual(epsilon);
     C.output(str_output);
     C.input(str_input);
     std::cout << std::endl;
     C.check_conditionality(epsilon, str_input);
     C.input(str_input);
     C.evaluate_conditionality(false, epsilon, str_input);
     C.check_conditionality(epsilon, str_input);
     C.check_multiply(epsilon, str_input);*/
     /* Hsole<float> C;
      double epsilon = pow(10, -16);
      std::string str_input = "D:\\input2_laba1.txt";
      std::string str_output = "D:\\output_laba1.txt";
      C.input(str_input);
      C.Gauss(epsilon, true);
      C.output(str_output);
      C.input(str_input);
      C.QR_solve(epsilon, true);
      C.output(str_output);
      C.input(str_input);
      matrix<double> Ainv = C.Inverse(epsilon, str_input);
      std::cout << "Inverse matrix" << std::endl;
      printmatrix(Ainv, Ainv.size());
      C.input(str_input);
      std::cout << std::endl;
      C.checkresidual(epsilon);
      C.output(str_output);
      C.input(str_input);
      std::cout << std::endl;
      C.check_conditionality(epsilon, str_input);
      C.input(str_input);
      C.evaluate_conditionality(false, epsilon, str_input);
      C.check_conditionality(epsilon, str_input);
      C.check_multiply(epsilon, str_input);*/
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.