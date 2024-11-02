#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

int PRECISION;
const int ITERATION_LIMIT = 10000;


class Matrix {
protected:
    int n, m;
public:
    vector<vector<double>> matrixData;
    Matrix(int n, int m) : n(n), m(m), matrixData(vector<vector<double>>(n, vector<double>(m, 0))) {}

    friend istream& operator>>(istream& in, Matrix& matrixToRead) {
        for (int i = 0; i < matrixToRead.n; ++i) {
            for (int j = 0; j < matrixToRead.m; ++j) {
                in >> matrixToRead.matrixData[i][j];
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const Matrix& matrix) {
        for (int i = 0; i < matrix.n; ++i) {
            for (int j = 0; j < matrix.m; ++j) {
                out << fixed << setprecision(PRECISION) << matrix.matrixData[i][j];
                if (j != matrix.m - 1) {
                    out << " ";
                }
            }
            out << endl;
        }
        return out;
    }
    Matrix operator=(Matrix& other) {
        n = other.n;
        m = other.m;
        matrixData = other.matrixData;
        return *this;
    }

    virtual Matrix* operator+(Matrix& other) {
        Matrix* newMatrix = new Matrix(n, m);
        if (n == other.n && m == other.m) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    newMatrix->matrixData[i][j] = matrixData[i][j] + other.matrixData[i][j];
                }
            }
            return newMatrix;
        }
        else {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
    }

    virtual Matrix* operator-(Matrix& other) {
        Matrix* newMatrix = new Matrix(n, m);
        if (n == other.n && m == other.m) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    newMatrix->matrixData[i][j] = matrixData[i][j] - other.matrixData[i][j];
                }
            }
            return newMatrix;
        }
        else {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
    }

    virtual Matrix* operator*(Matrix& other) {
        Matrix* newMatrix = new Matrix(n, other.m);
        if (m == other.n) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < other.m; j++) {
                    for (int p = 0; p < m; p++) {
                        newMatrix->matrixData[i][j] += matrixData[i][p] * other.matrixData[p][j];
                    }
                }
            }
            return newMatrix;
        }
        else {
            cout << "Error: the dimensional problem occurred" << endl;
            return nullptr;
        }
    }

    virtual Matrix* transpose() {
        Matrix* newMatrix = new Matrix(m, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                newMatrix->matrixData[j][i] = matrixData[i][j];
            }
        }
        return newMatrix;
    }
    void setDiagonal(Matrix& coef) {
        for (int i = 0; i < n; i++){
            this->matrixData[i][i] = coef.matrixData[i][0];
        }
    }
    Matrix* inverse();

    int getRows() const { return n; }

    int getColumns() const { return m; }
};

class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int n) : Matrix(n, n) {
        for (int i = 0; i < n; i++) {
            matrixData[i][i] = 1;
        }
    }
};

Matrix* concatenate(Matrix& A, Matrix& B, bool dim) {
    Matrix* newMatrix;
    if (dim){
        newMatrix = new Matrix(A.getRows(), A.getColumns() + B.getColumns());
        for (int i = 0; i < A.getRows(); i++) {
            for (int j = 0; j < A.getColumns(); j++) {
                newMatrix->matrixData[i][j] = A.matrixData[i][j];
            }
            for (int j = 0; j < B.getColumns(); j++) {
                newMatrix->matrixData[i][A.getColumns()+j] = B.matrixData[i][j];
            }
        }
    }
    else{
        newMatrix = new Matrix(A.getRows()+B.getRows(), A.getColumns());
        for (int i = 0; i < A.getRows(); i++) {
            for (int j = 0; j < A.getColumns(); j++) {
                newMatrix->matrixData[i][j] = A.matrixData[i][j];
            }
        }
        for (int i = 0; i < B.getRows(); i++) {
            for (int j = 0; j < B.getColumns(); j++){
                newMatrix->matrixData[i+A.getRows()][j] = B.matrixData[i][j];
            }
        }
    }

    return newMatrix;
}

Matrix* Matrix::inverse() {
    if (n != m) {
        cout << "Error: only square matrices have an inverse." << endl;
        return nullptr;
    }
    IdentityMatrix I(n);
    Matrix* augmented = concatenate(*this, I, true);

    for (int i = 0; i < n; i++) {
        // Проверяем, что диагональный элемент не равен нулю
        if (augmented->matrixData[i][i] == 0) {
            bool found = false;
            for (int j = i + 1; j < n; j++) {
                if (augmented->matrixData[j][i] != 0) {
                    swap(augmented->matrixData[i], augmented->matrixData[j]);
                    found = true;
                    break;
                }
            }
            if (!found) {
                cout << "Error: matrix is singular and cannot be inverted." << endl;
                delete augmented;
                return nullptr;
            }
        }

        // Нормализуем диагональный элемент
        double diagonalElement = augmented->matrixData[i][i];
        for (int j = 0; j < 2 * n; j++) {
            augmented->matrixData[i][j] /= diagonalElement;
        }

        // Обнуляем все элементы в текущем столбце, кроме диагонального
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double factor = augmented->matrixData[j][i];
                for (int k = 0; k < 2 * n; k++) {
                    augmented->matrixData[j][k] -= factor * augmented->matrixData[i][k];
                }
            }
        }
    }

    // Извлекаем правую часть, которая теперь является обратной матрицей
    Matrix* inverseMatrix = new Matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverseMatrix->matrixData[i][j] = augmented->matrixData[i][j + n];
        }
    }

    delete augmented;
    return inverseMatrix;
}

void makeBasicColumn(Matrix& table, int basic_row, int basic_col) {
    // Make basic element 1
    double basic_element = table.matrixData[basic_row][basic_col];
    for (int i = 0; i < table.getColumns(); i++) {
        table.matrixData[basic_row][i] /= basic_element;
    }


    // Nullify all elements except basic row
    for (int i = 0; i < table.getRows(); i++) {
        if (i == basic_row) {
            continue;
        }

        double coefficient = table.matrixData[i][basic_col];
        for (int j = 0; j < table.getColumns(); j++) {
            table.matrixData[i][j] -= coefficient * table.matrixData[basic_row][j];
        }

    }
}


class Answer {
public:
    bool solver_sate;
    Matrix solution;
    double z;

public:
    Answer(bool state, const Matrix& sol, float z): solver_sate(state), solution(sol), z(z) {}
};




Matrix getSolution(Matrix& A, Matrix& solution) {
    for (int j = 0; j < A.getColumns() - 1; j++) {
        int countOne = 0, oneInd = -1;
        bool flag = false;
        for (int i = 0; i < A.getRows(); i++) {
            if (A.matrixData[i][j] == 1 && countOne == 0) {
                countOne++;
                oneInd = i;
            }
            else if (A.matrixData[i][j] == 0) {
                continue;
            }
            else {
                flag = true;
            }
        }
        if (countOne == 1 && !flag) {
            solution.matrixData[0][j] = A.matrixData[oneInd][A.getColumns() - 1];
        }
    }
    return solution;
}



Answer solveLLP(Matrix& C, Matrix& A, Matrix& b, int eps) {
    PRECISION = eps;

    Matrix temp_col = *concatenate(*new Matrix(1, 1), b, false);

    Matrix tempMatrix = *new Matrix(A.getRows() + 1, A.getRows());
    for (int i = 0; i < A.getRows(); i++) {
        tempMatrix.matrixData[i+1][i] = 1;
    }
    Matrix table = *concatenate(*concatenate(C, A, false), tempMatrix, true);
    // cout << table << endl;

    for (int i = 0; i < A.getColumns(); i++) {
        int n = -1;
        int c = 0;
        if (table.matrixData[0][i] == 0) {
            for (int j = 1; j < table.getRows(); j++) {
                if (table.matrixData[j][i] == 1) {
                    c++;
                    n = j;
                }
            }
            // cout << c << ' ' << n << endl;
            if (c == 1) {
                for (int k = A.getColumns(); k < table.getColumns(); k++) {
                    if (table.matrixData[n][k] == 1) {
                        table.matrixData[n][k] = 0;
                        break;
                    }
                }
            }
            // cout << table << endl;
        }
    }

    table = *concatenate(table, temp_col, true);

    int n_var = table.getColumns() - 1;
    int n_constrains = table.getRows() -1;

    for (int i = 0; i < n_var; i++) {
        table.matrixData[0][i] *= -1;
    }
//    cout << table << endl;

    for (int _ = 0; _ < n_constrains; _++) {
        int basic_col = -1;
        int min_col = 0;
        for (int i = 0; i < n_var; i++) {
            if (table.matrixData[0][i] < min_col) {
                min_col = table.matrixData[0][i];
                basic_col = i;
            }
        }
        if (basic_col == -1) {
            break;
        };

        int basic_row = 0;
        int min_ratio = pow(10, 10);
        for (int i = 1; i < n_constrains; i++) {
            if (table.matrixData[i][n_var] / table.matrixData[i][basic_col] < min_ratio &&
                table.matrixData[i][n_var] / table.matrixData[i][basic_col] > 0) {
                basic_row = i;
                min_ratio = table.matrixData[i][n_var] / table.matrixData[i][basic_col];
            }
        }
        if (basic_row == 0) {
            return Answer(false, Matrix(1, n_var), -123123);
        }
        makeBasicColumn(table, basic_row, basic_col);

    }

    Matrix solution = getSolution(table, *new Matrix(1, C.getColumns()));
    for (int i = 0; i < solution.getColumns(); i++) {
        if (solution.matrixData[0][i] < 0) {
            return Answer(false, Matrix(1, n_var), -123123);
        }
    }

    return Answer(true, solution, table.matrixData[0][n_var]);

}

bool compareDoubleVectors(const std::vector<double>& vec1, const std::vector<double>& vec2, double epsilon) {
    if (vec1.size() != vec2.size()) {
        return false;
    }

    for (size_t i = 0; i < vec1.size(); ++i) {
        if (std::fabs(vec1[i] - vec2[i]) >= epsilon) {
            return false;
        }
    }

    return true;
}

bool compareDouble(double a, double b, double epsilon) {
    return std::fabs(a - b) < epsilon;
}

const double INF = 1e9;
const double EPS = 1e-9;

struct Simplex {
    int n, m;
    vector<vector<double>> A;
    vector<double> b, c, x;
    vector<int> B, N;
    double z;

    Simplex(const vector<vector<double>> &A_, const vector<double> &b_, const vector<double> &c_)
        : n(c_.size()), m(b_.size()), A(m, vector<double>(n + 1)), b(b_), c(n + 1), x(n), B(m), N(n), z(0) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                A[i][j] = A_[i][j];
            }
            A[i][n] = b[i]; // Добавляем столбец правых частей
            B[i] = n + i;   // Индексы базисных переменных
        }
        for (int j = 0; j < n; ++j) {
            c[j] = c_[j];
            N[j] = j;  // Индексы небазисных переменных
        }
        c[n] = 0; // Свободный член в целевой функции
    }

    // Поворотная операция (пивот)
    void pivot(int r, int s) {
        double inv = 1.0 / A[r][s];
        for (int j = 0; j <= n; ++j) {
            if (j != s) A[r][j] *= inv;
        }
        A[r][s] = inv;
        for (int i = 0; i < m; ++i) {
            if (i != r && fabs(A[i][s]) > EPS) {
                double coef = A[i][s];
                for (int j = 0; j <= n; ++j) {
                    if (j != s) A[i][j] -= coef * A[r][j];
                }
                A[i][s] = 0;
            }
        }
        double coef = c[s];
        for (int j = 0; j <= n; ++j) {
            if (j != s) c[j] -= coef * A[r][j];
        }
        c[s] = 0;
        z += coef * A[r][n];
        swap(B[r], N[s]);
    }

    // Решение задачи линейного программирования
    double solve() {
        while (true) {
            int s = -1;
            for (int j = 0; j < n; ++j) {
                if (c[j] > EPS) { // Находим положительный коэффициент в целевой функции
                    s = j;
                    break;
                }
            }
            if (s == -1) break; // Если положительных коэффициентов нет, решение найдено
            int r = -1;
            for (int i = 0; i < m; ++i) {
                if (A[i][s] > EPS) { // Проверяем допустимость шагов
                    if (r == -1 || A[i][n] / A[i][s] < A[r][n] / A[r][s]) {
                        r = i;
                    }
                }
            }
            if (r == -1) return INF; // Если нет допустимых шагов, задача не ограничена
            pivot(r, s);
        }
        x = vector<double>(n);
        for (int i = 0; i < m; ++i) {
            if (B[i] < n) x[B[i]] = A[i][n];
        }
        return z;
    }
};

struct InteriorPoint {
    int n, m;
    Matrix A;
    Matrix C, xInitial;
    vector<double> alpha;
    double accuracy;

    InteriorPoint(const Matrix& A_, const Matrix& B_, const Matrix& C_, double accuracy)
            : n(C_.getColumns()), m(B_.getColumns()), A(m, n + m), C(1, n), xInitial(n + m, 1) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                A.matrixData[i][j] = A_.matrixData[i][j];
            }
            A.matrixData[i][n + i] = 1;
        }
        for (int j = 0; j < n; ++j) {
            C.matrixData[0][j] = C_.matrixData[0][j];
        }
        this->accuracy = accuracy;
        alpha = {0.5, 0.9};
        
        for (int i = 0; i < n; i++) {
            xInitial.matrixData[i][0] = 1;
        }
        
        for (int i = 0; i < m; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                sum += A.matrixData[i][j] * xInitial.matrixData[j][0];
            }
            xInitial.matrixData[n + i][0] = B_.matrixData[0][i] - sum;
        }
    }
    
    vector<Matrix> solve() {
        vector<Matrix> ans{Matrix(1, n + m), Matrix(1, n + m)};
        IdentityMatrix I(n + m);
        Matrix secondInitial = xInitial;
        Matrix xPrev = xInitial;

        int iterations = 0;

        for (int i = 0; i < 2; i++) {
            while (true) {
                if (iterations > ITERATION_LIMIT) {
                    Matrix inf(1, 1);
                    inf.matrixData[0][0] = 1e9;
                    return {inf};
                    break;
                }

                iterations++;


                Matrix D = *(new Matrix(n + m, n + m));
                D.setDiagonal(xPrev);

                Matrix AA = *(A * D);
                Matrix cVec = Matrix(n + m, 1);
                
                for (int j = 0; j < n; j++) {
                    cVec.matrixData[j][0] = 1;
                }
                
                Matrix cc = *(D * cVec);
                Matrix AA_AAt_inv = *(*(AA * *AA.transpose())).inverse();
                Matrix AAt_AA_AAt_inv = *(*(AA.transpose())*AA_AAt_inv);
                Matrix AAt_AA_AAt_inv_AA = *(AAt_AA_AAt_inv*AA);
                Matrix P = *(I - AAt_AA_AAt_inv_AA);
                Matrix cp = *(P * cc);

                double alp = alpha[i];
                double v = getLargestAbsNegativeValue(cp);
                double mult = alp / v;
                cp = *multiplyCoef(cp, mult);
                Matrix ones(n+m, 1);
                for (int i = 0; i < n+m; i++) {
                    ones.matrixData[i][0] = 1;
                }
                Matrix xx = *(ones + cp);
                xInitial = *(D*xx);
                if (normDiff(xInitial, xPrev) < accuracy) {
                    ans[i] = xInitial;
                    break;
                }
                xPrev = xInitial;
            }
            xPrev = secondInitial;
            xInitial = secondInitial;
        }
        return ans;
    }

    double getLargestAbsNegativeValue(Matrix &A) {
        double minVal = 0;
        for (int i = 0; i < A.getRows(); i++) {
            minVal = min(minVal, A.matrixData[i][0]);
        }
        return abs(minVal);
    }

    Matrix* multiplyCoef(Matrix &A, double coef) {
        for (int i = 0; i < A.getRows(); i++) {
            for (int j = 0; j < A.getColumns(); j++) {
                A.matrixData[i][j] *= coef;
            }
        }
        return &A;
    }
    
    double normDiff(const Matrix &A, const Matrix &B) {
        double sum = 0;
        for (int i = 0; i < A.getRows(); i++) {
            sum += pow(A.matrixData[i][0] - B.matrixData[i][0], 2);
        }
        return sqrt(sum);
    }

    Matrix* addMatrices(const Matrix &A, const Matrix &B) {
        Matrix* result = new Matrix(A.getRows(), A.getColumns());
        for (int i = 0; i < A.getRows(); i++) {
            for (int j = 0; j < A.getColumns(); j++) {
                result->matrixData[i][j] = A.matrixData[i][j] + B.matrixData[i][j];
            }
        }
        return result;
    }

    void print() {
        cout << A;
        for (int i = 0; i < n + m; i++) {
            cout << xInitial;
        }
    }
};

/*
2 2
1 1
2 4
1 3
16 9
0.0001
*/

void runTests() {
//    Test 1

    Matrix C(1, 3);
    Matrix A(3, 3);
    Matrix b(3, 1);

    C.matrixData[0] = vector<double> {9, 10, 16};
    A.matrixData[0] = vector<double> {18, 15, 12};
    A.matrixData[1] = vector<double> {6, 4, 8};
    A.matrixData[2] = vector<double> {5, 3, 3};
    b.matrixData[0][0] = 360;
    b.matrixData[1][0] = 192;
    b.matrixData[2][0] = 180;

    cout << "Test 1:" << endl;

    Answer ans1 = solveLLP(C, A, b, 3);

    InteriorPoint interior1(A, *b.transpose(), C, 0.00001);
    vector<Matrix> p1 = interior1.solve();
    cout << "Alpha 0.5 x*: " << *p1[0].transpose();
    double ans1_z = 0;
    for (int i = 0; i < 3; i++) {
        ans1_z += C.matrixData[0][i] * p1[0].matrixData[i][0];
    }
    cout << "Z: " << ans1_z << endl;


    cout << "Alpha 0.9 x*: " << *p1[1].transpose();
    ans1_z = 0;
    for (int i = 0; i < 3; i++) {
        ans1_z += C.matrixData[0][i] * p1[1].matrixData[i][0];
    }
    cout << "Z: " << ans1_z << endl;

    cout << "Simplex:" << endl;
    cout << "x: " << ans1.solution;
    cout << "Z: " << ans1.z << endl << endl;

//    Test 2
    C = *(new Matrix(1, 5));
    A = *(new Matrix(3, 5));
    b = *(new Matrix(3, 1));

    C.matrixData[0] = vector<double> {1, 1, 0, 0, 0};
    A.matrixData[0] = vector<double> {2, 4, 1, 0, 0};
    A.matrixData[1] = vector<double> {-4, 2, 0, 1, 0};
    A.matrixData[2] = vector<double> {1, 3, 0, 0, 1};
    b.matrixData[0][0] = 16;
    b.matrixData[1][0] = 8;
    b.matrixData[2][0] = 9;
    
    cout << "Test 2:" << endl;

    Answer ans2 = solveLLP(C, A, b, 4);

    InteriorPoint interior2(A, *b.transpose(), C, 0.00001);
    vector<Matrix> p2 = interior2.solve();
    cout << "Alpha 0.5 x*: " << *p2[0].transpose();
    double ans2_z = 0;
    for (int i = 0; i < 3; i++) {
        ans2_z += C.matrixData[0][i] * p2[0].matrixData[i][0];
    }
    cout << "Z: " << ans2_z << endl;


    cout << "Alpha 0.9 x*: " << *p2[1].transpose();
    ans2_z = 0;
    for (int i = 0; i < 3; i++) {
        ans2_z += C.matrixData[0][i] * p2[1].matrixData[i][0];
    }
    cout << "Z: " << ans2_z << endl;

    cout << "Simplex:" << endl;
    cout << "x: " << ans2.solution;
    cout << "Z: " << ans2.z << endl << endl;

//    Test 3

    C.matrixData[0] = vector<double> {1, 1, 0, 0, 0};
    A.matrixData[0] = vector<double> {-2, -2, 1, 0, 0};
    A.matrixData[1] = vector<double> {-5, 3, 0, 1, 0};
    A.matrixData[2] = vector<double> {4, 6, 0, 0, 1};
    b.matrixData[0][0] = -14;
    b.matrixData[1][0] = 15;
    b.matrixData[2][0] = 24;

    cout << "Test 3:\n";

    Answer ans3 = solveLLP(C, A, b, 4);

    InteriorPoint interior3(A, *b.transpose(), C, 0.00001);
    vector<Matrix> p3 = interior3.solve();
    
    if (p3[0].matrixData[0][0] == 1e9) {
        cout << "Alpha 0.5 x*: The method is not applicable!" << endl;
        cout << "Alpha 0.9 x*: The method is not applicable!" << endl;
    } else {
        cout << "Alpha 0.5 x*: " << *p3[0].transpose();
        double ans3_z = 0;
        for (int i = 0; i < 3; i++) {
            ans3_z += C.matrixData[0][i] * p3[0].matrixData[i][0];
        }
        cout << "Z: " << ans3_z << endl;


        cout << "Alpha 0.9 x*: " << *p3[1].transpose();
        ans3_z = 0;
        for (int i = 0; i < 3; i++) {
            ans3_z += C.matrixData[0][i] * p3[1].matrixData[i][0];
        }
        cout << "Z: " << ans3_z << endl;
    }
    cout << "Simplex:" << endl;

    if (!ans3.solver_sate) {
        cout << "The method is not applicable!\n";
    } else {
        cout << "x: " << ans3.solution;
        cout << "Z: " << ans3.z << endl << endl;
    }



}


int main() {
    runTests();

    cout << "Input LLP:" << endl;
    
    
    int n_var,n_constr;
    cin >> n_var >> n_constr;
    Matrix C = *(new Matrix(1, n_var));
    Matrix A = *(new Matrix(n_constr, n_var));
    Matrix B = *(new Matrix(1, n_constr));
    double appr;
    cin >> C;
    cin >> A;
    cin >> B;
    cin >> appr;

    InteriorPoint interior(A, B, C, appr);
    vector<Matrix> p = interior.solve();

    cout << "Simplex:" << endl;
    Answer ans = solveLLP(C, A, *B.transpose(), 2);
    if (ans.solver_sate) {
         cout << "Solved: " << ans.solver_sate << endl;
         cout << ans.solution;
         cout << ans.z << endl;
    } else {
         cout << "The method is not applicable!" << endl;
    }

    cout << "Interior Point:\n";
    if (p[0].matrixData[0][0] == 1e9) {
        cout << "Alpha 0.5 x*: The method is not applicable!" << endl;
        cout << "Alpha 0.9 x*: The method is not applicable!" << endl;
    } else {
        cout << "Alpha 0.5 x*: " << *p[0].transpose();
        double ans3_z = 0;
        for (int i = 0; i < 3; i++) {
            ans3_z += C.matrixData[0][i] * p[0].matrixData[i][0];
        }
        cout << "Z: " << ans3_z << endl;


        cout << "Alpha 0.9 x*: " << *p[1].transpose();
        ans3_z = 0;
        for (int i = 0; i < 3; i++) {
            ans3_z += C.matrixData[0][i] * p[1].matrixData[i][0];
        }
        cout << "Z: " << ans3_z << endl;
    }
}
