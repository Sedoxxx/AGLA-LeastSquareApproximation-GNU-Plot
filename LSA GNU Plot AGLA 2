//Mohamed Aymen Daassi CS-04
//m.daassi@innopolis.university
#include "bits/stdc++.h"
using namespace std;
typedef string str;
typedef long long ll;
#define int ll
typedef double db;
typedef long double ld;
typedef pair<int,int> pi;
typedef vector<int> vi;
typedef vector<ld> vd;
typedef vector<str> vs;
typedef vector<pi> vpi;
#define pb push_back
#define endl "\n"
#define boost ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0)
const double EPSILON = 1e-10;
#ifndef LOCAL
#define cerr if(false) cerr
#endif

#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#endif

class Matrix {
public:
    int n,m;
    vector<vector<double>> mat;
    Matrix (int nn, int mm) {
        n = nn;
        m = mm;
    }

    Matrix operator+ (const Matrix &x){
        Matrix matri(n,m);
        for (int i = 0; i < n; i++){
            vector<double> v;
            for (int j = 0; j < m; j++){
                v.pb(mat[i][j] + x.mat[i][j]);
            }
            matri.mat.pb(v);
        }
        return matri;
    }

    Matrix operator* (const Matrix x){
        Matrix matri(n,x.m);
        for (int i = 0; i < n; i++){
            vector<double> v;
            for (int j = 0; j < x.m; j++){
                v.pb(0);
            }
            matri.mat.pb(v);
        }
        for (int i = 0; i < n; i++){
            for (int j = 0; j < x.m; j++){
                for (int k = 0; k < m; k++){
                    matri.mat[i][j] += mat[i][k] * x.mat[k][j];
                }
            }
        }
        return matri;
    }

    Matrix operator- (const Matrix x){
        Matrix matri(n,m);
        for (int i = 0; i < n; i++){
            vector<double> v;
            for (int j = 0; j < m; j++){
                v.pb(mat[i][j] - x.mat[i][j]);
            }
            matri.mat.pb(v);
        }
        return matri;
    }

    Matrix operator=(const Matrix x) {
        Matrix matri(x.n , x.m);
        for (int i = 0; i < x.n; i++){
            for (int j = 0; j < x.m; j++){
                matri.mat[i][j] = x.mat[i][j];
            }
        }
        return matri;
    }

    Matrix transpose (){
        Matrix matri(m,n);
        for (int i = 0; i < m; i++){
            vector<double> v;
            for (int j = 0; j < n; j++){
                v.pb(mat[j][i]);
            }
            matri.mat.pb(v);
        }
        return matri;
    }

    void Normalize (Matrix *x){
        for (int i = 0; i < x->n; i++) {
            double divisor = x->mat[i][i];
            for (int j = 0; j < x->n; j++) {
                mat[i][j] /= divisor;
            }
        }
    }

    void InitZero(){
        for (int i = 0; i < n;i++){
            vector<double> v;
            for (int j =0 ;j < m ; j++){
                v.pb(0.00);
            }
            mat.pb(v);
        }
    }

    friend ostream & operator << (ostream &out, const Matrix &A){
        for (int i = 0; i < A.n; i++){
            for (int j = 0; j < A.m; j++){
                double x = A.mat[i][j];
                if (fabs(x) < EPSILON){
                    out << "0.0000 " ;
                }
                else {
                    out << fixed << setprecision(4) <<  x <<" ";
                }
            }
            out << endl;
        }
        return out;
    }

    friend istream & operator >> (istream &in, Matrix *A){
        for (int i = 0; i < A->n; i++){
            vector<double> v;
            for (int j = 0; j <A->m; j++){
                double x; in>>x;
                v.pb(x);
            }
            A->mat.pb(v);
        }
        return in;
    }
    int getN (){
        return n;
    }
    int getM (){
        return m;
    }
};
class SquareMatrix : public Matrix {
public :
    SquareMatrix(int nn) : Matrix (nn , nn){
        n = nn;
        m = nn;
    }
};
class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int nn) : SquareMatrix(nn){
        n = nn;
        for (int i = 0; i < n; i++){
            vector<double> v;
            for (int j = 0; j < n; j++){
                if (i != j) v.pb(0);
                else v.pb(1);
            }
            mat.pb(v);
        }
    }

};
class EliminationMatrix : public IdentityMatrix{
public :
    EliminationMatrix (int nn , Matrix *x, int row, int col) : IdentityMatrix(nn){
        mat[row][col] = - x->mat[row][col] / x->mat[col][col];
    }
};

class PermutationMatrix : public IdentityMatrix {
public :
    PermutationMatrix (int nn, int i, int j) : IdentityMatrix(nn){
        for (int k = 0; k < nn; k++)
            swap(mat[i][k],mat[j][k]);
    }
};

class ColumnMatrix : public Matrix {
public :
    ColumnMatrix (int nn) : Matrix (nn, 1){}
};

class DataSetMatrix : public Matrix {
public:
    DataSetMatrix (int nn) : Matrix (nn, 2){
    }
};
class AugmentedMatrix : public Matrix {
public:
    AugmentedMatrix (Matrix *A, Matrix *B, int nn) : Matrix(nn, 2*nn){
        InitZero();
        for (int i = 0; i < n; i++){
            for (int j = 0; j < nn; j++){
                mat[i][j] = A->mat[i][j];
                mat[i][j+nn] = B->mat[i][j];
            }
        }
    }
};

class MatrixComputations : Matrix {
public :

    static int findPivotRow(int col, Matrix *x) {
        int pivotRow = col;
        for (int i = col + 1; i < x->mat.size(); i++) {
            if (fabs(x->mat[i][col]) > fabs(x->mat[pivotRow][col])) {
                pivotRow = i;
            }
        }
        return pivotRow;
    }
    static Matrix InverseMatrix (Matrix *x){
        Matrix *I = new IdentityMatrix(x->n);
        for (int col = 0; col < x->n; col++) {
            int pivotRow = findPivotRow(col, x);
            if (col != pivotRow){
                Matrix *P = new PermutationMatrix(x->n, pivotRow, col);
                Matrix A = *P * *x ;
                Matrix B = *P * *I ;
                x->mat = A.mat;
                I->mat = B.mat;
            }
            for (int row = col + 1; row < x->n; row++) {
                if (x-> mat[row][col] != 0.00) {
                    Matrix *E = new EliminationMatrix(x->n, x, row, col);
                    Matrix A = *E * *x;
                    Matrix B = *E * *I ;
                    x->mat = A.mat;
                    I->mat = B.mat;
                }
            }
        }
        for (int col = x->n - 1; col >= 0; col--) {
            for (int row = col -1 ; row >= 0; row--) {
                Matrix *E = new EliminationMatrix(x->n, x, row, col);
                Matrix A = *E * *x;
                Matrix B = *E * *I ;
                x->mat = A.mat;
                I->mat = B.mat;
            }
        }
        Matrix *E = new IdentityMatrix(x->n);
        E->Normalize(x);
        Matrix A = *E * *x;
        Matrix B = *E * *I;
        x->mat = A.mat;
        I->mat = B.mat;
        return *I;
    }

    static Matrix leastSquareApprox (Matrix *D, int poly){
        Matrix A(D->n, poly+1);
        A.InitZero();
        for (int i = 0; i < D->n ; i++){
            A.mat[i][0] = 1.0;
            for (int j = 1; j <= poly; j++){
                A.mat[i][j] = A.mat[i][j-1] * D->mat[i][0];
            }
        }
        Matrix A_T = A.transpose();
        cout <<"A:"<< endl;
        cout << A ;
        cout <<"A_T*A:"<< endl;
        Matrix B = A_T * A;
        cout << B;
        cout <<"(A_T*A)^-1:"<< endl;
        Matrix Inv = MatrixComputations::InverseMatrix(&B);
        cout << Inv;
        Matrix * b = new ColumnMatrix(D->n);
        b->InitZero();
        for (int i = 0; i < D->n ; i++){
            b->mat[i][0] = D->mat[i][1];
        }
        cout <<"A_T*b:"<< endl;
        cout<< A_T * *b;
        cout <<"x~:"<<endl;
        Matrix result = Inv * A_T * *b;
        cout << result;
        return result;
    }
};
int32_t main(){
    boost;
#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#endif
    int m;
    cin >> m;
    Matrix* D = new DataSetMatrix(m);
    cin >> D;
    int n;
    cin >> n;
    Matrix result = MatrixComputations::leastSquareApprox(D,n);
    if (pipe != NULL){
        double cof3 = result.mat[0][0];
        double cof2 = result.mat[1][0];
        double cof1 = result.mat[2][0];
        double cof0 = result.mat[3][0];
        fprintf(pipe, "plot [0 : 25] [0 : 25] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", cof0, cof1, cof2, cof3);
        for (int i = 0; i < m ; i++){
            double x = D->mat[i][0];
            double y = D->mat[i][1];
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe, "%s\n", "e");
        fflush(pipe);
#ifdef WIN32
        _pclose(pipe);
#endif
    }
    else{
        cout << "Could not open pipe" << endl;
    }
    return 0;
}
//Never stop eating nutella!
