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
    static double Norm(Matrix x){
        double sum = 0;
        for (int i = 0; i < x.n; i++){
            sum+= pow(x.mat[i][0],2);
        }
        return sqrt(sum);
    }
    friend ostream & operator << (ostream &out, const Matrix &A){
        for (int i = 0; i < A.n; i++){
            for (int j = 0; j < A.m; j++){
                double x = A.mat[i][j];
                if (fabs(x) < EPSILON){
                    out << "0.0000 " ;
                }
                else {
                    out << fixed << setprecision(4) << x <<" ";
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

class DiagonalMatrix : public SquareMatrix {
public :
    DiagonalMatrix (Matrix *A, int nn) : SquareMatrix(nn){
        InitZero();
        for (int i = 0; i < n; i++){
            mat[i][i] = A->mat[i][i];
        }
    }
};

class UpperTriangMatrix : public SquareMatrix {
public :
    UpperTriangMatrix (Matrix *A, int nn) : SquareMatrix(nn){
        InitZero();
        for (int i = 0; i < nn; i++){
            for (int j = 0; j < nn; j++){
                if (j > i)
                    mat[i][j] = A -> mat[i][j];
            }
        }
    }
};
class LowerTriangMatrix : public SquareMatrix {
public :
    LowerTriangMatrix (Matrix *A, int nn) : SquareMatrix(nn){
        InitZero();
        for (int i = 0; i < nn; i++){
            for (int j = 0; j < nn; j++){
                if (i >= j)
                    mat[i][j] = A -> mat[i][j];
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

    static void leastSquareApprox (Matrix *D, int poly){
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
        cout << Inv * A_T * *b ;
    }
    static bool Dominant (Matrix *A){
        for (int i = 0 ; i < A->n ; i++){
            double sum = 0;
            for (int j = 0; j < A->n; j++){
                if (i!=j) sum += A->mat[i][j];
            }
            if (fabs(A->mat[i][i]) < fabs(sum))
                return false;
        }
        return true;
    }
    static void jacobi (Matrix *A, Matrix *b, double acc){
        Matrix* D = new DiagonalMatrix(A, A->n);
        Matrix* I = new IdentityMatrix(A->n);
        Matrix D_1 = InverseMatrix(D);
        Matrix alpha = *I - (D_1 * *A);
        cout << "alpha:" << endl;
        cout << alpha;
        Matrix beta = D_1 * *b;
        cout << "beta:" << endl;
        cout << beta;
        Matrix x0 = beta;
        cout << "x(0):" << endl;
        cout << x0;
        Matrix xi = x0 + D_1 * (*b - *A * x0);
        double eps = Norm(xi - x0);
        cout << "e: " << eps << endl;
        cout << "x(1):" << endl;
        cout << xi;
        int count = 1;
        while (eps > acc){
            Matrix xiprev = xi;
            Matrix ITER =  xi + D_1 * (*b - *A * xi);
            xi.mat = ITER.mat;
            count++;
            eps = Norm(xi - xiprev);
            cout << "e: " << eps << endl;
            cout << "x("<<count<<"):" << endl;
            cout << xi;
        }
    }

    static void seidel (Matrix *A, Matrix *b, double acc){
        Matrix* D = new DiagonalMatrix(A, A->n);
        Matrix* I = new IdentityMatrix(A->n);
        Matrix D_1 = InverseMatrix(D);
        Matrix beta = D_1 * *b;
        cout << "beta:" << endl;
        cout << beta;
        Matrix alpha = *I - (D_1 * *A);
        cout << "alpha:" << endl;
        cout << alpha;
        Matrix* B = new LowerTriangMatrix(&alpha, A->n);
        Matrix* C = new UpperTriangMatrix(&alpha, A->n);
        Matrix Temp = *B;
        Matrix B_1 = InverseMatrix(&Temp);
        cout << "B:" << endl;
        cout << *B;
        cout << "C:" << endl;
        cout << *C;
        cout << "I-B:" << endl;
        Matrix IMinusB = *I - *B;
        cout << IMinusB;
        Matrix IMinusBInv = MatrixComputations::InverseMatrix(&IMinusB);
        cout<< "(I-B)_-1:" << endl;
        cout << IMinusBInv;
        B = new LowerTriangMatrix(A, A->n);
        C = new UpperTriangMatrix(A, A->n);
        Matrix Temp2 = *B;
        Matrix B_11 = InverseMatrix(&Temp2);
        Matrix x0 = beta;
        cout << "x(0):" << endl;
        cout << x0;
        Matrix xi = B_11 * (*b - *C * x0);
        double eps = Norm(xi - x0);
        cout << "e: " << eps << endl;
        cout << "x(1):" << endl;
        cout << xi;
        int count = 1;
        while (eps > acc){
            Matrix xiPrev = xi;
            Matrix ITER =  B_11 * (*b - *C * xi);
            xi.mat = ITER.mat;
            count++;
            eps = Norm(xi - xiPrev);
            cout << "e: " << eps << endl;
            cout << "x("<<count<<"):" << endl;
            cout << xi;
        }
    }

    static void PredPrey (double v0,double k0,
                          double alpha1, double beta1,
                          double alpha2, double beta2,
                          int time
                          ,int nbPoints,
                          vector<double> &victims,
                          vector<double> &killers,
                          vector<double> &times){
        double increment = 1.00 / (nbPoints / time);
        cout<<"t:"<<endl;
        cout<<fixed<<setprecision(2);
        for (double t = 0.00; t <= time; t+=increment){
            times.pb(t);
            cout<<t<<" ";
        }
        cout<<endl;
        v0 = v0 - (alpha2/beta2);
        k0 = k0 - (alpha1/beta1);
        cout<<"v:"<<endl;
        for (int i = 0; i < times.size(); i++){
            double x = v0*cos(sqrt(alpha1*alpha2)*times[i]);
            double y = ((k0*sqrt(alpha2)*beta1*sin(sqrt(alpha1*alpha2)*times[i]))/(beta2*sqrt(alpha1)));
            cout<< x - y + (alpha2/beta2) << " ";
            victims.pb(x - y + (alpha2/beta2));
        }
        cout<<endl;
        cout<<"k:"<<endl;
        for (int i = 0; i < times.size(); i++){
            double x = ((v0*sqrt(alpha1)*beta2*sin(sqrt(alpha1*alpha2)*times[i]))/(beta1*sqrt(alpha2)));
            double y = k0*cos(sqrt(alpha1*alpha2)*times[i]);
            cout << x + y + (alpha1/beta1) << " ";
            killers.pb(x + y + (alpha1/beta1));
        }
        cout<<endl;
    }
};

int32_t main(){
    boost;
#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#endif
    double v0,k0;
    double alpha1,beta1,alpha2,beta2;
    int time,nbPoints;
    cin>>v0>>k0;
    cin>>alpha1>>beta1>>alpha2>>beta2;
    cin>>time>>nbPoints;
    vector<double> times,victims,killers;
    MatrixComputations::PredPrey(v0,k0,alpha1,beta1,alpha2,beta2,time,nbPoints,victims,killers,times);
    if (pipe != NULL){
        fprintf(pipe, "plot '-' using 1:2 t 'victims' with lines, '-' using 1:2 t 'killers' with lines ");
        /*
         * In case you want a figure in good manual proportions, not automatic by GNU plot for my model:
          150
          50
          0.65
          0.0075
          0.45
          0.0045
          64
          512
          you can use this instead of the previous command
          fprintf(pipe, "plot [0:120] [0:320] '-' using 1:2 t 'victims' with lines, '-' using 1:2 t 'killers' with lines ");
         */
        for (int i = 0; i < times.size(); i++){
            double x = times[i];
            double y = victims[i];
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe, "%s\n", "e");
        for (int i = 0; i < times.size(); i++){
            double x = times[i];
            double y = killers[i];
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe, "%s\n", "e");
        /*
        ###################################################
        Dependency figure implementation Victims to Killers
        ###################################################
        fprintf(pipe, "plot  '-' using 1:2 t 'victims to killers' with lines ");
        for (int i = 0; i < times.size(); i++){
            double x = killers[i];
            double y = victims[i];
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe, "%s\n", "e");*/

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
