#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

using namespace std;

void showMatrix(double **t, int n, int m, int extraSize);
void showArray(double *t, int n);
double **readMatrix(double **t, int n, int m);
double **readArrayToMatrix(double **t, int n, int m);
double *gaussElimination(double **t, int n, int m);
double **generateMatrix1(double **t, int n, int m);
double **generateMatrix2(double **t, int n, int m);
double **generateMatrix3(double **t, int n, int m);
double *readArray(double *p, int n);
double *generateArray(double *p, int n);
double *multiplyMatrix(double **A1, int n, int m, double *A2);
double *getSolutionFromMatrix(double **A, int n, int m);
void showDifferenceEu(double *t1, double *t2, int n);
void showDifferenceMax(double *t1, double *t2, int n);
double *solve3DiagonalMatrix(double **t, int n);


// shift significant places before comma, cut everything after comma
// then shift everything back again
double getPrecision(double value, double precision)
{
    return (floor((value * pow(10, precision) + 0.5)) / pow(10, precision));
}


int n = 5, m;
const int prec = 15;     //significant places
const double eps = 1 / pow(10, prec);


int main() {

    // define size of the matrix
    //int n,m;
    cout << "Precision = " << eps << endl;
    cout << "A - matrix size: " << n <<endl;
    //cin >> n;
    m = n;

    // memory allocation
    double **A = new double*[n];
    for(int i = 0; i<n ; i++){
        A[i] = new double[m+1];
    }
    double *x = new double[m];
    double *B = new double[m];


    // Multiplying
    A = generateMatrix2(A, n, m);
    //readMatrix(A, n, m);
    //x = readArray(x, m);
    x = generateArray(x, m);
    showMatrix(A, n, m, 0);
    cout << "Multiply by (expected result) =============" << endl;
    showArray(x, m);
    cout << "===========================================" << endl;
    B = multiplyMatrix(A, n, m, x);
    cout << "Result:" << endl;
    showArray(B, m);
    cout<<endl;

    // only for testing
    // Resolving
    //A = readMatrix(A, n, m);
    //A = readArrayToMatrix(A, n, m);
    //or
    //generateMatrix(A, n, m);

    // copy the result to full A matrix -- only with multiplying part
    for(int i=0; i<m; i++){
        B[i] = getPrecision(B[i], prec);
        A[i][m] = B[i];
    }
    cout << "Result with precision cut: " << endl;
    showArray(B, m);

    //showMatrix(A, n, m+1, 1);
    double *solution = gaussElimination(A, n, m+1);

    if(solution != NULL) showDifferenceEu(x, solution, m);
    if(solution != NULL) showDifferenceMax(x, solution, m);

    double *solution3 = solve3DiagonalMatrix(A, n);

    showDifferenceEu(x, solution3, n);
    showDifferenceMax(x, solution3, n);

    return 0;
}


double **readArrayToMatrix(double **t, int n, int m){

    cout << "Array of size " << n << " (with spaces):" << endl;
    for(int i=0; i<n ; i++){
        cin >> t[i][m];
    }
    cout << endl;
    return t;
}


double **readMatrix(double **t, int n, int m){

    cout << "Matrix " << n << " x " << m << " (with spaces):" << endl;
    for(int i=0; i<n ; i++){
        for(int j=0; j<m; j++){
            cin >> t[i][j];
        }
    }
    cout << endl;
    return t;
}

// first task
double **generateMatrix1(double **t, int n, int m){

    for(int i=0; i<n ; i++){
        for(int j=0; j<m; j++){
            if(i == 0)  // i == 1, dif notation
                t[i][j] = 1.0;
            else
                t[i][j] = 1.0 / ((i+1) + (j+1) - 1);    // diffrent notation from 0 and from 1
            t[i][j] = getPrecision(t[i][j], prec);
        }
    }
    cout << endl;
    return t;
}


// second task
double **generateMatrix2(double **t, int n, int m){

    for(int i=0; i<n ; i++){
        for(int j=0; j<m; j++){
            if(j >= i)
                t[i][j] = double(2 * (i+1)) / (j+1);
            else
                t[i][j] = t[j][i];
            t[i][j] = getPrecision(t[i][j], prec);
        }
    }
    cout << endl;
    return t;
}

// 3 task
double **generateMatrix3(double **t, int n, int m){

    for(int i=0; i<n ; i++){
        for(int j=0; j<m; j++){
            if(j >= i)
                t[i][j] = double(2 * (i+1)) / (j+1);
            else
                t[i][j] = t[j][i];
            t[i][j] = getPrecision(t[i][j], prec);
        }
    }
    cout << endl;
    return t;
}


void showMatrix(double **t, int n, int m, int extraSize = 0){

    for(int i=0; i<n ; i++){
        for(int j=0; j<m; j++){
            cout << t[i][j] << "\t";
            if(j+2 == m && extraSize == 1) cout<< "|";
            cout << "\t";
        }
        cout << endl;
    }
    cout << endl;
}


void showArray(double *t, int n){

    for(int i=0; i<n ; i++){
        cout << setprecision(17) << t[i] << "\t\t";
    }
    cout << endl;
}


double *gaussElimination(double **t, int n, int m){
    double factor;
    int factor_column = 0;
    int factor_line = 0;
    int range = m-1;
    int conflict = 0;
    int identity = 0;

    // i = step number
    for(int i=1; i<n; i++){

        // k - multiplies every line
        for(int k = i; k<n; k++){

            if(t[factor_line][factor_column] == 0 ){
                // search for non zero element
                int l;
                for(l=k; l<n; l++){
                    if(t[l][factor_column] != 0) break;
                }

                if(l<n){
                    //element was found
                    double *tmp = t[factor_line];
                    t[factor_line] = t[l];
                    t[l] = tmp;
                    showMatrix(t, n, m, 1);
                } else {
                    factor_column++;
                    continue;
                }
            }
            factor = t[k][factor_column] / t[factor_line][factor_column];

            // j - subtract every column in k line
            int sum = 0;
            for(int j=factor_column; j<m; j++){
                t[k][j] = t[k][j] - factor * t[factor_line][j];
                sum += t[k][j];
            }
            if(sum - t[k][m-1] == 0) range--;
            if(sum - t[k][m-1] == 0 && t[k][m-1] != 0 ) conflict++;
            if(sum - t[k][m-1] == 0 && t[k][m-1] == 0) identity++;
        }
        factor_column++;
        factor_line++;
    }

    cout << "Rzad macierzy glownej wynosi " << range << endl;
    if(conflict != 0) cout << "Uklad sprzeczny!!!" << endl;
    else if(identity != 0) cout << "Uklad tozsamosciowy!!!" << endl;
    else {
        double *solution = getSolutionFromMatrix(t, n, m);
        cout << "============== ROZWIAZANIE ===============" << endl;
        showArray(solution, n);
        cout << "==========================================" << endl;
        return solution;
    }
    return NULL;
}


double *multiplyMatrix(double **A1, int n, int m, double *A2){
    double *p = new double[n];

    for(int i=0 ; i<n; i++){
        p[i] = 0;
        for(int j=0; j<m; j++){
            p[i] += A1[i][j] * A2[j];
        }
    }
    return p;
}


double *readArray(double *p, int n){
    cout << "Array x to multiply: " << endl;
    for(int i=0; i<n ; i++){
        cin >> p[i];
    }
    cout << endl;
    return p;
}

double *generateArray(double *p, int n){
    cout << "Array x to multiply: " << endl;
    srand(time(NULL));
    for(int i=0; i<n ; i++){
        int x = rand() % 2;
        if(x == 0)
            p[i] = 1.0;
        else
            p[i] = -1.0;
    }
    cout << endl;
    return p;
}

double *getSolutionFromMatrix(double **A, int n, int m){
    double *solution = new double[n];
    for(int i=0; i<n; i++) solution[i] = 0;

    for(int i=n-1; i>=0; i--){
        solution[i] = A[i][m-1];
        for(int j=m-2; j > i; j--){
            solution[i] -= A[i][j] * solution[j];
        }
        solution[i] /= A[i][i];
    }
    return solution;
}


void showDifferenceEu(double *t1, double *t2, int n){
    double *dif = new double[n];

    cout << endl << "======== Euclides difference: ===========" << endl;
    for(int i=0; i<n; i++){
        dif[i] = abs(t2[i] - t1[i]);
        cout << setprecision(32) << dif[i] << "\t";
    }
    cout << endl << "=========================================" << endl;
}

void showDifferenceMax(double *t1, double *t2, int n){

    cout << endl << "======== Maximum difference: ===========" << endl;
    double maxi = 0;
    double tmp;
    for(int i=0; i<n; i++){
        tmp = abs(t2[i] - t1[i]);
        if(tmp > maxi)
            maxi = tmp;
    }
    cout << setprecision(32) << maxi;
    cout << endl << "=========================================" << endl;
}



double *solve3DiagonalMatrix(double **t, int n){

    double *l = new double[n];
    double *u = new double[n];
    double *a = new double[n];
    double *c = new double[n];

    double *x = new double[n];
    double *y = new double[n];

    u[0] = t[0][0];
    c[0] = t[0][1];
    for(int i=1; i<n-1;i++){
        a[i] = t[i][i-1];
        c[i] = t[i][i+1];

        l[i] = (double)(a[i] / u[i-1]);
        u[i] = t[i][i] - (l[i] * c[i-1]);
    }
    a[n-1] = t[n-1][n-2];
    l[n-1] = (double)(a[n-1] / u[n-2]);
    u[n-1] = t[n-1][n-1] - (l[n-1] * c[n-2]);


    // L * y = b
    // U * x = y

    y[0] = t[0][n];
    for(int i=1; i<n; i++){
        y[i] = (double)(t[i][n] - (l[i] * y[i-1]));
    }


    x[n-1] = (double)(y[n-1] / u[n-1]);
    for(int i=n-2; i>=0 ; i--){
        x[i] = (double)((y[i] - (c[i] * x[i+1])) / u[i]);
    }

    cout << "============== ROZWIAZANIE2 ==============" << endl;
    showArray(x, n);
    cout << "==========================================" << endl;
    return x;
}

// juz dziala ale slabo liczy :/
