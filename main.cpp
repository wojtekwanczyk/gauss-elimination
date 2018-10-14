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
double **generateMatrix(double **t, int n, int m);
double *readArray(double *p, int n);
double *multiplyMatrix(double **A1, int n, int m, double *A2);
double *getSolutionFromArray(double **A, int n, int m);
void showDifferenceEu(double *t1, double *t2, int n);
void showDifferenceMax(double *t1, double *t2, int n);

double getPrecision(double value, double precision)
{
    return (floor((value * pow(10, precision) + 0.5)) / pow(10, precision));
}

const int prec = 15;     //significant places
const double eps = 1 / pow(10, prec);


int main() {

    // define size of the matrix
    int n,m;
    cout << "Precision = " << eps << endl;
    cout << "A - Height <space> Width" << endl;
    cin >> n >> m;
    //n = 3;
    //m = 3;

    // memory allocation
    double **A = new double*[n];
    for(int i = 0; i<n ; i++){
        A[i] = new double[m+1];
    }
    double *x = new double[m];
    double *B = new double[m];


    // Multiplying
    readMatrix(A, n, m);
    cout << "Array x to multiply: " << endl;
    x = readArray(x, m);
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
        A[i][m] = B[i];
    }

    showMatrix(A, n, m+1, 1);
    double *solution = gaussElimination(A, n, m+1);

    if(solution != NULL) showDifferenceEu(x, solution, m);
    if(solution != NULL) showDifferenceMax(x, solution, m);

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


double **generateMatrix(double **t, int n, int m){

    srand(time(NULL));

    cout << "Matrix " << n << " x " << m << " (with spaces):" << endl;
    for(int i=0; i<n ; i++){
        for(int j=0; j<m; j++){
            t[i][j] = rand()%20;
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
        cout <<t[i] << "\t\t";
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

            if( abs(t[factor_line][factor_column]) < eps ){         // == 0
                // search for non zero element
                int l;
                for(l=k; l<n; l++){
                    if( abs(t[l][factor_column]) > 0) break;        // != 0
                }

                if(l<n){
                    //element was found
                    double *tmp = t[factor_line];
                    t[factor_line] = t[l];
                    t[l] = tmp;
                    cout << "Line change!\n";
                    showMatrix(t, n, m, 1);
                } else {
                    cout << "NOT FOUND";
                    //range--;
                    factor_column++;
                    continue;
                }
            }
            factor = t[k][factor_column] / t[factor_line][factor_column];
            factor = getPrecision(factor, prec);


            // j - subtract every column in k line
            int sum = 0;
            for(int j=factor_column; j<m; j++){
                t[k][j] = t[k][j] - factor * t[factor_line][j];
                t[k][j] = getPrecision(t[k][j], prec);
                sum += t[k][j];
            }
            if(abs(sum - t[k][m-1]) < eps) range--;
            if(abs(sum - t[k][m-1]) < eps && abs(t[k][m-1]) > eps ) conflict++;
            if(abs(sum - t[k][m-1]) < eps && abs(t[k][m-1]) < eps ) identity++;
            showMatrix(t, n ,m, 1);
        }
        //if(t[factor_line][factor_column] < eps) range--;    // == 0
        factor_column++;
        factor_line++;
    }

    cout << "Rzad macierzy glownej wynosi " << range << endl;
    if(conflict != 0) cout << "Uklad sprzeczny!!!" << endl;
    else if(identity != 0) cout << "Uklad tozsamosciowy!!!" << endl;
    else {
        double *solution = getSolutionFromArray(t, n, m);
        cout << "=== ROZWIAZANIE ==========================" << endl;
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
            p[i] = getPrecision(p[i], prec);
        }
    }

    return p;
}


double *readArray(double *p, int n){
    for(int i=0; i<n ; i++){
        cin >> p[i];
    }
    cout << endl;
    return p;
}


double *getSolutionFromArray(double **A, int n, int m){
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

