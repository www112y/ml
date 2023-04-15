#include <math.h>
#include <stdio.h>

#define EPSILON 0.000001

typedef double (* matrix_fptr)(int m, int n);
typedef double (* b_fptr)(int n);

void print_matrix(double* M, int m, int n)
{
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      printf(" %7f ", M[i*m+j]);
    }
    printf("\n");
  }
}

void print_x(double x[], int n)
{
  printf("\n  ( ");

  for (int l=0; l<n; l++) {
    printf(" %f, ", x[l]);
  }
  printf(")\n");
}

double * trans(double* A, matrix_fptr M, int m, int n)
{
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      A[i*m+j] = M(j, i);
    }
  }
  return A;
}

void mult(double* M, matrix_fptr M1, matrix_fptr M2, int m1, int n1, int n2)
{
  for (int row1=0; row1<m1; row1++) {
    for (int col2=0; col2<n2; col2++) {
      M[row1*m1+col2] = 0;
      for (int i=0; i<n1; i++) {
        M[row1*m1+col2] += M1(row1,i)*M2(i,col2);
      }
    }
  }
}

double * buffer(int slot)
{
  static double x1[100];
  static double x2[100];
  static double x3[100];

  int i = 0;
  switch(slot) {
    case 1:
      for(i=0; i<100; i++)
        x1[i] = 0;
      return x1; 
    case 2:
      for(i=0; i<100; i++)
        x2[i] = 0;
      return x2; 
    case 3:
      for(i=0; i<100; i++)
        x3[i] = 0;
      return x3;
    default:
      for(i=0; i<100; i++)
        x1[i] = 0;
      return x1; 
  }
}

double * buffer_matrix(matrix_fptr A1, double* A2, int m, int n)
{
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      A2[i*m+j] = A1(i, j);
    }
  }

  return A2;
}

/* Ax = b*/
void root(double * A, double * b, int n, double x[])
{
  /* A = LU*/
  /* Ux = L^(-1)b */
  int i = 0;
  int j = 0;
  int k = 0;
  for (i=0; i<n; i++) {
    for (j=i+1; j<n; j++) {
      double c = A[j*n+i]/A[i*n+i];

      for (k=i+1; k<n; k++) {
        //A[j][k] = A[j][k] - A[j][k]*A[i][i]/A[j][i];
	A[j*n+k] = -1.0*c*A[i*n+k] + A[j*n+k];
      }
      //b[j] -= b[j]*A[i][i]/A[j][i];
      b[j*n] = -1.0*c*b[i*n] + b[j*n];

      //A[j][i]=0;
      A[j*n+i]=0;
    }
  }

  /* x = U^(-1)L^(-1)b */
  for (i=n-1; i>=0; i--) {
    x[i] = b[i*n]/A[i*n+i];
 
    for (j=n-1; j>i; j--) {
      x[i] -= x[j]*A[i*n+j]/A[i*n+i];
    }
  }
}

/* LS */
void ls_root(matrix_fptr A0, matrix_fptr At, matrix_fptr b0, int m, int n, double x[])
{
  double * A = buffer(1);
  mult(A, At, A0, n, m, n);

  printf("----- The Matrix A -----\n");
  print_matrix(buffer_matrix(A0, buffer(3), m, n), m, n);

  double * b = buffer(2);
  mult(b, At, b0, n, m, 1);

  printf("----- The vector b -----\n");
  print_matrix(buffer_matrix(b0, buffer(3), m, 1), m, 1);

  root(A, b, n, x);
}

double matrix(int m, int n) {
  static double a[3][2];

  a[0][0] = 1; 
  a[0][1] = 1;
  a[1][0] = 1; 
  a[1][1] = -1;
  a[2][0] = 1; 
  a[2][1] = 1;

  return  a[m][n];
}

double trans_matrix(int m, int n) {
  static double a[2][3];

  a[0][0] = 1;
  a[1][0] = 1;
  a[0][1] = 1;
  a[1][1] = -1;
  a[0][2] = 1;
  a[1][2] = 1;

  return  a[m][n];
}

double b(int m, int n) {
  static double v[3];

  (void) n;

  v[0] = 2;
  v[1] = 1;
  v[2] = 3;

  return v[m];
}

int main()
{
  double x[3];
  ls_root(matrix, trans_matrix, b, 3, 2, x);

  printf("----- The fitting solution -----\n");
  printf(" (%f,  %f)\n", x[0], x[1]);
  return 0;
}
