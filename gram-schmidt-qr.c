#include <math.h>
#include <stdio.h>

#define EPSILON 0.000001

void print_matrix(double* M, int m, int n)
{
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      printf(" %7f ", M[i*n+j]);
    }
    printf("\n");
  }
}

double * trans(double* At, double* A, int m, int n)
{
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      At[i*m+j] = A[j*m+i];
    }
  }
  return A;
}

void mult(double* M, double* M1, double* M2, int m1, int n1, int n2)
{
  for (int r1=0; r1<m1; r1++) {
    for (int c2=0; c2<n2; c2++) {
      M[r1*m1+c2] = 0;
      for (int i=0; i<n1; i++) {
        M[r1*m1+c2] += M1[r1*m1+i]*M2[i*m1+c2];
      }
    }
  }
}

double col_norm(double* v, int m, int n,int col)
{
  double sum = 0.0;
  for (int i=0; i<m; i++) {
    sum += v[i*n+col]*v[i*n+col];
  }

  return sqrt(sum);
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

void col_copy(double* M, double *M0, int m, int n, int col)
{
  for (int r=0; r<m; r++) {
    M[r*n + col] = M0[r*n + col];
  }
}

void col_mul_sub(double* M1, double* M2, int m, int n, int col1, int col2, double v)
{
  for (int r=0; r<m; r++) {
    M1[r*n + col1] -= M2[r*n + col2]*v;
  }
}

void col_mult(double* M, int m, int n, int col, double v)
{
  for (int r=0; r<m; r++) {
    M[r*n + col] *= v;
  }
}

void col_div(double* M, int m, int n, int col, double v)
{
  for (int r=0; r<m; r++) {
    M[r*n + col] /= v;
  }
}

double col_dot_mult(double* M1, double * M2, int m, int n, int col1, int col2)
{
  double sum = 0.0;
  for (int r=0; r<m; r++) {
    sum += M1[r*n + col1] * M2[r*n +col2];
  }
  return sum;
}

void qr(double* A, int m, int n, double* Q, double *R)
{
  double c_jj = 1.0;
  double c_ij = 1.0;

  for (int j=0; j<n; j++) {
    col_copy(Q, A, m, n, j);

    for(int i=0; i<j; i++) {
      c_ij = col_dot_mult(Q, Q, m, n, i, j);
      col_mul_sub(Q, Q, m, n, j, i, c_ij);

      R[i*n+j] = c_ij; 
    }
    
    c_jj = col_norm(Q, m, n, j);
    col_div(Q, m, n, j, c_jj);

    R[j*n+j] = c_jj; 
  }
}

double* matrix() {
  static double a[3][2];

  a[0][0] = 1; 
  a[0][1] = -4;
  a[1][0] = 2; 
  a[1][1] = 3;
  a[2][0] = 2; 
  a[2][1] = 2;

  return a;
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
  double *Q = buffer(1);
  double *R = buffer(2);

  printf("\n----- The A matrix -----\n");
  print_matrix(matrix(), 3, 2);

  qr(matrix(), 3, 2, Q, R);

  printf("\n----- The Q matrix -----\n");
  print_matrix(Q, 3, 2);

  printf("\n----- The R matrix -----\n");
  print_matrix(R, 2, 2);

  return 0;
}
