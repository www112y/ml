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

void print_vec(double* v, int n)
{
  for (int j=0; j<n; j++) {
    printf(" %7f \n", v[j]);
  }
}

double* trans(double* At, double* A, int m, int n)
{
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      At[j*m+i] = A[i*n+j];
    }
  }
  return A;
}

void mul(double* M1, double* M2, int m1, int n1, int n2, double* M)
{
  for (int r=0; r<m1; r++) {
    for (int c=0; c<n2; c++) {
      double sum = 0.0;
      for (int r2=0; r2<n1; r2++) {
        sum += M1[r*n1+r2]*M2[r2*n1+c];
      }
      M[r*n1+c] = sum;
    }
  }
}

void copy(double* M1, double* M2, int m, int n)
{
  for (int r=0; r<m; r++) {
    for (int c=0; c<n; c++) {
      M1[r*n+c] = M2[r*n+c];
    }   
  }
}

void sub_mul_v(double* M1, double* M2, int m, int n, double v, double* M)
{
  for (int r=0; r<m; r++) {
    for (int c=0; c<n; c++) {
      M[r*n+c] = M1[r*n+c] - v*M2[r*n+c];
    }
  }
}

void mult_v(double* M, int m, int n, double v)
{
  for (int r=0; r<m; r++) {
    for (int c=0; c<n; c++) {
      M[r*n+c] *= v;
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
  static double x4[100];

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
        x4[i] = 0;
      return x4; 
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

void col_sub(double* M1, double* M2, int m, int n, int col1, int col2)
{
  for (int r=0; r<m; r++) {
    M1[r*n + col1] -= M2[r*n + col2];
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

void col_cross_mult(double* M1, double* M2, int m, int n, int col1, int col2, double* M)
{
  for (int r=0; r<m; r++) {
    for (int c=0; c<m; c++) {
      M[r*m + c] = M1[r*n + col1] * M2[c*n + col2];
    }
  }
}

void matrix_vec_mul(double* M, int m, int n, double* v, double* v_out)
{
  for (int r=0; r<m; r++) {
    double sum = 0.0;
    for (int c=0; c<n; c++) {
      sum += M[r*n+c]*v[c];
    }
    v_out[r] = sum;
  }
}

// modified gram-schmidt
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

double* I(int n)
{
  static double buffer[100] = {0};

  for (int j=0; j<100; j++) {
    buffer[j] = 0;
  }

  for (int i=0; i<n; i++) {
    buffer[i*n + i] = 1;
  }

  return buffer;
}
/*
void householder(double* A, int m, int n, double* Q, double *R)
{
  double H[50] = {0};
  double V[50] = {0};

  for(int j=0; j<n; j++) {
    double r_ii = col_norm(A, m, n, c);
 
    // w = (r_ii, 0, ..., 0)
    V[j*n + j] = r_ii;

    //( w - Aj )/2
    col_sub(V, A, m, n, 0, j);
  
    //H = I - 2 *(v * vt)/(vt * v)
    double vtv = col_norm(V, m, n, c);
    col_cross_mult(V, V, m, n, c, c, H);
    sub_mul_v(I(m), m, H, 2.0/vtv, H);

    //
    copy(Q);
  }
}
*/

double* matrix() {
  static double a[3][2];

  a[0][0] = 1; 
  a[0][1] = -4;
  a[1][0] = 2; 
  a[1][1] = 3;
  a[2][0] = 2; 
  a[2][1] = 2;

  return (double *) a;
}

double* b() {
  static double v[3];

  v[0] = -3;
  v[1] = 15;
  v[2] = 9;

  return v;
}

void Uxb(double* M, int n, double* b, double* x) {
  for (int i=n-1; i>=0; i--) {
    double c = M[i*n+i];
    b[i] /= c;
    for (int j=i+1; j<n; j++) {
      b[i] -= M[i*n+j]/c*x[j];
    }
    x[i] = b[i];
  }
}

int main()
{
  double *Q = buffer(1);
  double *R = buffer(2);
  double *Qt = buffer(3);
  double *b_out = buffer(4);

  printf("\n----- The A matrix -----\n");
  print_matrix(matrix(), 3, 2);

  qr(matrix(), 3, 2, Q, R);

  printf("\n----- The Q matrix -----\n");
  print_matrix(Q, 3, 2);

  printf("\n----- The R matrix -----\n");
  print_matrix(R, 2, 2);

  printf("\n----- The Qt matrix -----\n");
  trans(Qt, Q, 3, 2);
  print_matrix(Qt, 2, 3);

  printf("\n----- The Qt * b --------\n");
  matrix_vec_mul(Qt, 2, 3, b(), b_out);
  print_vec(b_out, 2);

  printf("\n--- The resolutions of Rx = Qt * b ---\n");
  double x[2] = { 0 };
  Uxb(R, 2, b_out, x);
  print_vec(x, 2);

  return 0;
}
