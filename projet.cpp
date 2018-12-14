#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>

void data_export(const char *_file_name,
                 const int _Nx,
                 const int _Ny,
                 const double *_x,
                 const double *_y,
                 const double **_f);

void upwind(const int Nx, const int Ny, double a_x, double a_y, const double CFL);

void alternate_upwind(const int Nx, const int Ny, double a_x, double a_y, const double CFL);

inline double circle(const double _x, const double _y) {
  if (_x*_x + _y*_y < 0.01) {
    return 1;
  } else {
    return 0;
  }
}

int main() {
	upwind(100, 100, 1, 1, 0.5);
  alternate_upwind(100, 100, 1, 1, 0.5);
	return 0;
}

void upwind(const int Nx, const int Ny, double a_x, double a_y, const double CFL) {

  const double x_min = -1.;
  const double x_max = 1.;
  const double y_min = -1.;
  const double y_max = 1.;

  double *x;
  x = new double[Nx + 2];

  for (int i = 0; i < Nx + 2; i++) {
    x[i] = x_min + (x_max - x_min)/double(Nx)*(double(i) + 0.5); 
  }

  double *y;
  y = new double[Ny + 2];

  for (int j = 0; j < Ny + 2; j++) {
    y[j] = y_min + (y_max - y_min)/double(Ny)*(double(j) + 0.5);
  }

  const double t_max = 0.1;
  const double step_x = (x_max - x_min)/double(Nx);
  const double step_y = (y_max - y_min)/double(Ny);

  const double step = 1/(fabs(a_x)/step_x + fabs(a_y)/step_y) * CFL;
  const int N = int(t_max/step);

  double **u0;
  double **u1;
  double **u_star_x;
  double **u_star_y;
  int test_x;
  int test_y;

  u0 = new double *[Nx + 2];
  u1 = new double *[Nx + 2];

  u_star_x = new double *[Nx + 1];
  u_star_y = new double *[Nx + 2];

  for (int i = 0; i < Nx + 2; i++) {
    u_star_y[i] = new double[Ny + 1]; 
  } 

  for (int fi = 0; fi < Nx + 1; fi++) {
    u_star_x[fi] = new double[Ny + 2];
  }

  for (int i = 0; i < Nx + 2; i++) {
    u0[i] = new double[Ny + 2];
    u1[i] = new double[Ny + 2];
    for (int j = 0; j < Ny + 2; j++) {
      u0[i][j] = circle(x[i], y[j]);
      u1[i][j] = 0.0;
    }
  }

    if (a_x > 0) {
      test_x = 1;
    } else {
      test_x = 0;
    }

    if (a_y > 0) {
      test_y = 1;
    } else {
      test_y = 0;
    }

  data_export("upwind0", Nx, Ny, x, y, (const double **)u0);
  for (int t = 0; t < N; t++) {
    for (int fi = 0; fi < Nx + 1; fi++) {
      for (int j = 0; j < Ny + 2; j++) {
        u_star_x[fi][j] = (1 - test_x)*u0[fi + 1][j] + test_x*u0[fi][j];
      }
    }

    for (int i = 0; i < Nx + 2; i++) {
      for (int fj = 0; fj < Ny + 1; fj++) {
        u_star_y[i][fj] = (1 - test_y)*u0[i][fj + 1] + test_y*u0[i][fj];
      }
    }

  
    for (int i = 1; i < Nx + 1; i++) {
      for (int j = 1; j < Ny + 1; j++) {
        u1[i][j] = u0[i][j] - (u_star_x[i][j] - u_star_x[i - 1][j])*(step/step_x) - (u_star_y[i][j] - u_star_y[i][j - 1])*(step/step_y);
        u0[i][j] = u1[i][j];
      }
    }
  }
  data_export("upwind1", Nx, Ny, x, y, (const double **)u0);

  for (int nx  = 0; nx < Nx + 2; nx++) {
  delete [] u0[nx];
  delete [] u1[nx];
  }
  delete [] u0;
  delete [] u1;
  for (int fi = 0; fi < Nx + 1; fi++) {
    delete [] u_star_x[fi];
  }
  delete [] u_star_x;
  for (int fi = 0; fi < Nx + 2; fi++) {
    delete [] u_star_y[fi];
  }
  delete [] u_star_y;

  delete [] x;
  delete [] y;
}

void alternate_upwind(const int Nx, const int Ny, double a_x, double a_y, const double CFL) {

  const double x_min = -1.;
  const double x_max = 1.;
  const double y_min = -1.;
  const double y_max = 1.;

  double *x;
  x = new double[Nx + 2];

  for (int i = 0; i < Nx + 2; i++) {
    x[i] = x_min + (x_max - x_min)/double(Nx)*(double(i) + 0.5); 
  }

  double *y;
  y = new double[Ny + 2];

  for (int j = 0; j < Ny + 2; j++) {
    y[j] = y_min + (y_max - y_min)/double(Ny)*(double(j) + 0.5);
  }

  const double t_max = 0.1;
  const double step_x = (x_max - x_min)/double(Nx);
  const double step_y = (y_max - y_min)/double(Ny);

  const double step = 0.001;
  const int N = int(t_max/step);

  double **u0;
  double **u1_tilde;
  double **u1;

  double **u_star_x;
  double **u_star_tilde_y;

  int test_x;
  int test_y;

  u0 = new double *[Nx + 2];
  u1_tilde = new double *[Nx + 2];
  u1 = new double *[Nx + 2];

  u_star_x = new double *[Nx + 1];
  u_star_tilde_y = new double *[Nx + 2];

  for (int i = 0; i < Nx + 2; i++) {
    u_star_tilde_y[i] = new double[Ny + 1];
  } 

  for (int fi = 0; fi < Nx + 1; fi++) {
    u_star_x[fi] = new double[Ny + 2];
  }

  for (int i = 0; i < Nx + 2; i++) {
    u0[i] = new double[Ny + 2];
    u1[i] = new double[Ny + 2];
    u1_tilde[i] = new double[Ny + 2];

    for (int j = 0; j < Ny + 2; j++) {
      u0[i][j] = circle(x[i], y[j]);
      u1[i][j] = 0.0;
      u1_tilde[i][j] = 0.0;
    }
  }

  if (a_x > 0) {
      test_x = 1;
    } else {
      test_x = 0;
    }

    if (a_y > 0) {
      test_y = 1;
    } else {
      test_y = 0;
    }

  data_export("alternate_upwind0", Nx, Ny, x, y, (const double **)u0);
  for (int t = 0; t < N; t++) {
    for (int fi = 0; fi < Nx + 1; fi++) {
      for (int j = 0; j < Ny + 2; j++) {
        u_star_x[fi][j] = (1 - test_x)*u0[fi + 1][j] + test_x*u0[fi][j];
      }
    }

    for (int i = 1; i < Nx + 1; i++) {
      for (int j = 1; j < Ny + 1; j++) {
        u1_tilde[i][j] = u0[i][j] - (u_star_x[i][j] - u_star_x[i - 1][j])*(step/step_x);
      }
    }

    for (int i = 0; i < Nx + 2; i++) {
      for (int fj = 0; fj < Ny + 1; fj++) {
        u_star_tilde_y[i][fj] = (1 - test_y)*u1_tilde[i][fj + 1] + test_y*u1_tilde[i][fj];
      }
    }


    for (int i = 1; i < Nx + 1; i++) {
      for (int j = 1; j < Ny + 1; j++) {
        u1[i][j] = u1_tilde[i][j] - (u_star_tilde_y[i][j] - u_star_tilde_y[i][j - 1])*(step/step_y);
        u0[i][j] = u1[i][j];
      }
    }
  }
  data_export("alternate_upwind1", Nx, Ny, x, y, (const double **)u0);

  for (int nx  = 0; nx < Nx + 2; nx++) {
  delete [] u0[nx];
  delete [] u1_tilde[nx];
  delete [] u1[nx];
  }
  delete [] u0;
  delete [] u1_tilde;
  delete [] u1;

  for (int fi = 0; fi < Nx + 1; fi++) {
    delete [] u_star_x[fi];
  }
  delete [] u_star_x;

  for (int i = 0; i < Nx + 2; i++) {
    delete [] u_star_tilde_y[i];
  }
  delete [] u_star_tilde_y;

  delete [] x;
  delete [] y;
}

void data_export(const char *_file_name,
                 const int _Nx,
                 const int _Ny,
                 const double *_x,
                 const double *_y,
                 const double **_f) {
  
  FILE *file_pointer;
  
  file_pointer = fopen((const char*) _file_name, "w+");

  fprintf(file_pointer, "# x \t y \t f(x,y) \n");

  for (int nx = 1; nx < _Nx + 1; nx++) {
    for (int ny = 1; ny < _Ny + 1; ny++) {
      fprintf(file_pointer, "%e \t %e \t %e \n", _x[nx], _y[ny], _f[nx][ny]);
    }
    fprintf(file_pointer, "\n");
  }
  fprintf(file_pointer, "\n");

  fclose(file_pointer);
}