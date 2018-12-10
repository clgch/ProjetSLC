#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

void data_export(const char *_file_name,
								 const int _Nx,
								 const int _Ny,
								 const double *_x,
								 const double *_y,
								 const double **_f,
								 int t);

inline double fct_f( const double _x, const double _y) {
	const double pi = M_PI;
	return sin(2. * pi * _x ) * exp(-_y*_y);
}

int main() {

const char* nom_de_fichier = "out.txt";

const double x_min = 0.;
const double x_max = 1.;
const double y_min = 0.;
const double y_max = 1.;

const int Nx = 50;
const int Ny = 50;

//definition du maillage 

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

// definition d'un champ scalaire

double **f;

f = new double *[Nx + 2];

for (int i = 0; i < Nx + 2; i++) {
	f[i] = new double[Ny + 2];
	for (int j = 0; j < Ny + 2; j++) {
		f[i][j] = fct_f(x[i], y[j]);
	}
}

// schema upwind, ax = 1 et ay = 1

const int N = 100;
const int t_max = 1;
const double step = double(t_max/N);
const double step_x = (x_max - x_min)/double(Nx);
const double step_y = (y_max - y_min)/double(Ny);

double **u0;
double **u1;

u0 = new double *[Nx + 2];
u1 = new double *[Nx + 2];

for (int i = 0; i < Nx + 2; i++) {
	u0[i] = new double[Ny + 2];
	u1[i] = new double[Ny + 2];
	for (int j = 0; j < Ny + 2; j++) {
		u0[i][j] = fct_f(x[i], y[j]);
		u1[i][j] = 0.0;
	}
}

for (int t = 0; t < Nx; t++) {
	data_export(nom_de_fichier, Nx, Ny, x, y, (const double **)u0, t);
	for (int i = 1; i < Nx + 1; i++) {
		for (int j = 1; j < Ny + 1; j++) {
			u1[i][j] = u0[i][j] - (u0[i][j] - u0[i - 1][j])*(step/step_x) - (u0[i][j] - u0[i][j-1])*(step/step_y);
			u0[i][j] = u1[i][j];
		}
	}
}


//data_export(nom_de_fichier, Nx, Ny, x, y, (const double **)f);


//desalloc

delete [] x;
delete [] y;
for (int nx = 0; nx < Nx + 2; nx++) {
	delete [] f[nx];
}
delete [] f;
for(int nx  = 0; nx < Nx + 2; nx++) {
	delete [] u0[nx];
	delete [] u1[nx];
}
delete [] u0;
delete [] u1;

return 0;
}

// fonction d'exportation du tableau

void data_export(const char *_file_name,
								 const int _Nx,
								 const int _Ny,
								 const double *_x,
								 const double *_y,
								 const double **_f,
								 int t) {
	
	FILE *file_pointer;

	file_pointer = fopen((const char*)(_file_name), "a");

	if (t == 0) {
			fprintf(file_pointer, "# t \t x \t y \t f(x,y) \n");
	}

	for (int nx = 1; nx < _Nx + 1; nx++) {
		for (int ny = 1; ny < _Ny + 1; ny++) {
			fprintf(file_pointer, "%i \t %e \t %e \t %e \n",t , _x[nx], _y[ny], _f[nx][ny]);
		}
		fprintf(file_pointer, "\n");
	}
	fprintf(file_pointer, "\n");

	fclose(file_pointer);
	
}