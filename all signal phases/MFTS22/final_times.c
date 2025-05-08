#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdbool.h>

#define Power(a, b) pow(a, b)

double XMAX = 300.0;
double XMIN = 0.0;
double VMAX = 16.0;
double VMIN = 0.0;

double x_space = 300.0 + 50.0;

int numsteps = 60;
double T = 1.0;

/* analytic parameters  */
double ve = 16.0;
double xe = 300.0 + 70.0;
double w = 0.1;
double tf = 38.0;
double wopt = 1.0;

double u_min = -3.0, u_max = 3.0;

/* Newton parameters */
int newtonIntervals = 1;
double newton_add_term = 5.0;
int negSpeedInitCon;

double du = 0.125;

double ***J_opt, ***te_opt, ***J_opt_w;

struct {
	int NX, *NXlive_min, *NXlive_max;
	double dx, **X;
	double *xmin, *xmax;

	int NV, *NVlive_min, *NVlive_max;
	double *vmin, *vmax, dv, **V;

	int NU;
	double umin, umax, du;
}D = { 0 };

void initialize() {
    int i, k, ix, iv;

    D.NV = (int)ceil((VMAX - VMIN) / du) + 1;
    D.NX = (int)ceil((x_space - XMIN) / (du/2.0)) + 1;
    D.NU = (int)ceil((u_max - u_min) / du) + 1;

    fprintf(stderr, "HELLO %d \n", D.NX);

    assert((D.V = (double**)calloc(sizeof(double*), D.NV)));
	for (i = 0; i < D.NV; i++)
		assert((D.V[i] = (double*)calloc(sizeof(double), numsteps)));
	assert((D.X = (double**)calloc(sizeof(double*), D.NX)));
	for (i = 0; i < D.NX; i++)
		assert((D.X[i] = (double*)calloc(sizeof(double), numsteps)));

    D.dv = du;
    for (k = 0; k < numsteps; k++){
		for (iv = 0; iv < D.NV; iv++) {
			D.V[iv][k] = VMIN + D.dv * iv;
            // fprintf(stderr, "D.V[%d][%d]: %.4f \n", iv, k, D.V[iv][k]);
		}
	}

    D.dx = 0.5 * du;
    for (k = 0; k < numsteps; k++) {
		for (ix = 0; ix < D.NX; ix++) {
			D.X[ix][k] = XMIN + D.dx * ix;
            // fprintf(stderr, "D.X[%d][%d]: %.4f \n", ix, k, D.X[ix][k]);
		}
	}


    assert((D.NXlive_min = (int*)calloc(sizeof(int), numsteps)));
	assert((D.NXlive_max = (int*)calloc(sizeof(int), numsteps)));
	assert((D.NVlive_min = (int*)calloc(sizeof(int), numsteps)));
	assert((D.NVlive_max = (int*)calloc(sizeof(int), numsteps)));
	
	for (k = 0; k < numsteps; k++) {
		D.NXlive_min[k] = 0;
		D.NXlive_max[k] = D.NX;

		D.NVlive_min[k] = 0;
		D.NVlive_max[k] = D.NV;
	}

    
	assert((J_opt = (double ***)calloc(sizeof(double*), numsteps)));
    assert((J_opt_w = (double ***)calloc(sizeof(double*), numsteps)));
	assert((te_opt = (double ***)calloc(sizeof(double*), numsteps)));
	for (k = 0; k < numsteps; k++) {
		assert((J_opt[k] = (double **)calloc(sizeof(double*), D.NX)));
        assert((J_opt_w[k] = (double **)calloc(sizeof(double*), D.NX)));
		assert((te_opt[k] = (double **)calloc(sizeof(double*), D.NX)));
		for (ix = 0; ix < D.NX; ix++) {
			assert((J_opt[k][ix] = (double *)calloc(sizeof(double), D.NV)));
            assert((J_opt_w[k][ix] = (double *)calloc(sizeof(double), D.NV)));
			assert((te_opt[k][ix] = (double *)calloc(sizeof(double), D.NV)));
		}
	}

}

double UP_control(double t, double t0, double x0, double v0, double te) {
		
	double res = (5.999999999999999*t*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)));

	return res;
}

double UP_position(double t, double t0, double x0, double v0, double te) {

	return (0.9999999999999998*Power(t,3)*(1.*Power(t0,2)*v0 - 2.*t0*te*v0 + 1.*Power(te,2)*v0 + 
        1.*Power(t0,2)*ve - 2.*t0*te*ve + 1.*Power(te,2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe
        ))/(1.*Power(t0,4) - 3.9999999999999996*Power(t0,3)*te + 
      5.999999999999999*Power(t0,2)*Power(te,2) - 3.9999999999999996*t0*Power(te,3) + 1.*Power(te,4))\
    - (0.9999999999999999*Power(t,2)*(1.*Power(t0,4)*v0 - 1.*Power(t0,3)*te*v0 - 
        3.*Power(t0,2)*Power(te,2)*v0 + 5.*t0*Power(te,3)*v0 - 2.*Power(te,4)*v0 + 2.*Power(t0,4)*ve - 
        5.*Power(t0,3)*te*ve + 3.*Power(t0,2)*Power(te,2)*ve + 1.*t0*Power(te,3)*ve - 
        1.*Power(te,4)*ve - 3.*Power(t0,3)*x0 + 3.*Power(t0,2)*te*x0 + 3.*t0*Power(te,2)*x0 - 
        3.*Power(te,3)*x0 + 3.*Power(t0,3)*xe - 3.*Power(t0,2)*te*xe - 3.*t0*Power(te,2)*xe + 
        3.*Power(te,3)*xe))/
    ((1.*t0 - 1.*te)*(1.*Power(t0,4) - 3.9999999999999996*Power(t0,3)*te + 
        5.999999999999999*Power(t0,2)*Power(te,2) - 3.9999999999999996*t0*Power(te,3) + 1.*Power(te,4))
      ) + (1.*t*(1.9999999999999998*Power(t0,4)*te*v0 - 4.999999999999999*Power(t0,3)*Power(te,2)*v0 + 
        2.9999999999999996*Power(t0,2)*Power(te,3)*v0 + 1.*t0*Power(te,4)*v0 - 1.*Power(te,5)*v0 + 
        1.*Power(t0,5)*ve - 1.*Power(t0,4)*te*ve - 2.9999999999999996*Power(t0,3)*Power(te,2)*ve + 
        4.999999999999999*Power(t0,2)*Power(te,3)*ve - 1.9999999999999998*t0*Power(te,4)*ve - 
        5.999999999999999*Power(t0,3)*te*x0 + 11.999999999999998*Power(t0,2)*Power(te,2)*x0 - 
        5.999999999999999*t0*Power(te,3)*x0 + 5.999999999999999*Power(t0,3)*te*xe - 
        11.999999999999998*Power(t0,2)*Power(te,2)*xe + 5.999999999999999*t0*Power(te,3)*xe))/
    ((1.*t0 - 1.*te)*(1.*Power(t0,4) - 3.9999999999999996*Power(t0,3)*te + 
        5.999999999999999*Power(t0,2)*Power(te,2) - 3.9999999999999996*t0*Power(te,3) + 1.*Power(te,4))
      ) - (1.*(0.9999999999999999*Power(t0,4)*Power(te,2)*v0 - 
        2.9999999999999996*Power(t0,3)*Power(te,3)*v0 + 3.*Power(t0,2)*Power(te,4)*v0 - 
        1.*t0*Power(te,5)*v0 + 1.*Power(t0,5)*te*ve - 3.*Power(t0,4)*Power(te,2)*ve + 
        2.9999999999999996*Power(t0,3)*Power(te,3)*ve - 
        0.9999999999999999*Power(t0,2)*Power(te,4)*ve - 
        2.9999999999999996*Power(t0,3)*Power(te,2)*x0 + 6.999999999999999*Power(t0,2)*Power(te,3)*x0 - 
        5.*t0*Power(te,4)*x0 + 1.*Power(te,5)*x0 - 1.*Power(t0,5)*xe + 5.*Power(t0,4)*te*xe - 
        6.999999999999999*Power(t0,3)*Power(te,2)*xe + 2.9999999999999996*Power(t0,2)*Power(te,3)*xe))/
    ((1.*t0 - 1.*te)*(1.*Power(t0,4) - 3.9999999999999996*Power(t0,3)*te + 
        5.999999999999999*Power(t0,2)*Power(te,2) - 3.9999999999999996*t0*Power(te,3) + 1.*Power(te,4))
      );
}

double UP_speed(double t, double t0, double x0, double v0, double te) {
	return (2.9999999999999996*Power(t,2)*(1.*Power(t0,2)*v0 - 2.*t0*te*v0 + 1.*Power(te,2)*v0 + 
        1.*Power(t0,2)*ve - 2.*t0*te*ve + 1.*Power(te,2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe
        ))/(1.*Power(t0,4) - 3.9999999999999996*Power(t0,3)*te + 
      5.999999999999999*Power(t0,2)*Power(te,2) - 3.9999999999999996*t0*Power(te,3) + 1.*Power(te,4))\
    - (1.9999999999999998*t*(1.*Power(t0,4)*v0 - 1.*Power(t0,3)*te*v0 - 
        3.*Power(t0,2)*Power(te,2)*v0 + 5.*t0*Power(te,3)*v0 - 2.*Power(te,4)*v0 + 2.*Power(t0,4)*ve - 
        5.*Power(t0,3)*te*ve + 3.*Power(t0,2)*Power(te,2)*ve + 1.*t0*Power(te,3)*ve - 
        1.*Power(te,4)*ve - 3.*Power(t0,3)*x0 + 3.*Power(t0,2)*te*x0 + 3.*t0*Power(te,2)*x0 - 
        3.*Power(te,3)*x0 + 3.*Power(t0,3)*xe - 3.*Power(t0,2)*te*xe - 3.*t0*Power(te,2)*xe + 
        3.*Power(te,3)*xe))/
    ((1.*t0 - 1.*te)*(1.*Power(t0,4) - 3.9999999999999996*Power(t0,3)*te + 
        5.999999999999999*Power(t0,2)*Power(te,2) - 3.9999999999999996*t0*Power(te,3) + 1.*Power(te,4))
      ) + (1.*(1.9999999999999998*Power(t0,4)*te*v0 - 4.999999999999999*Power(t0,3)*Power(te,2)*v0 + 
        2.9999999999999996*Power(t0,2)*Power(te,3)*v0 + 1.*t0*Power(te,4)*v0 - 1.*Power(te,5)*v0 + 
        1.*Power(t0,5)*ve - 1.*Power(t0,4)*te*ve - 2.9999999999999996*Power(t0,3)*Power(te,2)*ve + 
        4.999999999999999*Power(t0,2)*Power(te,3)*ve - 1.9999999999999998*t0*Power(te,4)*ve - 
        5.999999999999999*Power(t0,3)*te*x0 + 11.999999999999998*Power(t0,2)*Power(te,2)*x0 - 
        5.999999999999999*t0*Power(te,3)*x0 + 5.999999999999999*Power(t0,3)*te*xe - 
        11.999999999999998*Power(t0,2)*Power(te,2)*xe + 5.999999999999999*t0*Power(te,3)*xe))/
    ((1.*t0 - 1.*te)*(1.*Power(t0,4) - 3.9999999999999996*Power(t0,3)*te + 
        5.999999999999999*Power(t0,2)*Power(te,2) - 3.9999999999999996*t0*Power(te,3) + 1.*Power(te,4))
      );
}

double te_func(double t0, double x0, double v0, double te)
{
	double res = 0.5*w + 0.5*Power((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 +
		2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 +
			5.000000000000001*t0*Power(te, 3)*v0 - 2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve +
			3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve + 0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve -
			3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 + 3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 +
			3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe - 3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))), 2) +
		((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe -
			2.*te*xe)) / (1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
			(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 +
				5.000000000000001*t0*Power(te, 3)*v0 - 2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve +
				3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve + 0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve -
				3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 + 3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 +
				3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe - 3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))*
		((-5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 +
			2.*t0*xe - 2.*te*xe)) / (1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
			(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 +
				5.000000000000001*t0*Power(te, 3)*v0 - 2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve +
				3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve + 0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve -
				3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 + 3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 +
				3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe - 3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)))) +
		(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe -
			2.*te*xe)*((2.9999999999999996*Power(te, 2)*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve -
				2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
				(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
				(1.9999999999999996*te*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 +
					5.000000000000001*t0*Power(te, 3)*v0 - 2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve +
					3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve + 0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve -
					3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 + 3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 +
					3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe - 3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
					((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
				(1.*(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 +
					1.0000000000000004*t0*Power(te, 4)*v0 - 1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve -
					2.999999999999999*Power(t0, 3)*Power(te, 2)*ve + 4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve -
					5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 - 5.999999999999999*t0*Power(te, 3)*x0 +
					5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
					((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))) /
				(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4));

	return res;
}

double dte_func(double t0, double x0, double v0, double te)
{
	
	double res = 1.*((5.999999999999999*te*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(5.999999999999999*te*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2) +
		(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*(-0.9999999999999998*Power(t0, 3)*v0 - 6.000000000000001*Power(t0, 2)*te*v0 + 15.000000000000004*t0*Power(te, 2)*v0 - 8.000000000000002*Power(te, 3)*v0 -
			5.000000000000001*Power(t0, 3)*ve + 6.000000000000001*Power(t0, 2)*te*ve + 2.999999999999999*t0*Power(te, 2)*ve - 3.9999999999999996*Power(te, 3)*ve +
			3.0000000000000004*Power(t0, 2)*x0 + 6.000000000000001*t0*te*x0 - 9.000000000000002*Power(te, 2)*x0 - 3.0000000000000004*Power(t0, 2)*xe - 6.000000000000001*t0*te*xe +
			9.000000000000002*Power(te, 2)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.9999999999999996*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))*
		((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
			(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
				2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
				0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
				3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
				3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)))) +
			((-5.999999999999999*te*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
				(5.999999999999999*te*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
				(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
				Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2) -
				(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
				(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
				(1.9999999999999996*(-0.9999999999999998*Power(t0, 3)*v0 - 6.000000000000001*Power(t0, 2)*te*v0 + 15.000000000000004*t0*Power(te, 2)*v0 - 8.000000000000002*Power(te, 3)*v0 -
					5.000000000000001*Power(t0, 3)*ve + 6.000000000000001*Power(t0, 2)*te*ve + 2.999999999999999*t0*Power(te, 2)*ve - 3.9999999999999996*Power(te, 3)*ve +
					3.0000000000000004*Power(t0, 2)*x0 + 6.000000000000001*t0*te*x0 - 9.000000000000002*Power(te, 2)*x0 - 3.0000000000000004*Power(t0, 2)*xe - 6.000000000000001*t0*te*xe +
					9.000000000000002*Power(te, 2)*xe)) /
					((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) -
				(1.9999999999999996*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
				(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
					2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
					0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
					3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
					3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
					((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) +
				(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
					2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
					0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
					3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
					3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
					(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))*
				((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
					(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
						2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
						0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
						3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
						3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
						((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)))) +
					((5.999999999999999*te*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
						(5.999999999999999*te*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
						(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
						Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2) +
						(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
						(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
						(1.9999999999999996*(-0.9999999999999998*Power(t0, 3)*v0 - 6.000000000000001*Power(t0, 2)*te*v0 + 15.000000000000004*t0*Power(te, 2)*v0 - 8.000000000000002*Power(te, 3)*v0 -
							5.000000000000001*Power(t0, 3)*ve + 6.000000000000001*Power(t0, 2)*te*ve + 2.999999999999999*t0*Power(te, 2)*ve - 3.9999999999999996*Power(te, 3)*ve +
							3.0000000000000004*Power(t0, 2)*x0 + 6.000000000000001*t0*te*x0 - 9.000000000000002*Power(te, 2)*x0 - 3.0000000000000004*Power(t0, 2)*xe - 6.000000000000001*t0*te*xe +
							9.000000000000002*Power(te, 2)*xe)) /
							((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.9999999999999996*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
		((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
		(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))*
		((-5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)))) +
			(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)*
		((2.9999999999999996*Power(te, 2)*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
				(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(2.9999999999999996*Power(te, 2)*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2) +
		(5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*te*(-0.9999999999999998*Power(t0, 3)*v0 - 6.000000000000001*Power(t0, 2)*te*v0 + 15.000000000000004*t0*Power(te, 2)*v0 - 8.000000000000002*Power(te, 3)*v0 -
			5.000000000000001*Power(t0, 3)*ve + 6.000000000000001*Power(t0, 2)*te*ve + 2.999999999999999*t0*Power(te, 2)*ve - 3.9999999999999996*Power(te, 3)*ve +
			3.0000000000000004*Power(t0, 2)*x0 + 6.000000000000001*t0*te*x0 - 9.000000000000002*Power(te, 2)*x0 - 3.0000000000000004*Power(t0, 2)*xe - 6.000000000000001*t0*te*xe +
			9.000000000000002*Power(te, 2)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.*(1.9999999999999996*Power(t0, 4)*v0 - 9.999999999999996*Power(t0, 3)*te*v0 + 8.999999999999996*Power(t0, 2)*Power(te, 2)*v0 + 4.000000000000002*t0*Power(te, 3)*v0 -
			5.000000000000001*Power(te, 4)*v0 - 1.0000000000000004*Power(t0, 4)*ve - 5.999999999999998*Power(t0, 3)*te*ve + 14.999999999999995*Power(t0, 2)*Power(te, 2)*ve -
			7.999999999999997*t0*Power(te, 3)*ve - 5.999999999999999*Power(t0, 3)*x0 + 23.999999999999996*Power(t0, 2)*te*x0 - 17.999999999999996*t0*Power(te, 2)*x0 +
			5.999999999999999*Power(t0, 3)*xe - 23.999999999999996*Power(t0, 2)*te*xe + 17.999999999999996*t0*Power(te, 2)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.9999999999999996*te*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) -
		(1.9999999999999996*te*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) -
		(1.*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 + 1.0000000000000004*t0*Power(te, 4)*v0 -
			1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve - 2.999999999999999*Power(t0, 3)*Power(te, 2)*ve +
			4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve - 5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 -
			5.999999999999999*t0*Power(te, 3)*x0 + 5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) +
		(1.*(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 + 1.0000000000000004*t0*Power(te, 4)*v0 -
			1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve - 2.999999999999999*Power(t0, 3)*Power(te, 2)*ve +
			4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve - 5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 -
			5.999999999999999*t0*Power(te, 3)*x0 + 5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
			(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
		(5.999999999999999*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)*
		((2.9999999999999996*Power(te, 2)*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe -
		2.*te*xe)) / (1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*te*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.*(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 + 1.0000000000000004*t0*Power(te, 4)*v0 -
			1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve - 2.999999999999999*Power(t0, 3)*Power(te, 2)*ve +
			4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve - 5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 -
			5.999999999999999*t0*Power(te, 3)*x0 + 5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(5.999999999999999*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)*
		((2.9999999999999996*Power(te, 2)*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe -
			2.*te*xe)) / (1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
			(1.9999999999999996*te*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
				2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
				0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
				3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
				3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
			(1.*(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 + 1.0000000000000004*t0*Power(te, 4)*v0 -
				1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve - 2.999999999999999*Power(t0, 3)*Power(te, 2)*ve +
				4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve - 5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 -
				5.999999999999999*t0*Power(te, 3)*x0 + 5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))) /
		Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2);

			return res;
}

double newton_raphson(double t0, double x0, double v0, double te) {
	
	/* Newton - Raphson Method */
	int iter = 0;
	double h = te_func(t0, x0, v0, te) / dte_func(t0, x0, v0, te);
	while (fabs(h) >= 0.001 && iter < 100)
	{
		h = te_func(t0, x0, v0, te) / dte_func(t0, x0, v0, te);
		  
		te = te - h;
		iter += 1;
	}

	return te;
}

void optimal_times_cost(double x0, double v0, int t0, int ix, int iv) {
    double te;
    
    double te_init, temp_te, check_te, temp_cost = DBL_MAX;
    double checkSpeed = v0;
    int negativeSpeed = 0;
    
    /* calculate max te assuming the min speed is 1 */
    double max_te = (xe - x0) / 1.0;
    te_init = t0 + 1;
    for (int i = 0; i < newtonIntervals; i++) {
        double te_cost = 0.0;

        check_te = newton_raphson(t0, x0, v0, te_init);
        if (check_te > 2.0*max_te || isnan(check_te) || check_te < 0.0)
            continue;
        
        for (int i = 0; i < numsteps; i++) {
            double t = (i == 0) ? 0.0:(i*(check_te - t0) / (numsteps - 1) + t0);
            te_cost += 0.5*pow(UP_control(t, t0, x0, v0, check_te), 2);
            /* skip optimal te that leads to negative speed trajectories */
            checkSpeed = UP_speed(t, t0, x0, v0, check_te);
            if (checkSpeed < 0.0) {
                negativeSpeed = 1;
            }
        }

        if ((check_te > t0) && (te_cost < temp_cost) && (negativeSpeed != 1)) {
            temp_cost = te_cost;
            temp_te = check_te;
        }
        te_init += newton_add_term;
    }
    te = temp_te;
    te_opt[t0][ix][iv] = te;

    // calculate the optimal control cost
    double step;
    double result = 0.0;
    for (int i = 0; i < 50000; i++) {
        step = (te - t0) / (50000 - 1);
        double t = (i*step + t0);
        // fprintf(stderr, "%.4f  ", t);
        result += 0.5* (pow(UP_control(t, t0, x0, v0, te), 2));
    }
    result *= step;
    J_opt[t0][ix][iv] = wopt * result;

    result += 0.5 * w * te;
    J_opt_w[t0][ix][iv] = wopt * result;
    // fprintf(stderr, "\n J_opt[%d][%d][%d]: %.4f - J_opt_w: %.4f -- te: %.4f -- x: %.4f -- v: %.4f \n", t0, ix, iv, J_opt[t0][ix][iv], J_opt_w[t0][ix][iv], te_opt[t0][ix][iv], x0, v0);
}

int main() {
    int k, ix, iv;
    double x, v;
    initialize();

    FILE *f1, *f2, *f3;
    f1 = fopen("j_opt_te=16.0.txt", "w");
    f2 = fopen("j_opt_w_te=16.0.txt", "w");
    f3 = fopen("te_opt_te=16.0.txt", "w");

    fputs("var Data = {\n", f1);
    fputs("var Data = {\n", f2);
    fputs("var Data = {\n", f3);
    // for (k = 0; k < numsteps; k++) {
    for (k = 0; k < 1; k++) {
		for (ix = D.NXlive_min[k]; ix < D.NXlive_max[k]; ix++) {
        // for (ix = 0; ix < 2; ix++) {
			x = D.X[ix][k];
			for (iv = D.NVlive_min[k]; iv < D.NVlive_max[k]; iv++) {
            // for (iv = 0; iv < 2; iv++) {
				v = D.V[iv][k];

                optimal_times_cost(x, v, k, ix, iv);

                fprintf(f1, "\"J_opt(%d,%d,%d)\":%.4f,\n", k, ix, iv, J_opt[k][ix][iv]);
                fprintf(f2, "\"J_opt_w(%d,%d,%d)\":%.4f,\n", k, ix, iv, J_opt_w[k][ix][iv]);
                fprintf(f3, "\"te_opt(%d,%d,%d)\":%.4f,\n", k, ix, iv, te_opt[k][ix][iv]);
            }
        }
    }
    fprintf(f1, "};\n\n");
    fprintf(f2, "};\n\n");
    fprintf(f3, "};\n\n");

    fclose(f1);
    fclose(f2);
    fclose(f3);

    return 0;
}