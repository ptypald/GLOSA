#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// #include "opt_costs.h"

#define Power(a, b) pow(a, b)
#define sscanf_s sscanf
#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iv]


/* modes */
#define STOCHASTIC 1
#define OUTPUTS 1
#define DISPLAY 1
#define OBSTACLES 0
#define ARRB 1
#define IDM 1

double TIME_STEP = 1.0;
double du = 0.125;
double du_min = 0.125, dx_min = 0.0625;

double X0 = 0.0, V0 = 5.0;

double XMAX = 300.0;
double XMIN =  0.0;
double VMAX = 16.0;
double VMIN = 0.0;

double *x_min, *x_max;
double u_min = -2.0, u_max = 1.0;

/* analytic parameters  */
double ve = 11.0, xe = XMAX + 70.0;
double w = 0.1;

int newtonIntervals = 1;
double newton_add_term = 5.0;
int negSpeedInitCon;

typedef struct node {
	int k;
	int ix, iv;
	double x, v;
	double _x, _v;
} node_t;

typedef struct {
	double u;
} ctrl_t;

double ***J, ***J1, ***J_opt, ***te_opt, ***matrix_te, ***matrix_cost_te;
double *escape_idm;
// double **escape_x, **escape_v, **escape_a;
// double **escape_x_idm, **escape_v_idm, **escape_a_idm;
double *escape_arrb, *escape_arrb_idm;
double *arrb_fuel, *arrb_fuel_idm;

ctrl_t ***C, ***C1;
double *p_red, *p_green;
double *opt_x, *opt_v, *opt_u;
double *idm_x, *idm_v, *idm_u;
FILE *f;

/* UP problem */
double te_up;
/* Dv - Dx for the corridor size */
double Dv, Dx;

/* problem */
struct {
	int n;
	double T;
	double x0, y0, v0;
	int numsteps;
	double ***opt_te;

    // read data
    int m_k, m_ix, m_iv;
    int m_k_2, m_ix_2, m_iv_2;
}P = { 0 };

/* traffic signal variables */
struct signal{
	double Tg_min, Tg_max, Tr_min, Tr_max, Tg_s;
}s = {10.0, 30.0, 40.0, 60.0, 30.0};

/* state variable domains */
struct {
	int *NX, *NXlive_min, *NXlive_max;
	double dx, **X;
	double *xmin, *xmax;

	int *NV, *NVlive_min, *NVlive_max;
	double *vmin, *vmax, dv, **V;

	int NU;
	double umin, umax, du;

	int NX_max, NV_max;
}D = { 0 };

/* IDM */
struct {
    double aMax; 
    double s0;
    double b;
    int d;
}idm = { 0 };

static node_t * getnext(node_t here, ctrl_t c, int print);
static int discretizexv(node_t *state, int k);

static void init_D(double T) {

	int ix, iv;
	int i, k;
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* u domain: du = dv/T */
	D.umin = u_min;
	D.umax = u_max;

	D.du = du;
	D.NU = (int)ceil((D.umax - D.umin) / D.du) + 1;

	/* v domain */
	assert((D.vmax = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.vmin = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.NV = (int*)calloc(sizeof(int), P.numsteps)));

	D.dv = D.du;
	int temp_nv_max = 0;
	for (k = 0; k < P.numsteps; k++){
		D.vmin[k] = 0.0;
		D.vmax[k] = VMAX;

		D.NV[k] = (int)ceil((D.vmax[k] - D.vmin[k]) / D.dv) + 1;
		if (D.NV[k] > temp_nv_max)
			temp_nv_max = D.NV[k];
	}
	D.NV_max = temp_nv_max;

	assert((D.V = (double**)calloc(sizeof(double*), D.NV_max)));
	for (i = 0; i < D.NV_max; i++)
		assert((D.V[i] = (double*)calloc(sizeof(double), P.numsteps)));
	for (k = 0; k < P.numsteps; k++) {
		for (iv = 0; iv < D.NV[k]; iv++)
			D.V[iv][k] = D.vmin[k] + D.dv * iv;
	}

	/* x domain */
	assert((D.xmax = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.xmin = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.NX = (int*)calloc(sizeof(int), P.numsteps)));

	D.dx = 0.5 * D.du;
	int temp_nx_max = 0;
	for (k = 0; k < P.numsteps; k++) {
        if (k >= 0 && k <= s.Tg_max){
			D.xmin[k] = 0.0;
			D.xmax[k] = XMAX + (VMAX);

            D.NX[k] = (int)ceil((D.xmax[k] - D.xmin[k]) / D.dx) + 1;
            if (D.NX[k] > temp_nx_max) {
                temp_nx_max = D.NX[k];
            }
        }
        else {
			D.xmin[k] = 0.0;
			D.xmax[k] = XMAX;
			
            D.NX[k] = (int)ceil((D.xmax[k] - D.xmin[k]) / D.dx) + 1;
            if (D.NX[k] > temp_nx_max)
                temp_nx_max = D.NX[k];
        }
	}
	D.NX_max = temp_nx_max;
	
	assert((D.X = (double**)calloc(sizeof(double*), D.NX_max)));
	for (i = 0; i < D.NX_max; i++)
		assert((D.X[i] = (double*)calloc(sizeof(double), P.numsteps)));
	for (k = 0; k < P.numsteps; k++) {
		for (ix = 0; ix < D.NX[k]; ix++) 
			D.X[ix][k] = D.xmin[k] + D.dx * ix;
	}

	assert((D.NXlive_min = (int*)calloc(sizeof(int), P.numsteps)));
	assert((D.NXlive_max = (int*)calloc(sizeof(int), P.numsteps)));
	assert((D.NVlive_min = (int*)calloc(sizeof(int), P.numsteps)));
	assert((D.NVlive_max = (int*)calloc(sizeof(int), P.numsteps)));

	for (k = 0; k < P.numsteps; k++) {
		D.NXlive_min[k] = 0;
		D.NXlive_max[k] = D.NX[k];
		D.NVlive_min[k] = 0;
		D.NVlive_max[k] = D.NV[k];
	}
}

static void init_JC(void) {
	assert((J = (double ***)calloc(sizeof(double**), P.numsteps)));
	assert((J1 = (double ***)calloc(sizeof(double**), P.numsteps)));
	assert((C = (ctrl_t ***)calloc(sizeof(ctrl_t**), P.numsteps)));
	assert((C1 = (ctrl_t ***)calloc(sizeof(ctrl_t**), P.numsteps)));
	assert((J_opt = (double ***)calloc(sizeof(double**), P.numsteps)));

	assert((matrix_cost_te = (double ***)calloc(sizeof(double**), P.numsteps)));
    assert((matrix_te = (double ***)calloc(sizeof(double**), P.numsteps)));
    assert((te_opt = (double ***)calloc(sizeof(double**), P.numsteps)));
	int k;
    for (k = 0; k < P.numsteps; k++) {
		assert((J[k] = (double **)calloc(sizeof(double*), D.NX_max)));
		assert((J1[k] = (double **)calloc(sizeof(double*), D.NX_max)));
		assert((C[k] = (ctrl_t **)calloc(sizeof(ctrl_t*), D.NX_max)));
		assert((C1[k] = (ctrl_t **)calloc(sizeof(ctrl_t*), D.NX_max)));
		assert((J_opt[k] = (double **)calloc(sizeof(double*), D.NX_max)));

		assert((matrix_cost_te[k] = (double **)calloc(sizeof(double*), D.NX_max)));
        assert((matrix_te[k] = (double **)calloc(sizeof(double*), D.NX_max)));
		assert((te_opt[k] = (double **)calloc(sizeof(double*), D.NX_max)));
		int ix;
		for (ix = 0; ix < D.NX_max; ix++) {
			assert((J[k][ix] = (double *)calloc(sizeof(double), D.NV_max)));
			assert((J1[k][ix] = (double *)calloc(sizeof(double), D.NV_max)));
			assert((C[k][ix] = (ctrl_t *)calloc(sizeof(ctrl_t), D.NV_max)));
			assert((C1[k][ix] = (ctrl_t *)calloc(sizeof(ctrl_t), D.NV_max)));
			assert((J_opt[k][ix] = (double *)calloc(sizeof(double), D.NV_max)));

			assert((matrix_cost_te[k][ix] = (double *)calloc(sizeof(double), D.NV_max)));
		    assert((matrix_te[k][ix] = (double *)calloc(sizeof(double), D.NV_max)));
			assert((te_opt[k][ix] = (double *)calloc(sizeof(double), D.NV_max)));
		}
	}
	assert((p_red = (double *)calloc(sizeof(double), P.numsteps)));
	assert((p_green = (double *)calloc(sizeof(double), P.numsteps)));

	assert((opt_x = (double*)calloc(sizeof(double), P.numsteps)));
	assert((opt_v = (double*)calloc(sizeof(double), P.numsteps)));
	assert((opt_u = (double*)calloc(sizeof(double), P.numsteps)));

	assert((idm_x = (double*)calloc(sizeof(double), P.numsteps)));
	assert((idm_v = (double*)calloc(sizeof(double), P.numsteps)));
	assert((idm_u = (double*)calloc(sizeof(double), P.numsteps)));

	assert((escape_idm = (double *)calloc(sizeof(double), P.numsteps)));

	assert((escape_arrb = (double *)calloc(sizeof(double), P.numsteps)));
	assert((escape_arrb_idm = (double *)calloc(sizeof(double), P.numsteps)));
	assert((arrb_fuel = (double *)calloc(sizeof(double), P.numsteps)));
	assert((arrb_fuel_idm = (double *)calloc(sizeof(double), P.numsteps)));
}

static void free_JC(void) {
	int k, ix, iv;
	for (k = 0; k < P.numsteps; k++) {
		for (ix = 0; ix < D.NX_max; ix++) {
			free(J[k][ix]);
			free(J1[k][ix]);
			free(C[k][ix]);
			free(C1[k][ix]);
			free(J_opt[k][ix]);

			free(te_opt[k][ix]);
			free(matrix_te[k][ix]);
			free(matrix_cost_te[k][ix]);
		}
	}

	for (k = 0; k < P.numsteps; k++) {
		free(J[k]);
		free(J1[k]);
		free(C[k]);
		free(C1[k]);
		free(J_opt[k]);

		free(te_opt[k]);
		free(matrix_te[k]);
		free(matrix_cost_te[k]);
	}
	free(J);
	free(J1);
	free(C);
	free(C1);
	free(J_opt);

	free(te_opt);
	free(matrix_te);
	free(matrix_cost_te);
	//....................................//
	free(D.vmax);
    free(D.vmin);
    free(D.NV);
    free(D.xmax);
    free(D.xmin);
    free(D.NX);
    free(p_red);
	free(p_green);

	free(opt_x);
	free(opt_v);
	free(opt_u);

	free(idm_x);
	free(idm_v);
	free(idm_u);

    for (ix = 0; ix < D.NX_max; ix++)
        free(D.X[ix]);
    free(D.X);

    for (iv = 0; iv < D.NV_max; iv++)
        free(D.V[iv]);
    free(D.V);
}

static int discretizexv(node_t *state, int k) {
	int ix, iv;

	ix = (int)round((state->x - D.X[0][k]) / D.dx);
	iv = (int)round((state->v - D.V[0][k]) / D.dv);
	// fprintf(stderr, "ix: %d -- x: %.3f -- D.X[0][%d]: %.3f \n", ix, state->x, k, D.X[0][k]);
	if (!(ix >= 0 && ix < D.NX[k]))
		return 1;

	if (!(iv >= 0 && iv < D.NV[k]))
		return 2;

	state->x = D.X[(state->ix = ix)][k];
	state->v = D.V[(state->iv = iv)][k];
	//fprintf(stderr, "next.k: %d -- next.x: %.3f -- next.v: %.3f -- D.X[%d][%d]: %.3f \n", k, state->x, state->v, ix, k, D.X[ix][k]);
	return 0;
}

static node_t *getnext(node_t here, ctrl_t c, int print) {
	static node_t next;
	int err;

	next.k = here.k + 1;
	next.x = next._x = here.x + here.v*P.T + 0.5*c.u*pow(P.T, 2);
	next.v = next._v = here.v + c.u*P.T;
	if ((err = discretizexv(&next, next.k))) {
		if (print == 1)
			fprintf(stderr, "err: %d (x: %.4f | v: %.4f)\n", err, next.x, next.v);
		return NULL;
	}

	return &next;
}

static double UP_control(double t, double t0, double x0, double v0, double te) {

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

double te_func(double t0, double x0, double v0, double te) {
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

double dte_func(double t0, double x0, double v0, double te) {
	
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

double arrb(double *v, double *a, int K, double T) {

    double vk, ak, Rt, fuel;
    double m = 1600.0, alpha = 0.666, beta1 = 0.0717, beta2 = 0.0344, b1 = 0.269, b2 = 0.0171, b3 = 0.000672;

    fuel = 0.0;
    for (int k = 0; k < K; k++) {
        vk = v[k]; ak = a[k];
        Rt = b1 + b2 * vk + b3 * pow(vk, 2) + m * ak / 1000.0;
        if (Rt > 0.0) {
            if (ak > 0.0)
                fuel += alpha + beta1 * Rt * vk + (beta2 * m * pow(ak, 2) * vk / 1000.0);
            else
                fuel += alpha + beta1 * Rt * vk;
        }
        else{
            fuel += alpha;
        }
		// fprintf(stderr, "(%d) arrb fuel: %.4f \n", k, fuel * P.T);
    }
    return fuel * T;

}

double arrb_instant(double v, double a) {

    double vk, ak, Rt, fuel;
    double m = 1600.0, alpha = 0.666, beta1 = 0.0717, beta2 = 0.0344, b1 = 0.269, b2 = 0.0171, b3 = 0.000672;

    fuel = 0.0;
    
	vk = v; ak = a;
	Rt = b1 + b2 * vk + b3 * pow(vk, 2) + m * ak / 1000.0;
	if (Rt > 0.0) {
		if (ak > 0.0)
			fuel += alpha + beta1 * Rt * vk + (beta2 * m * pow(ak, 2) * vk / 1000.0);
		else
			fuel += alpha + beta1 * Rt * vk;
	}
	else{
		fuel += alpha;
	}
		
    return fuel;

}

double idm_accel(double x, double xl, double v, double vl, double v0, double ts, int leaderID) {
    
    double s = x - xl;
    double sStar = idm.s0 + fmax(v*ts + ((v * (v - vl)) / (2.0 * sqrt(idm.aMax * idm.b))), 0.0);
    double aFree = idm.aMax * (1.0 - pow(v/v0, idm.d));

    double accIDM, accIDM_x, accIDM_v;
    if (leaderID == -1)
        accIDM = aFree;
    else if (leaderID >= 0)
        accIDM = aFree - (idm.aMax * pow(sStar/s, 4));
    else if (leaderID == -2) {
		double a1 = 1.0 - pow(x/xl, 4);
		double a2 = 1.0 - pow(v/v0, 4);
		
		if (a1 == 0.0 && a2 == 0.0)
			accIDM = 0.0;
		else if (a1 == 0.0 && a2 != 0.0)
			accIDM = idm.aMax * a2;
		else if (a1 != 0.0 && a2 == 0.0)
			accIDM = idm.aMax * a1;
		else
			accIDM = idm.aMax * fmax(a1, a2);
    }
    else
        fprintf(stderr, "Unknown Leader Type. \n");
    
    
    return accIDM;
}

static void printsol(node_t initial, int it) {

	node_t *nextptr;
	int k = 0;
	node_t *state;
	assert((state = (node_t*)calloc(sizeof(node_t), P.numsteps+1)));

	/* trace sequence of states */
	state[0] = initial;
	ctrl_t c, c1;

	#if (OUTPUTS)
		f = fopen("output_mfts.txt","wb");
	#endif
	/* ~~~~~~~~~~~~~~~~~~~~~~~~ */
	if (DISPLAY) {
		fprintf(stderr, "\nk \t Position \t Speed    \t Control  \t Cost \t \t Hybrid  \t J_opt  \t te \t \t p_red \t \t p_green \t (1/2u^2) \n");
		fprintf(stderr, "----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	}

	int actual_T_max;
	double hybrid_cost = 0.0;
	double real_cost = 0.0;
	double real_escape_cost = 0.0; // 0.5u^2 + the escape cost in each stage
	double fuel_sum = 0.0;
	double accel;
	for (k = 0; k < P.numsteps; k++) {
		c = ACCESS(C, state[k]);
		c1 = ACCESS(C1, state[k]);

		accel = c.u;
		// fprintf(stderr, "(%d) u: %.4f \t u1: %.4f \t J: %.4f \t J1: %.4f \n", k, c.u, c1.u, ACCESS(J, state[k]), ACCESS(J1, state[k]));

		/* real cost*/
		real_cost += 0.5 * pow(c.u, 2);
		real_escape_cost = real_cost + ACCESS(J_opt, state[k]);
		
		fuel_sum += arrb_instant(state[k].v, c.u);
		arrb_fuel[k] = fuel_sum;

		if (k < s.Tg_max || k >= s.Tr_min) {
			
			double te;

			double x0 = state[k].x;
			double v0 = state[k].v;
			double t0 = k;
			int numsteps = P.numsteps;
			
			te = ACCESS(te_opt, state[k]);

			// calculate the optimal control cost
			double step;
			double esc_arrb = 0.0;
			for (int i = 0; i < 50000; i++) {
				step = (te - t0) / (50000 - 1);
				double t = (i*step + t0);
				
				esc_arrb += arrb_instant(UP_speed(t, t0, x0, v0, te), UP_control(t, t0, x0, v0, te));
			}
			// fprintf(stderr, "esc_cost: %.4f\n", esc_arrb*step);
			escape_arrb[k] = esc_arrb*step;
		}

		if (k == P.numsteps - 1)
		    ACCESS(J, state[k + 1]) = ACCESS(J, state[k]);

		//fprintf(stderr, "J: %.4f -- J1: %.4f\n", ACCESS(J, state[k]), ACCESS(J1, state[k]));

		if (DISPLAY) {
			fprintf(stderr, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n",
				k,
				state[k].x,
				state[k].v,
				accel,
				ACCESS(J, state[k]),
				hybrid_cost,
				ACCESS(J_opt, state[k]),
				ACCESS(te_opt, state[k]),
				p_red[k],
				p_green[k],
                0.5 * pow(c.u, 2)            
				);
		}
		#if (OUTPUTS)
			fprintf(f, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n",
				k,
				state[k].x,
				state[k].v,
				c.u,
				ACCESS(J, state[k]),
				hybrid_cost,
				ACCESS(J_opt, state[k]),
				ACCESS(te_opt, state[k]),
				p_red[k],
				p_green[k],
				real_cost,
				real_escape_cost,
				s.Tg_min,
				s.Tg_max,
				s.Tr_min,
				s.Tr_max,
				s.Tg_s,
				arrb_fuel[k],
				escape_arrb[k]
			);
		#endif

		actual_T_max = k;
		/* break at green pass */
		if ( k <= s.Tr_max && state[k].x > XMAX) {
			opt_x[k] = state[k].x;
			opt_v[k] = state[k].v;
			opt_u[k] = c.u;
			
			break;
		}

		if (k < P.numsteps - 1) {
			if ((nextptr = getnext(state[k], c, 1)) == NULL) {
				fprintf(stderr, "invalid \n");
				exit(1);
			}
			state[k + 1] = *nextptr;
		}

		opt_x[k] = state[k].x;
		opt_v[k] = state[k].v;
		opt_u[k] = c.u;
	}
	if (DISPLAY)
		fprintf(stderr, "----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");

	FILE *f_det;
	f_det = fopen("output_deter_mfts.txt", "w");
	double up_a, up_v, up_x, step, time;
	double fuel_det = 0.0;
	// fprintf(stderr, "Dt: %.4f - te: %.4f - st: %.4f \n", (ACCESS(te_opt, state[actual_T_max]) - (double)(actual_T_max)), ACCESS(te_opt, state[actual_T_max]), (double)(actual_T_max));
	for (int i = 0; i < P.numsteps; i++) {
		step = (ACCESS(te_opt, state[actual_T_max]) - (double)(actual_T_max)) / (double)(P.numsteps - 1);
		time = (double)(actual_T_max) + i*step;
 		
		// fprintf(stderr, "step: %.4f - time: %.4f \n", step, time);
		up_a = UP_control(time, actual_T_max, state[actual_T_max].x, state[actual_T_max].v, ACCESS(te_opt, state[actual_T_max]));
		up_v = UP_speed(time, actual_T_max, state[actual_T_max].x, state[actual_T_max].v, ACCESS(te_opt, state[actual_T_max]));
		up_x = UP_position(time, actual_T_max, state[actual_T_max].x, state[actual_T_max].v, ACCESS(te_opt, state[actual_T_max]));
		fuel_det += arrb_instant(up_v, up_a);
		// fprintf(stderr, "(%.4f) UP_x: %.4f - UP_v: %.4f - UP_a: %.4f \n", time, up_x, up_v, up_a);
		fprintf(f_det, "%.4f \t %.4f \t %.4f \t %.4f\n", time, up_x, up_v, up_a);
	}
	fuel_det *= step;
	double fuel = arrb(opt_v, opt_u, actual_T_max, P.T);
	fprintf(stderr, "fuel SDP: %.4f -- fuel Det: %.4f -- total: %.4f  with step: %.4f\n", fuel, fuel_det, fuel + fuel_det, step);

	fclose(f_det);
	free(state);
}

static void dp(void) {

	int k, ix, iv;
	double x, v;
	double te;

	fprintf(stderr, "NU: %d -- NV: %d -- NX: %d \n", D.NU, D.NV[2], D.NX[2]);
	double phi, phi1;
	/* for each stage k */
	for (k = P.numsteps - 1; k >= 0; k--) {
		if (DISPLAY)
			fprintf(stderr, "%d ", k);
		
		/* for each possible x */
		for (ix = D.NXlive_min[k]; ix < ((k == 0) ? 1 : D.NXlive_max[k]); ix++) {
			x = ((k == 0) ? P.x0 : D.X[ix][k]);
			/* for each possible v */
			for (iv = D.NVlive_min[k]; iv < ((k == 0) ? 1 : D.NVlive_max[k]); iv++) {
				v = ((k == 0) ? P.v0 : D.V[iv][k]);

				node_t here = { 0 };
				here.k = k; 
				here.x = x;
				here.v = v;
				here.ix = ((k == 0) ? 0 : ix); here.iv = ((k == 0) ? 0 : iv);

				if (here.ix > D.NX[k] || here.iv > D.NV[k]) {
					fprintf(stderr, "invalid state on the grid: k %d x %.1f/%.1f v %.1f/%.1f\n",
						here.k, here.x, D.X[D.NX[k] - 1][k], here.v, D.V[D.NV[k] - 1][k]);
					exit(1);
				}

				node_t next = { 0 };
				double u;
				double *Jhere, *Jhere1, *Jnext, *Jnext1;
				ctrl_t *Chere, *Chere1;

				if (here.k == 0)
					assert(here.ix == 0 && here.iv == 0);

				double result = 0.0;
                /* The 4 phases of the problem */
				//Part 1: uncertain red phase
				if (s.Tr_min <= here.k && here.k <= s.Tr_max) {
                    int ix2 = (int) (here.x / (du * 0.5));
                    int iv2 = (int) (here.v / du);

                    te_opt[here.k][here.ix][here.iv] = matrix_te[0][ix2][iv2] + here.k;
                    result = matrix_cost_te[0][ix2][iv2];
					J_opt[here.k][here.ix][here.iv] = result;

					/* uniform distribution */
					p_red[here.k] = 1.0 / ((s.Tr_max - here.k) + 1.0);
					p_green[here.k] = 0.0;
				}
				//Part 2: red phase
				else if(s.Tg_max < here.k && here.k < s.Tr_min) {					
					// fprintf(stderr, "(%d) phase 2 \n", here.k);
					p_red[here.k] = 0.0;
					p_green[here.k] = 0.0;
				}
				//Part 3.1: last k of uncertain green phase
				else if (here.k == s.Tg_max) {
					// fprintf(stderr, "(%d) phase 3.1 \n", here.k);
                    p_green[here.k] = 0.0;
                    p_red[here.k] = 0.0;
				}
				//Part 3.2: uncertain green phase
				else if (s.Tg_min - 1 <= here.k && here.k < s.Tg_max) {
					// fprintf(stderr, "(%d) phase 3.2 \n", here.k);
                    int ix2 = (int) (here.x / (du * 0.5));
                    int iv2 = (int) (here.v / du);

                    te_opt[here.k][here.ix][here.iv] = matrix_te[0][ix2][iv2] + here.k;
                    result = matrix_cost_te[0][ix2][iv2];
					J_opt[here.k][here.ix][here.iv] = result;

                    /* uniform distribution */
					p_green[here.k] = 1.0 / (((s.Tg_max - 1) - here.k) + 1.0);
					p_red[here.k] = 0.0;
                }
				//Part 4: green phase
                else if (here.k < s.Tg_min - 1) {
					// fprintf(stderr, "(%d) phase 4 \n", here.k);
                    int ix2 = (int) (here.x / (du * 0.5));
                    int iv2 = (int) (here.v / du);

                    te_opt[here.k][here.ix][here.iv] = matrix_te[0][ix2][iv2] + here.k;
                    result = matrix_cost_te[0][ix2][iv2];
					J_opt[here.k][here.ix][here.iv] = result;

                    p_green[here.k] = 0.0;
                    p_red[here.k] = 0.0;
                }
				// unknown phase
				else {
					fprintf(stderr, "Error with time windows during probability calculations\n");
					exit(1);
				}

				/* Calculate Cost */
				/* final stage costs */
				if (here.k == P.numsteps - 1) {
					J[here.k][here.ix][here.iv] = J_opt[here.k][here.ix][here.iv];
					continue;
				}
				
				/* intermediate state */
				Jhere = &J[here.k][here.ix][here.iv];
				Jhere1 = &J1[here.k][here.ix][here.iv];
				Chere = &C[here.k][here.ix][here.iv];
				Chere1 = &C1[here.k][here.ix][here.iv];

				*Jhere = DBL_MAX;
				*Jhere1 = DBL_MAX;
				next.k = here.k + 1;

				for (u = D.umin; u <= D.umax; u += D.du) {
					next.x = here.x + here.v * P.T + 0.5 * u * pow(P.T, 2);
					next.v = here.v + u * P.T;

					if (discretizexv(&next, next.k) != 0) {
						continue;
					}

					Jnext = &J[next.k][next.ix][next.iv];
					Jnext1 = &J1[next.k][next.ix][next.iv];
					
					/* check this */
					double g_offset = s.Tg_max - s.Tg_s;
					// double g_offset = 0.0;
					if ( (next.k < s.Tg_max - g_offset) && (next.x > XMAX) ) {
						*Jnext = J_opt[next.k][next.ix][next.iv];
						*Jnext1 = 0.0;
					}
					
					// cost of red phases
					if (s.Tg_max - g_offset <= here.k && here.k <= s.Tr_max)
						phi = ((1.0 - p_red[here.k]) * (0.5*pow(u, 2) + *Jnext)) + (p_red[here.k] * J_opt[here.k][here.ix][here.iv]);
					else if (here.k < s.Tg_max - g_offset) {
						phi = 0.5*pow(u, 2) + J_opt[here.k][here.ix][here.iv] + p_green[here.k] * (*Jnext1) + (1 - p_green[here.k]) * (*Jnext);
						phi1 = 0.5*pow(u, 2) + *Jnext1;
						if (here.x > XMAX)
							phi1 = DBL_MAX;
					}
					else {
						fprintf(stderr, "Error: Unknown cost calculation\n");
						exit(1);
					}

					if (phi < *Jhere) {
						*Jhere = phi;
						Chere->u = u;
						
						if ((here.k == s.Tg_max - g_offset)) {
							*Jhere1 = phi;
							Chere1->u = u;
						}
					}
					if ( (phi1 < *Jhere1) && (here.k < s.Tg_max - g_offset)) {
						*Jhere1 = phi1;
						Chere1->u = u;
					}

				}
			}
		}
	}
}

void read_data() {
	// int k, ix, iv;
	// double x, v;

	int i;
	int buffSize = 512;
	char buff[512];
    char buff2[512];
	// const char* line;
	double dval;
	int ival, ival2, ival3;

	FILE *f1, *f2, *f3;
	f1 = fopen("j_opt_mfts.txt", "r");
    f2 = fopen("te_opt_mfts.txt", "r");

	while (fgets(buff, sizeof(buff), f1)) {
		// fprintf(stderr, "parsing: %s\n", buff);
		if (buff[0] == '\n')
			break;

		/* read initial states */
		if (sscanf_s(buff, "\"J_opt(%d,%d,%d)\":%lf", &ival, &ival2, &ival3, &dval) == 4) {
			if (ival2 + 1 <= D.NX_max) {
				matrix_cost_te[ival][ival2][ival3] = dval;
				P.m_k = MAX(P.m_k, ival + 1);
				P.m_ix = MAX(P.m_ix, ival2 + 1);
				P.m_iv = MAX(P.m_iv, ival3 + 1);
			}
		}
	}

    while (fgets(buff2, sizeof(buff2), f2)) {
		// fprintf(stderr, "parsing: %s\n", buff2);
		if (buff2[0] == '\n')
			break;

		/* read initial states */
		if (sscanf_s(buff2, "\"te_opt(%d,%d,%d)\":%lf", &ival, &ival2, &ival3, &dval) == 4) {
			if (ival2 + 1 <= D.NX_max) {
				matrix_te[ival][ival2][ival3] = dval;
				P.m_k_2 = MAX(P.m_k_2, ival + 1);
				P.m_ix_2 = MAX(P.m_ix_2, ival2 + 1);
				P.m_iv_2 = MAX(P.m_iv_2, ival3 + 1);
			}
		}
	}

    fclose(f1);
    fclose(f2);
}

int main(int argc, char **argv) {
	int k;

	P.n = 1;
	P.numsteps = (int)s.Tr_max + 1;
	P.T = TIME_STEP;
	P.x0 = X0;
	P.v0 = V0;

	node_t initial = { 0 };
	initial.k = 0;
	initial.x = X0;
	initial.v = V0;
	initial.ix = initial.iv = 0;

	init_D(P.T);
	init_JC();

    read_data();

	/* run dp */
	clock_t start_it, end_it;
	start_it = clock();

	dp();

	end_it = clock();
	double cpu_time_used_it = ((double)(end_it - start_it)) / CLOCKS_PER_SEC;
	
	printsol(initial, 0);

	#if (IDM)
		FILE *f_idm;
		f_idm = fopen("output_idm_mfts.txt", "w");

		idm.aMax = 1.0; idm.b = 1.5; idm.d = 4; 
		idm.s0 = 3.0;

		assert((idm_x = (double*)calloc(sizeof(double), 5*P.numsteps)));
		assert((idm_v = (double*)calloc(sizeof(double), 5*P.numsteps)));
		assert((idm_u = (double*)calloc(sizeof(double), 5*P.numsteps)));

		idm_x[0] = P.x0; idm_v[0] = P.v0;
		double idm_fuel_sum = 0.0;
		int leader = 0;
		for (int k = 0; k < 5*P.numsteps; k++) {
			// calculate IDM trajectory 
			
			if (s.Tg_s == 20.0) {
				// switch at k = 20.0
				idm_u[k] = idm_accel(idm_x[k], XMAX + idm.s0, idm_v[k], 0.0, 12.0, 0.0, 0);
				if (k > P.numsteps - 1) {
					idm_u[k] = idm_accel(idm_x[k], XMAX + 70.0 + idm.s0, idm_v[k], 0.0, 12.0, 0.0, -2);
					// P.T = 0.141;
				}
			}
			else if (s.Tg_s == 30.0) {
				// switch at k = 30.0
				idm_u[k] = idm_accel(idm_x[k], XMAX + 70.0 + idm.s0, idm_v[k], 0.0, 14.0, 0.0, -1);
				if (idm_x[k] > XMAX - 14.0) {
					idm_u[k] = idm_accel(idm_x[k], XMAX + 70.0 + idm.s0, idm_v[k], 0.0, 11.0, 0.0, -1);
					// P.T = 0.141;
				}
			}
			
			// states
			idm_v[k + 1] = idm_v[k] + P.T * idm_u[k];
			idm_x[k + 1] = idm_x[k] + P.T * idm_v[k] + 0.5 * pow(P.T, 2) * idm_u[k];
			double aCon;
			if (idm_v[k + 1] <= 0.0) {
				idm_u[k] = 0.0; aCon = idm.b;
				idm_x[k + 1] = idm_x[k] + pow(idm_v[k], 2) / (2.0 * aCon);
				idm_v[k + 1] = 0.0;
			}

			idm_fuel_sum += arrb_instant(idm_v[k], idm_u[k]) * P.T;
			fprintf(stderr, "idm_x: %.4f \t idm_v: %.4f \t idm_a: %.4f\n", idm_x[k], idm_v[k], idm_u[k]);
			fprintf(f_idm, "%.4f \t %.4f \t %.4f\n", idm_x[k], idm_v[k], idm_u[k]);

			if (idm_x[k + 1] >= XMAX + 70.0 - 1.0) {
				fprintf(f_idm, "%.4f \t %.4f \t %.4f\n", idm_x[k + 1], idm_v[k + 1], 0.0);
				fprintf(stderr, "idm_x: %.4f \t idm_v: %.4f \t idm_a: %.4f\n", idm_x[k + 1], idm_v[k + 1], 0.0);
				break;
			}
		}
		fclose(f_idm);

		fprintf(stderr, "idm fuel: %.4f \n", idm_fuel_sum);
		P.n++;
	#endif

	#if (OUTPUTS)
		fclose(f);
	#endif

	free_JC();
	fprintf(stderr, "OK \n");
	return 0;
}
