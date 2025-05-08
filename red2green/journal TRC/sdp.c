#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <time.h>


#define Power(a, b) pow(a, b)

#define sscanf_s sscanf
#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iv]

#define FIXEDTIME 0
#define STOCHASTIC 1
#define TRIANGULAR 0

#define TIME_STEP 1.0
#define DU 0.125

/* scenario 1 */
double X0 = 0.0, V0 = 11.0;							// stochactic DP initial states

double T_min = 10.0, T_max = 30.0, T_tl = 0.0;  	// possible change times from T_min to T_max
double x_min = 0.0, x_max = 150.0;					// x domain
double v_min = 0.0, v_max = 16.0;					// v domain
double u_min = -3.0, u_max = 3.0;					// u domain
double ve = 11.0, xe = 220.0;						// optimal control final states

double w = 0.1;
double tf = 30.0;
double wopt = 1.0;

double distC = 30.0;
int newtonIntervals = 5;

typedef struct node {
	int k;
	int ix, iv;
	double x, v;
	double _x, _v;
} node_t;

typedef struct {
	double u;
} ctrl_t;

double ***J;
ctrl_t ***C;
double ***J_opt;
double ***escape;
double ***te_opt;
double *p;

/* problem */
struct {
	double T;
	double x0, y0, v0;
	int numsteps, n;
	double ***opt_te;
}P = { 0 };

/* state variable domains */
struct {
	int NX, *NXlive_min, *NXlive_max;
	double dx, *X;

	int NV, *NVlive_min, *NVlive_max;
	double vmin, vmax, dv, *V;

	int NU;
	double umin, umax, du;
}D = { 0 };

static node_t * getnext(node_t here, ctrl_t c);

static int discretizexv(node_t *state);

static void init_D(double T) {

	int ix, iv;
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* v domain: dv is chosen arbitrarily small (i.e. 0.5) */
	D.umin = -3;
	D.umax = +3;

	// D.dv = P.T; 
	D.dv = DU;
	D.vmin = v_min;
	D.vmax = v_max;
	D.NV = (int)ceil(D.vmax / D.dv)+1;
	
	assert((D.V = (double*)calloc(sizeof(double), D.NV)));
	for (iv = 0; iv < D.NV; iv++) {
		D.V[iv] = D.vmin + D.dv * iv;
	}

	/* u domain: du = dv/T */
	D.du = D.dv;
	D.NU = (int)ceil((D.umax - D.umin) / D.du) + 1;

	//D.dx = 0.5 * D.du * pow(T, 2);
	D.dx = 0.5 * D.du;
	double xmax = x_max;
	D.NX = (int)ceil(xmax/D.dx)+1;
	assert((D.X = (double*)calloc(sizeof(double), D.NX)));
	for (ix = 0; ix < D.NX; ix++) {
		D.X[ix] = x_min + D.dx * ix;
	}

	fprintf(stderr, "NV: %d, NU: %d, NX: %d, numsteps: %d \n", D.NV, D.NU, D.NX, P.numsteps);
	int k;
	D.NXlive_min = (int*)calloc(sizeof(int), P.numsteps);
	D.NXlive_max = (int*)calloc(sizeof(int), P.numsteps);

	D.NVlive_min = (int*)calloc(sizeof(int), P.numsteps);
	D.NVlive_max = (int*)calloc(sizeof(int), P.numsteps);
	
	for (k=0; k < P.numsteps; k++) {
		D.NXlive_min[k] = 0;
		D.NXlive_max[k] = D.NX;

		D.NVlive_min[k] = 0;
		D.NVlive_max[k] = D.NV;
	}
}

static void init_JC(void) {

	assert((J = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((C = (ctrl_t ***)calloc(sizeof(ctrl_t*), P.numsteps)));
	assert((J_opt = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((te_opt = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((escape = (double ***)calloc(sizeof(double*), P.numsteps)));
	int k;
	for (k = 0; k < P.numsteps; k++) {
		assert((J[k] = (double **)calloc(sizeof(double*), D.NX)));
		assert((C[k] = (ctrl_t **)calloc(sizeof(ctrl_t*), D.NX)));
		assert((J_opt[k] = (double **)calloc(sizeof(double*), D.NX)));
		assert((escape[k] = (double **)calloc(sizeof(double*), D.NX)));
		assert((te_opt[k] = (double **)calloc(sizeof(double*), D.NX)));
		int ix;
		for (ix = 0; ix < D.NX; ix++) {
			assert((J[k][ix] = (double *)calloc(sizeof(double), D.NV)));
			assert((C[k][ix] = (ctrl_t *)calloc(sizeof(ctrl_t), D.NV)));
			assert((J_opt[k][ix] = (double *)calloc(sizeof(double), D.NV)));
			assert((escape[k][ix] = (double *)calloc(sizeof(double), D.NV)));
			assert((te_opt[k][ix] = (double *)calloc(sizeof(double), D.NV)));
		}
	}
	assert((p = (double *)calloc(sizeof(double), P.numsteps)));
}

static void free_JC(void) {
	int k, ix;
	for (k = 0; k < P.numsteps; k++) {
		for (ix = 0; ix < D.NX; ix++) {
			free(J[k][ix]);
			free(C[k][ix]);
			free(J_opt[k][ix]);
			free(escape[k][ix]);
			free(te_opt[k][ix]);
		}
	}

	for (k = 0; k < P.numsteps; k++) {
		free(J[k]);
		free(C[k]);
		free(J_opt[k]);
		free(escape[k]);
		free(te_opt[k]);
	}
	free(J);
	free(C);
	free(J_opt);
	free(escape);
	free(te_opt);
}

static int discretizexv(node_t *state) {
	int ix, iv;

	ix = (int)round((state->x - D.X[0]) / D.dx);
	iv = (int)round((state->v - D.V[0]) / D.dv);

	if (!(ix >= 0 && ix < D.NX)) 
		return 1;
	
	if (!(iv >= 0 && iv < D.NV))  
		return 2;

	state->x = D.X[(state->ix = ix)];
	state->v = D.V[(state->iv = iv)];

	return 0;
}

static node_t * getnext(node_t here, ctrl_t c) {
	static node_t next;
	int err;

	next.k = here.k + 1;
	next.x = next._x = here.x + here.v*P.T + 0.5*c.u*pow(P.T, 2);
	next.v = next._v = here.v + c.u*P.T;
	if ((err = discretizexv(&next))) 
		return NULL;

	return &next;
}

/* ----------------------- Optimal Control Functions ---------------------------- */
static double oc_control(double t, double t0, double x0, double v0, double te) {
		
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

static double te_func(double t0, double x0, double v0, double te)
{
	return 0.5*w + 0.5*Power((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 +
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
}

static double dte_func(double t0, double x0, double v0, double te)
{
	
	return 1.*((5.999999999999999*te*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
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

}

static double newton_raphson(double t0, double x0, double v0, double te) {
	
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
/* -------------------------------------------------------------------------- */

static void printsol(node_t initial) {
	
	node_t *nextptr;
	int k = 0;
	node_t *state;

	assert((state = (node_t*)calloc(sizeof(node_t), P.numsteps)));

	/* trace sequence of states */
	state[0] = initial;
	ctrl_t c;

	FILE *output; // *output2;
	char buf[0x100];
	if (STOCHASTIC) {
		snprintf(buf, sizeof(buf), "results_sdp/stochastic x0=%.1f - v0=%.1f.txt", P.x0, P.v0);
		output = fopen(buf,"w");
		// output = fopen("results_sdp/output_stochastic.txt", "w+");
		// output2 = fopen("results_sdp/output_oc_stochastic.txt", "w+");
	}
	else {
		snprintf(buf, sizeof(buf), "results_sdp/deterministic x0=%.1f - v0=%.1f.txt", P.x0, P.v0);
		output = fopen(buf,"w");
		// output = fopen("results_sdp/output_deterministic.txt", "w+");
		// output2 = fopen("results_sdp/output_oc_deterministic.txt", "w+");
	}

	fprintf(stderr, "\nk \t Position \t Speed    \t Control  \t Cost \t \t Hybrid  \t J_opt  \t te \t Probability \t Real Cost (1/2u^2) \t Real+escape \t Escape \n");
	fprintf(stderr, "-----------------------------------------------------------------------------------------------------------------\n");

	double hybrid_cost = 0.0;
	double real_cost = 0.0;
	double real_escape_cost = 0.0; // 0.5u^2 + the escape cost in each stage

	//to keep optimal trajectories
	double opt_x[P.numsteps+1], opt_v[P.numsteps+1], opt_u[P.numsteps+1];
	for (k = 0; k < P.numsteps; k++) {
		c = ACCESS(C, state[k]);

		

		/* hybrid cost */
		if (k < T_min){
			hybrid_cost += 0.5 * pow(c.u, 2);
			real_escape_cost += 0.5 * pow(c.u, 2);
		}
		else{
			hybrid_cost += (1.0 - p[k]) * 0.5 * pow(c.u, 2) + p[k] * ACCESS(J_opt, state[k]);
			if (STOCHASTIC) {
				real_escape_cost = real_cost + ACCESS(J_opt, state[k]);
			}	
			else {
				if (k == P.numsteps - 1) {
					real_escape_cost = real_cost + ACCESS(J_opt, state[k]);
				}	
				else {
					real_escape_cost = real_cost + ACCESS(escape, state[k]);
				}
			}
		}
		
		/* real cost*/
		real_cost += 0.5 * pow(c.u, 2);

		if (k == (T_max - T_tl)){
			real_cost += ACCESS(J_opt, state[k]);
			// fprintf(stderr, "\nreal cost: %lf\n", real_cost);
		}
		opt_x[k] = state[k].x;
		opt_v[k] = state[k].v;
		opt_u[k] = c.u;
		fprintf(stderr, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n",
			k, 
			state[k].x,
			state[k].v,
			c.u,
			ACCESS(J, state[k]),
			hybrid_cost,
			ACCESS(J_opt, state[k]),
			ACCESS(te_opt, state[k]),
			p[k],
			real_cost,
			real_escape_cost,
			ACCESS(escape, state[k])
			);
		fprintf(output, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n",
			k, 
			state[k].x,
			state[k].v,
			c.u,
			ACCESS(J, state[k]),
			hybrid_cost,
			ACCESS(J_opt, state[k]),
			ACCESS(te_opt, state[k]),
			p[k],
			real_cost,
			real_escape_cost,
			ACCESS(escape, state[k])
			);

		if (k < P.numsteps - 1) {
			if ((nextptr = getnext(state[k], c)) == NULL) {
				fprintf(stderr, "invalid \n");
				exit(1);
			}
			state[k + 1] = *nextptr;
		}
	}

	double t0 = T_max;
	double te;
	double x0 = state[P.numsteps-1].x, v0 = state[P.numsteps-1].v;
	if (FIXEDTIME) {
		te = tf;
	}
	else {
		te = ACCESS(te_opt, state[P.numsteps-1]);
	}

	double xk = x0, vk = v0;
	int opNumsteps = (int)(te - t0) + 1;
	//double step = (te - t0) / (opNumsteps);
	int count = 0;
	double a_op[100];

	//for (double t = t0; t <= (te); t = t + step) {
	double step = (te - t0) / opNumsteps;
	for (int i = 0; i < opNumsteps; i++) {
		double t = (i*(te - t0) / (opNumsteps - 1) + t0);
		a_op[count] = oc_control(t, t0, x0, v0, te);
		// fprintf(output2, "%.4f \t", a_op[count]);
		// fprintf(stderr, "t: %.4f(%.4f) -- count: %d \n", t, te, count);
		count++;
	}
	
	for (int i = 0; i < count; i++) {
		xk += vk * step + 0.5*pow(step, 2)*a_op[i];
		vk += step * a_op[i];

		// fprintf(output2, "\n %.4f", xk);
		// fprintf(output2, "\t %.4f", vk);
	}

	for (int i = 0; i < P.numsteps; i++)
		fprintf(stderr, "%.4f, ", opt_x[i]);
	fprintf(stderr, "\n\n");

	for (int i = 0; i < P.numsteps; i++)
		fprintf(stderr, "%.4f, ", opt_v[i]);
	fprintf(stderr, "\n\n");

	for (int i = 0; i < P.numsteps; i++)
		fprintf(stderr, "%.4f, ", opt_u[i]);
	fprintf(stderr, "\n\n");

	free(state);
	
	fclose(output);
	// fclose(output2);
}

static void dp(void) {
	
	int k, ix, iv;
	double x, v;
	double te;
	if (STOCHASTIC)
		fprintf(stderr, "Stochastic DP \n");
	else
		fprintf(stderr, "Deterministic DP \n");
	/* for each stage k */
	for (k = P.numsteps - 1; k >= 0; k--) {
		fprintf(stderr, "%d ", k);
		/* for each possible x */
		for (ix = D.NXlive_min[k]; ix < ((k == 0) ? 1 : D.NXlive_max[k]); ix++) {
			x = ((k == 0) ? P.x0 : D.X[ix]);
		// #pragma omp parallel for
			/* for each possible v */
			for (iv = D.NVlive_min[k]; iv < ((k == 0) ? 1 : D.NVlive_max[k]); iv++) {
				v = ((k == 0) ? P.v0 : D.V[iv]);

				node_t here = { 0 };
				here.k = k;
				here.x = x;
				here.v = v;
				here.ix = ((k == 0) ? 0 : ix);
				here.iv = ((k == 0) ? 0 : iv);

				if (here.ix > D.NX || here.iv > D.NV) {
					fprintf(stderr, "invalid state on the grid: k %d x %.1f/%.1f v %.1f/%.1f\n",
						here.k, here.x, D.X[D.NX - 1], here.v, D.V[D.NV - 1]);
					exit(1);
				}

				node_t next = { 0 };
				double u, phi;
				double *Jhere, *Jnext;
				ctrl_t *Chere;
				
				if (here.k == 0)
					assert(here.ix == 0 && here.iv == 0);

				double result = 0.0;
				if ((T_min <= here.k && here.k <= T_max) && (STOCHASTIC || k == P.numsteps - 1)) {
					double te_init, temp_te, check_te, temp_cost = DBL_MAX;

					double max_te = 2.0*(xe - here.x)/0.1; // assuming min speed is 0.1
					te_init = here.k+1;

					for (int i = 0; i < newtonIntervals; i++) {
						double te_cost = 0.0;

						check_te = newton_raphson(here.k, here.x, here.v, te_init);
						double checkSpeed = here.v;
						int negativeSpeed = 0;
						if (check_te > 3.0*max_te || isnan(check_te) || check_te < 0.0)
							continue;

						// fprintf(stderr, "%.3f \n", check_te);
						for (int i = 0; i < P.numsteps; i++) {
							double step = (te - here.k) / (P.numsteps);
							double t = (i*(te - here.k) / (P.numsteps - 1) + here.k);

							te_cost += 0.5*pow(oc_control(t, here.k, here.x, here.v, check_te), 2);
							/* skip optimal te that leads to negative speed trajectories */
							checkSpeed +=  oc_control(t, here.k, here.x, here.v, check_te)*P.T;
							if (checkSpeed < 0)
								negativeSpeed = 1;
						}
						if ((check_te > here.k) && (te_cost < temp_cost) && (negativeSpeed != 1)) {
							temp_cost = te_cost;
							temp_te = check_te;
						}
						te_init += 5.0;
					}
					te = temp_te;
					te_opt[here.k][here.ix][here.iv] = temp_te;

					// calculate the optimal control cost
					double step;
					for (int i = 0; i < P.numsteps; i++) {
						step = (te - here.k) / (P.numsteps);
						double t = (i*(te - here.k) / (P.numsteps - 1) + here.k);
						result += 0.5* (pow(oc_control(t, here.k, here.x, here.v, te), 2));
					}
					result *= step;

					if ((!STOCHASTIC && k == P.numsteps - 1) || STOCHASTIC) {
						J_opt[here.k][here.ix][here.iv] = wopt*result;
					}
					else {
						escape[here.k][here.ix][here.iv] = wopt*result;
						J_opt[here.k][here.ix][here.iv] = 0.0;
					}
		
					if (TRIANGULAR){
						/* triangular distribution */
						if ( here.k >= T_min && here.k < distC )
							p[here.k] = ( 2.0 * (here.k - T_min) ) / ( (T_max - T_min)*(distC - T_min) );
						else if (here.k > distC && here.k <= T_max)
							p[here.k] = ( 2.0 * (T_max - here.k) ) / ( (T_max - T_min)*(T_max - distC) );
						else if (here.k == distC)
							p[here.k] = 2.0 / (T_max - T_min);
					}
					else{
						/* uniform distribution */
						p[here.k] = 1.0 / ((T_max - here.k) + 1.0);
					}

					if (!STOCHASTIC) {
						if (here.k == (int)T_max) {
							p[here.k] = 1.0;
						}
						else {
							p[here.k] = 0.0;
						}	
					}
				}
				else {
					p[here.k] = 0.0;
				}

				if (here.k == P.numsteps - 1) {
				/* final stage costs */
					J[here.k][here.ix][here.iv] = J_opt[here.k][here.ix][here.iv];
					continue;
				}
				
				/* intermediate state */
				Jhere = &J[here.k][here.ix][here.iv];
				Chere = &C[here.k][here.ix][here.iv];
			
				*Jhere = DBL_MAX;
				next.k = here.k + 1;

				/* for each possible u */
				for (u = D.umin; u <= D.umax; u += D.du) {
					next.x = here.x + here.v * P.T + 0.5 * u * pow(P.T, 2);
					next.v = here.v + u * P.T;
						
					if (discretizexv(&next))
						continue;			

					Jnext = &J[next.k][next.ix][next.iv];

					phi = ((1.0 - p[here.k]) * (0.5*pow(u, 2) + *Jnext)) + (p[here.k] * J_opt[here.k][here.ix][here.iv]);
					// phi = 0.5*pow(u, 2) + ((1.0 - p[here.k]) * (*Jnext)) + (p[here.k] * J_opt[here.k][here.ix][here.iv]);
					
					if (phi < *Jhere) {
						*Jhere = phi;
						Chere->u = u;
					}
				}
			}
		}
	}
}

int main(int argc, char **argv) {
	
	P.numsteps = (int)((T_max - T_tl)) + 1;
	P.T = TIME_STEP;

	// double scenario_x[] = { 0.0,  10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 };
  	// double scenario_v[] = { 1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0 };
	
	double scenario_x[] = { 30.0 };
  	double scenario_v[] = { 9.0, 10.0, 11.0 };

	for (int xi = 0; xi < 1; xi++) {
		for (int vi = 0; vi < 3; vi++) {
			
			P.x0 = scenario_x[xi]; P.v0 = scenario_v[vi];
			// P.x0 = X0; P.v0 = V0;
			
			init_D(P.T);
			init_JC();

			node_t initial = { 0 };
			initial.k = 0;
			initial.x = P.x0;
			initial.v = P.v0;
			initial.ix = initial.iv = 0;
			
			/* run dp */
			double start, end;
			start = clock();
				
			dp();
			
			end = clock();
			double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
			double complexity = (P.numsteps*D.NV*D.NX*D.NU) * pow(10, -6);
			
			/* prints */
			printsol(initial);
			fprintf(stderr, "\nDU: %.2f -- DV: %.2f DX: %.2f \n", D.du, D.dv, D.dx);
			fprintf(stderr, "CPU time: %.4f \n", cpu_time_used);
			fprintf(stderr, "Complexity: %.4f \n", complexity);
			free_JC();
		}
	}
	return 0;
}
