#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdbool.h>


#define Power(a, b) pow(a, b)
// #define fopen fopen_s
#define sscanf_s sscanf
//#define fopen_s fopen
#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iv]

#define STOCHASTIC 1
#define TRIANGULAR 0
#define FIXED_TIME 0

#define TIME_STEP 1.0

#define DELTA_V 4.0
#define DELTA_X_mult 5.0

#define TERM_CRITERION 0.1
#define TERM_CRITERION_2 0.00001
#define TERM_CRITERION_DU 1

#define ALLOC_STATES 90

/* modes */
#define DISPLAY 0
#define OUTPUT 0
#define OBSTACLES 0
#define IDM 0
#define ARRB 0

double du, init_du = 0.5;

double XMAX = 150.0, XMIN = 0.0;
double VMAX = 16.0,  VMIN = 0.0;
double T_min = 10.0, T_max = 30.0, T_tl = 0.0;  	// possible change times from T_min to T_max
double *x_min, *x_max;								// x domain (difDP)
double u_min = -3.0, u_max = 3.0;					// u domain

/* analytic parameters  */
double ve = 11.0, xe = 220.0;						// analytic optimal control final states
double w = 0.1;
double tf = 38.0;
double wopt = 1.0;

int newtonIntervals = 1;
double newton_add_term = 5.0;
int negSpeedInitCon;

#define SCENARIO 0
#if (SCENARIO == 1)
/* scenario 1 x0=0 - v0=5 (analytical solution) */
double X0 = 0.0, V0 = 5.0;
double initial_v[] = { 5.0000000000000355, 4.811738999540683, 4.643294946498108, 4.494667840872311, 4.365857682663291, 4.25686447187105, 4.167688208495586, 4.098328892536901, 4.048786523994992, 4.019061102869862, 4.009152629161509, 4.019061102869934, 4.048786523995137, 4.098328892537117, 4.167688208495876, 4.2568644718714115, 4.365857682663725, 4.4946678408728165, 4.643294946498687, 4.811738999541333, 5.0000000000007585, 5.20807794787696, 5.435972843169941, 5.683684685879699, 5.9512134760062345, 6.238559213549548, 6.5457218985096395, 6.872701530886509, 7.219498110680156, 7.586111637890581, 7.972542112517783 };
double initial_x[] = { 0., 4.904218087485628, 9.630083648220293, 14.197413629620769, 18.62602497910384, 22.935734644086278, 27.146359571984867, 31.277716710216378, 35.34962300619759, 39.38189540734529, 43.39435086107624, 47.406806314807234, 51.439078715955034, 55.51098501193643, 59.64234215016819, 63.85296707806712, 68.16267674304994, 72.59128809253347, 77.1586180739345, 81.88448363466978, 86.78870172215609, 91.89108928381022, 97.21146326704894, 102.76964061928904, 108.58543828794726, 114.67867322044042, 121.06916236418527, 127.77672266659863, 134.82117107509723, 142.22232453709785, 150.00000000001734 };
double initial_a[] = { -0.19816947416774147, -0.17835252675096372, -0.15853557933418594, -0.13871863191740819, -0.11890168450063043, -0.09908473708385268, -0.07926778966707491, -0.05945084225029715, -0.039633894833519395, -0.019816947416741643, 0.0, 0.01981694741681389, 0.03963389483359164, 0.05945084225036942, 0.07926778966714718, 0.09908473708392493, 0.11890168450070268, 0.13871863191748043, 0.1585355793342582, 0.17835252675103594, 0.1981694741678137, 0.2179864215845915, 0.23780336900136925, 0.257620316418147, 0.27743726383492473, 0.2972542112517025, 0.31707115866848035, 0.3368881060852581, 0.35670505350203585, 0.3765220009188136, 0.39633894833559136 };
#elif (SCENARIO == 2)
/* scenario 2 x0=0 - v0=11 (analytical solution) */
double X0 = 0.0, V0 = 11.0;
double initial_v[] = { 11, 10.1558, 9.35835, 8.60767, 7.90374, 7.24657, 6.63616, 6.0725, 5.55561, 5.08547, 4.66209, 4.28547, 3.95561,3.6725,3.43616,3.24657,3.10374,3.00767,2.95835,2.9558,3,3.09096,3.22868,3.41316,3.64439,3.92239,4.24714,4.61865,5.03691,5.50194,6.01372};
double initial_x[] = { 0,10.574,20.3272,29.3063,37.5581,45.1294,52.0668,58.4173,64.2274,69.5441,74.4139,78.8838,83.0005,86.8106,90.3611,93.6985,96.8698,99.9216,102.901,105.854,108.828,111.869,115.025,118.342,121.867,125.647,129.728,134.157,138.981,144.246,150};
double initial_a[] = { -0.867582,-0.820823,-0.774065,-0.727307,-0.680549,-0.633791,-0.587033,-0.540274,-0.493516,-0.446758,-0.4,-0.353242,-0.306484,-0.259726,-0.212967,-0.166209,-0.119451,-0.0726928,-0.0259347,0.0208235,0.0675817,0.11434,0.161098,0.207856,0.254614,0.301372,0.348131,0.394889,0.441647,0.488405,0.535163};
#elif (SCENARIO == 3)
/* scenario 3 x0=50 - v0=11 (analytical solution) */
double X0 = 50.0, V0 = 11.0;
double initial_v[] = {11., 9.96957, 8.99381, 8.07271, 7.20627, 6.39451, 5.6374, 4.93496, 4.28719, 3.69408, 3.15564, 2.67186, 2.24274, 1.8683, 1.54851, 1.28339, 1.07294, 0.917152, 0.816029, 0.769571, 0.777778, 0.84065, 0.958187, 1.13039, 1.35726, 1.63879, 1.97498, 2.36585, 2.81137, 3.31157, 3.86642};
double initial_x[] = {50., 60.4802, 69.9574, 78.4861, 86.121, 92.9168, 98.9282, 104.21, 108.816, 112.802, 116.223, 119.132, 121.585, 123.636, 125.34, 126.751, 127.925, 128.915, 129.777, 130.565, 131.334, 132.139, 133.034, 134.074, 135.313, 136.806, 138.609, 140.775, 143.359, 146.416, 150.};
double initial_a[] = {-1.05776, -1.0031, -0.948431, -0.893766, -0.839101, -0.784436, -0.729771, -0.675106, -0.620441, -0.565776, -0.511111, -0.456446, -0.401781, -0.347116, -0.292451, -0.237786, -0.183121, -0.128456, -0.0737907, -0.0191257, 0.0355393, 0.0902044, 0.144869, 0.199534, 0.2542, 0.308865, 0.36353, 0.418195, 0.47286, 0.527525, 0.58219};
#endif

typedef struct node {
	int k;
	int ix, iv;
	double x, v;
	double _x, _v;
} node_t;

typedef struct {
	double u;
} ctrl_t;

double ***J, ***J_opt, ***escape, ***te_opt;
ctrl_t ***C;
double *p;

/* new */
double *opt_x, *opt_v, *opt_u;
double *idm_x, *idm_v, *idm_u;
double opt_J_it;
double *init_x, *init_v, *init_u;
double *initial_x, *initial_v, *initial_u;
// double *initial_x, *initial_v, *initial_a; //fix all the vectors, check for un-needed
/* UP problem */
double te_up;
double Dv, Dx;

/* problem */
struct {
	double T;
	double x0, y0, v0;
	int numsteps;
	double ***opt_te;

	/* obstacles */
	double safety;
	int n;
}P = { 0 };

/* state variable domains */
struct {
	int *NX, *NXlive_min, *NXlive_max;
	double dx, **X;
	double *xmin, *xmax;

	int *NV, *NVlive_min, *NVlive_max;
	double *vmin, *vmax, dv, **V;

	int NU;
	double umin, umax, du;
}D = { 0 };

/* obstacles */
struct {
	double **obst_x;
	double **obst_y;
	double **obst_v;
}O = { 0 };

/* IDM */
struct {
    double aMax; 
    double s0;
    double b;
    int d;
}idm;

static node_t * getnext(node_t here, ctrl_t c);

static int discretizexv(node_t *state, int k);

static double crash(node_t state) {
	int i;
	double dx;
	double sgap;
	double cost = 0;
	double L = 5.0;

	for (i = 0; i < P.n; i++) {
		dx = fabs(O.obst_x[state.k][i] - state.x);

		if (O.obst_x[state.k][i] - state.x >= 0) {
			sgap = (state.v * (P.safety) + L);
		}
		else {
			sgap = (O.obst_v[state.k][i] * (P.safety) + L);
		}
		if (dx <= sgap){
			cost += 1.0;
		}

	}
	return cost;
}

static void init_S(int numsteps, double T, int n, double obst_x, double obst_v)
{
	int i, k;

	/* allocate */
	assert((O.obst_x = (double**)calloc(sizeof(double*), numsteps + 1)));
	assert((O.obst_v = (double**)calloc(sizeof(double*), numsteps + 1)));
	for (k = 0; k < numsteps + 1; k++) {
		assert((O.obst_x[k] = (double*)calloc(sizeof(double), n)));
		assert((O.obst_v[k] = (double*)calloc(sizeof(double), n)));
	}

	for (i = 0; i < n; i++) {
		O.obst_x[0][i] = obst_x;
		O.obst_v[0][i] = obst_v;
	}

	/* compute */
	for (k = 0; k < P.numsteps; k++) {
		for (i = 0; i < n; i++) {
			O.obst_x[k + 1][i] = O.obst_x[k][i] + O.obst_v[k][i] * P.T;
			O.obst_v[k + 1][i] = O.obst_v[k][i];
		}
	}
}

static void init_D(double T, int iteration) {

	int ix, iv;
	int i, k;
  
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	assert((D.V = (double**)calloc(sizeof(double*), ALLOC_STATES)));
	for (i = 0; i < ALLOC_STATES; i++)
		assert((D.V[i] = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.X = (double**)calloc(sizeof(double*), ALLOC_STATES)));
	for (i = 0; i < ALLOC_STATES; i++)
		assert((D.X[i] = (double*)calloc(sizeof(double), P.numsteps)));
  
	assert((D.vmax = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.vmin = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.xmax = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.xmin = (double*)calloc(sizeof(double), P.numsteps)));

	assert((D.NV = (int*)calloc(sizeof(int), P.numsteps)));
	assert((D.NX = (int*)calloc(sizeof(int), P.numsteps)));
  
	/* u domain: du = dv/T */
	D.umin = u_min;
	D.umax = u_max;

	D.du = du;
	D.NU = (int)ceil((D.umax - D.umin) / D.du) + 1;
  
	/* v domain */
	D.dv = D.du;
	Dv = DELTA_V;

	for (k = 0; k < P.numsteps; k++){
		if (iteration == 0) {
			D.vmin[k] = floor(MAX(init_v[k] - Dv * D.dv, VMIN));
			D.vmax[k] = ceil(MIN(init_v[k] + Dv * D.dv, VMAX));

		}
		else {
			D.vmin[k] = MAX(init_v[k] - Dv * D.dv, VMIN);
			D.vmax[k] = MIN(init_v[k] + Dv * D.dv, VMAX);
	
		}
		D.NV[k] = (int)ceil((D.vmax[k] - D.vmin[k]) / D.dv) + 1;
		
		for (iv = 0; iv < D.NV[k]; iv++) {
			D.V[iv][k] = D.vmin[k] + D.dv * iv;
		// fprintf(stderr, "D.V[%d][%d]: %.3f \n", iv, k, D.V[iv][k]);
		}
	}
	
	/* x domain */
	D.dx = 0.5 * D.du;
	Dx = DELTA_V * DELTA_X_mult;
	for (k = 0; k < P.numsteps; k++) {
		if (iteration == 0) {
			D.xmin[k] = floor(MAX(init_x[k] - Dx * D.dx, XMIN));
			D.xmax[k] = ceil(MIN(init_x[k] + Dx * D.dx, XMAX));
	
		}
		else {
			D.xmin[k] = MAX(init_x[k] - Dx * D.dx, XMIN);
			D.xmax[k] = MIN(init_x[k] + Dx * D.dx, XMAX);

		}
		D.NX[k] = (int)ceil((D.xmax[k] - D.xmin[k]) / D.dx) + 1;

		for (ix = 0; ix < D.NX[k]; ix++) {
			D.X[ix][k] = D.xmin[k] + D.dx * ix;
		// fprintf(stderr, "(%d) D.X[%d][%d]: %.3f \n", D.NX[k], ix, k, D.X[ix][k]);
		}
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

	assert((J = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((C = (ctrl_t ***)calloc(sizeof(ctrl_t*), P.numsteps)));
	assert((J_opt = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((te_opt = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((escape = (double ***)calloc(sizeof(double*), P.numsteps)));
	int k;
	for (k = 0; k < P.numsteps; k++) {
		assert((J[k] = (double **)calloc(sizeof(double*), ALLOC_STATES)));
		assert((C[k] = (ctrl_t **)calloc(sizeof(ctrl_t*), ALLOC_STATES)));
		assert((J_opt[k] = (double **)calloc(sizeof(double*), ALLOC_STATES)));
		assert((escape[k] = (double **)calloc(sizeof(double*), ALLOC_STATES)));
		assert((te_opt[k] = (double **)calloc(sizeof(double*), ALLOC_STATES)));
		int ix;
		for (ix = 0; ix < ALLOC_STATES; ix++) {
			assert((J[k][ix] = (double *)calloc(sizeof(double), ALLOC_STATES)));
			assert((C[k][ix] = (ctrl_t *)calloc(sizeof(ctrl_t), ALLOC_STATES)));
			assert((J_opt[k][ix] = (double *)calloc(sizeof(double), ALLOC_STATES)));
			assert((escape[k][ix] = (double *)calloc(sizeof(double), ALLOC_STATES)));
			assert((te_opt[k][ix] = (double *)calloc(sizeof(double), ALLOC_STATES)));
		}
	}
	assert((p = (double *)calloc(sizeof(double), P.numsteps)));
}

static void free_JC(void) {
	int k, ix;
	for (k = 0; k < P.numsteps; k++) {
		for (ix = 0; ix < ALLOC_STATES; ix++) {
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

static int discretizexv(node_t *state, int k) {
	int ix, iv;

	ix = (int)round((state->x - D.X[0][k]) / D.dx);
	iv = (int)round((state->v - D.V[0][k]) / D.dv);
	//fprintf(stderr, "ix: %d -- x: %.3f -- D.X[0][%d]: %.3f \n", ix, state->x, k, D.X[0][k]);
	if (!(ix >= 0 && ix < D.NX[k])) 
		return 1;
	
	if (!(iv >= 0 && iv < D.NV[k]))  
		return 2;

	state->x = D.X[(state->ix = ix)][k];
	state->v = D.V[(state->iv = iv)][k];
	//fprintf(stderr, "next.k: %d -- next.x: %.3f -- next.v: %.3f -- D.X[%d][%d]: %.3f \n", k, state->x, state->v, ix, k, D.X[ix][k]);
	return 0;
}

static node_t * getnext(node_t here, ctrl_t c) {
	static node_t next;
	int err;

	next.k = here.k + 1;
	next.x = next._x = here.x + here.v*P.T + 0.5*c.u*pow(P.T, 2);
	next.v = next._v = here.v + c.u*P.T;
	if ((err = discretizexv(&next, next.k))) 
		return NULL;

	return &next;
}

/* ----------------------- Optimal Control Functions ---------------------------- */
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

static double UP_position(double t, double t0, double x0, double v0, double te) {

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

static double UP_speed(double t, double t0, double x0, double v0, double te) {
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

static double te_func(double t0, double x0, double v0, double te)
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

static double dte_func(double t0, double x0, double v0, double te)
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

static double te_func_cp(double t0, double x0, double v0, double te, double t1, double xred) {

	double res = 0.5*w + 0.5*Power((0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 
           10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 
           30.*Power(t0,3)*Power(t1,4)*te*v0 + 30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.08333333333333333*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))),2) + 
   ((0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 10.*Power(t0,2)*Power(t1,6)*v0 + 
           5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 30.*Power(t0,3)*Power(t1,4)*te*v0 + 
           30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.08333333333333333*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))))*
    ((-0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 10.*Power(t0,2)*Power(t1,6)*v0 + 
           5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 30.*Power(t0,3)*Power(t1,4)*te*v0 + 
           30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
      (0.08333333333333333*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)))) - 
   (0.08333333333333333*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 
        5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
        40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
        30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
        30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
        40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
        4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
        10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
        5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
        4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
        1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
        5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
        5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
        15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
        42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
        1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
        65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 2.*Power(t0,5)*Power(te,4)*ve - 
        10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
        10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
        18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
        72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
        72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 
        12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 
        48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
        18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
        1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 
        25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 
        80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 
        1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
        100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 
        4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
        20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
        2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
        2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
        18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
        27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
        2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
        4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
        32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
        3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
        12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred)*
      ((0.05555555555555555*(-t1 + te)*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 
             10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 
             30.*Power(t0,3)*Power(t1,4)*te*v0 + 30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 
             3.*Power(t0,5)*t1*Power(te,2)*v0 - 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
             30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 
             1.*Power(t0,5)*Power(te,3)*v0 + 5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
             10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 
             1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 
             10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 
             30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 
             3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 
             30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 
             5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
             5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
             18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
             54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
             36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 
             9.*Power(t1,5)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
             12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
             15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
             6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
             30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
             30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 
             3.*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 
             12.*Power(t0,2)*Power(t1,5)*xred + 3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 
             24.*Power(t0,3)*Power(t1,3)*te*xred + 6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 
             3.*Power(t0,5)*Power(te,2)*xred - 6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
             24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 
             12.*Power(t0,3)*t1*Power(te,3)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
        (0.041666666666666664*Power(-t1 + te,2)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
             10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
             40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
             60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
             20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
             20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
             5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
             4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
             1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
             5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
             5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
             15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
             42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
             1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
             65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
             2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 
             20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 
             12.*Power(t0,3)*Power(t1,5)*x0 + 18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 
             48.*Power(t0,3)*Power(t1,4)*te*x0 - 72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 
             18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 
             72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 
             20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 
             30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 
             10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
             100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 
             11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 
             40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 
             8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 
             18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 
             1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 
             28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 
             7.*Power(t1,6)*Power(te,2)*xred + 4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 
             8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 
             28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 3.*Power(t0,4)*Power(te,4)*xred + 
             12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
           (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
             0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
             0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
             0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
             0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
             0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
             0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
        (0.01388888888888889*(1.*Power(t0,5)*Power(t1,4)*v0 - 5.*Power(t0,4)*Power(t1,5)*v0 + 10.*Power(t0,3)*Power(t1,6)*v0 - 10.*Power(t0,2)*Power(t1,7)*v0 + 
             5.*t0*Power(t1,8)*v0 - 1.*Power(t1,9)*v0 - 3.9999999999999996*Power(t0,5)*Power(t1,3)*te*v0 + 20.*Power(t0,4)*Power(t1,4)*te*v0 - 
             40.*Power(t0,3)*Power(t1,5)*te*v0 + 40.*Power(t0,2)*Power(t1,6)*te*v0 - 20.*t0*Power(t1,7)*te*v0 + 3.999999999999999*Power(t1,8)*te*v0 + 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 - 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 + 59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 - 
             60.00000000000001*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 + 30.000000000000004*t0*Power(t1,6)*Power(te,2)*v0 - 
             6.000000000000001*Power(t1,7)*Power(te,2)*v0 - 3.9999999999999996*Power(t0,5)*t1*Power(te,3)*v0 + 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 - 
             40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 + 40.00000000000001*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 - 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 + 
             4.000000000000001*Power(t1,6)*Power(te,3)*v0 + 1.*Power(t0,5)*Power(te,4)*v0 - 5.*Power(t0,4)*t1*Power(te,4)*v0 + 
             10.*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 - 10.*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 + 5.*t0*Power(t1,4)*Power(te,4)*v0 - 
             1.*Power(t1,5)*Power(te,4)*v0 + 1.*Power(t0,6)*Power(t1,3)*ve - 6.*Power(t0,5)*Power(t1,4)*ve + 14.999999999999998*Power(t0,4)*Power(t1,5)*ve - 
             19.999999999999996*Power(t0,3)*Power(t1,6)*ve + 14.999999999999998*Power(t0,2)*Power(t1,7)*ve - 6.*t0*Power(t1,8)*ve + 1.*Power(t1,9)*ve - 
             3.*Power(t0,6)*Power(t1,2)*te*ve + 18.*Power(t0,5)*Power(t1,3)*te*ve - 44.99999999999999*Power(t0,4)*Power(t1,4)*te*ve + 
             59.99999999999999*Power(t0,3)*Power(t1,5)*te*ve - 44.99999999999999*Power(t0,2)*Power(t1,6)*te*ve + 18.*t0*Power(t1,7)*te*ve - 3.*Power(t1,8)*te*ve + 
             3.*Power(t0,6)*t1*Power(te,2)*ve - 18.*Power(t0,5)*Power(t1,2)*Power(te,2)*ve + 44.99999999999999*Power(t0,4)*Power(t1,3)*Power(te,2)*ve - 
             59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*ve + 44.99999999999999*Power(t0,2)*Power(t1,5)*Power(te,2)*ve - 18.*t0*Power(t1,6)*Power(te,2)*ve + 
             3.*Power(t1,7)*Power(te,2)*ve - 1.*Power(t0,6)*Power(te,3)*ve + 6.*Power(t0,5)*t1*Power(te,3)*ve - 
             14.999999999999998*Power(t0,4)*Power(t1,2)*Power(te,3)*ve + 19.999999999999996*Power(t0,3)*Power(t1,3)*Power(te,3)*ve - 
             14.999999999999998*Power(t0,2)*Power(t1,4)*Power(te,3)*ve + 6.*t0*Power(t1,5)*Power(te,3)*ve - 1.*Power(t1,6)*Power(te,3)*ve - 
             3.*Power(t0,4)*Power(t1,4)*x0 + 12.*Power(t0,3)*Power(t1,5)*x0 - 18.*Power(t0,2)*Power(t1,6)*x0 + 12.*t0*Power(t1,7)*x0 - 3.*Power(t1,8)*x0 + 
             12.*Power(t0,4)*Power(t1,3)*te*x0 - 48.*Power(t0,3)*Power(t1,4)*te*x0 + 72.*Power(t0,2)*Power(t1,5)*te*x0 - 48.*t0*Power(t1,6)*te*x0 + 
             12.*Power(t1,7)*te*x0 - 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 + 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 - 
             108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 + 72.*t0*Power(t1,5)*Power(te,2)*x0 - 18.*Power(t1,6)*Power(te,2)*x0 + 12.*Power(t0,4)*t1*Power(te,3)*x0 - 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 + 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,4)*Power(te,3)*x0 + 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 - 3.*Power(t0,4)*Power(te,4)*x0 + 12.*Power(t0,3)*t1*Power(te,4)*x0 - 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 + 12.*t0*Power(t1,3)*Power(te,4)*x0 - 3.*Power(t1,4)*Power(te,4)*x0 + 3.*Power(t0,6)*Power(t1,2)*xe - 
             18.*Power(t0,5)*Power(t1,3)*xe + 45.*Power(t0,4)*Power(t1,4)*xe - 60.*Power(t0,3)*Power(t1,5)*xe + 45.*Power(t0,2)*Power(t1,6)*xe - 18.*t0*Power(t1,7)*xe + 
             3.*Power(t1,8)*xe - 6.*Power(t0,6)*t1*te*xe + 36.*Power(t0,5)*Power(t1,2)*te*xe - 90.*Power(t0,4)*Power(t1,3)*te*xe + 120.*Power(t0,3)*Power(t1,4)*te*xe - 
             90.*Power(t0,2)*Power(t1,5)*te*xe + 36.*t0*Power(t1,6)*te*xe - 6.*Power(t1,7)*te*xe + 3.*Power(t0,6)*Power(te,2)*xe - 18.*Power(t0,5)*t1*Power(te,2)*xe + 
             45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe - 60.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe + 45.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe - 
             18.*t0*Power(t1,5)*Power(te,2)*xe + 3.*Power(t1,6)*Power(te,2)*xe - 3.*Power(t0,6)*Power(t1,2)*xred + 18.*Power(t0,5)*Power(t1,3)*xred - 
             42.*Power(t0,4)*Power(t1,4)*xred + 48.*Power(t0,3)*Power(t1,5)*xred - 26.999999999999996*Power(t0,2)*Power(t1,6)*xred + 
             5.999999999999999*t0*Power(t1,7)*xred + 6.*Power(t0,6)*t1*te*xred - 36.*Power(t0,5)*Power(t1,2)*te*xred + 78.*Power(t0,4)*Power(t1,3)*te*xred - 
             72.*Power(t0,3)*Power(t1,4)*te*xred + 18.*Power(t0,2)*Power(t1,5)*te*xred + 11.999999999999998*t0*Power(t1,6)*te*xred - 
             5.999999999999999*Power(t1,7)*te*xred - 3.*Power(t0,6)*Power(te,2)*xred + 18.*Power(t0,5)*t1*Power(te,2)*xred - 
             26.999999999999996*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 12.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred + 63.*Power(t0,2)*Power(t1,4)*Power(te,2)*xred - 
             54.00000000000001*t0*Power(t1,5)*Power(te,2)*xred + 15.*Power(t1,6)*Power(te,2)*xred - 12.*Power(t0,4)*t1*Power(te,3)*xred + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred + 48.00000000000001*t0*Power(t1,4)*Power(te,3)*xred - 
             12.000000000000002*Power(t1,5)*Power(te,3)*xred + 3.*Power(t0,4)*Power(te,4)*xred - 12.*Power(t0,3)*t1*Power(te,4)*xred + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred - 12.*t0*Power(t1,3)*Power(te,4)*xred + 3.*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)))))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
      (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
        0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
        0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
        0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
        0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
        0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
        0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)));

	return res;
}

static double dte_func_cp(double t0, double x0, double v0, double te, double t1, double xred) {

	double res = 1.*((0.05555555555555555*(-3.*Power(t0,5)*Power(t1,2)*v0 + 15.*Power(t0,4)*Power(t1,3)*v0 - 30.*Power(t0,3)*Power(t1,4)*v0 + 30.*Power(t0,2)*Power(t1,5)*v0 - 
           15.*t0*Power(t1,6)*v0 + 3.*Power(t1,7)*v0 + 6.*Power(t0,5)*t1*te*v0 - 30.*Power(t0,4)*Power(t1,2)*te*v0 + 60.*Power(t0,3)*Power(t1,3)*te*v0 - 
           60.*Power(t0,2)*Power(t1,4)*te*v0 + 30.000000000000004*t0*Power(t1,5)*te*v0 - 6.000000000000001*Power(t1,6)*te*v0 - 3.*Power(t0,5)*Power(te,2)*v0 + 
           15.*Power(t0,4)*t1*Power(te,2)*v0 - 30.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*v0 + 30.000000000000007*Power(t0,2)*Power(t1,3)*Power(te,2)*v0 - 
           15.000000000000005*t0*Power(t1,4)*Power(te,2)*v0 + 3.0000000000000013*Power(t1,5)*Power(te,2)*v0 + 3.*Power(t0,5)*Power(t1,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*ve + 30.*Power(t0,3)*Power(t1,4)*ve - 30.*Power(t0,2)*Power(t1,5)*ve + 15.*t0*Power(t1,6)*ve - 3.*Power(t1,7)*ve - 
           6.*Power(t0,5)*t1*te*ve + 30.*Power(t0,4)*Power(t1,2)*te*ve - 60.*Power(t0,3)*Power(t1,3)*te*ve + 60.*Power(t0,2)*Power(t1,4)*te*ve - 
           30.*t0*Power(t1,5)*te*ve + 6.*Power(t1,6)*te*ve + 3.*Power(t0,5)*Power(te,2)*ve - 15.*Power(t0,4)*t1*Power(te,2)*ve + 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*ve - 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*ve + 15.*t0*Power(t1,4)*Power(te,2)*ve - 3.*Power(t1,5)*Power(te,2)*ve + 
           9.*Power(t0,4)*Power(t1,2)*x0 - 36.*Power(t0,3)*Power(t1,3)*x0 + 54.*Power(t0,2)*Power(t1,4)*x0 - 36.*t0*Power(t1,5)*x0 + 9.*Power(t1,6)*x0 - 
           18.*Power(t0,4)*t1*te*x0 + 72.*Power(t0,3)*Power(t1,2)*te*x0 - 108.*Power(t0,2)*Power(t1,3)*te*x0 + 72.*t0*Power(t1,4)*te*x0 - 18.*Power(t1,5)*te*x0 + 
           9.*Power(t0,4)*Power(te,2)*x0 - 36.*Power(t0,3)*t1*Power(te,2)*x0 + 54.*Power(t0,2)*Power(t1,2)*Power(te,2)*x0 - 
           36.00000000000001*t0*Power(t1,3)*Power(te,2)*x0 + 9.000000000000004*Power(t1,4)*Power(te,2)*x0 + 6.*Power(t0,5)*t1*xe - 30.*Power(t0,4)*Power(t1,2)*xe + 
           60.*Power(t0,3)*Power(t1,3)*xe - 60.*Power(t0,2)*Power(t1,4)*xe + 30.*t0*Power(t1,5)*xe - 6.*Power(t1,6)*xe - 6.*Power(t0,5)*te*xe + 
           30.*Power(t0,4)*t1*te*xe - 60.*Power(t0,3)*Power(t1,2)*te*xe + 60.*Power(t0,2)*Power(t1,3)*te*xe - 30.*t0*Power(t1,4)*te*xe + 6.*Power(t1,5)*te*xe - 
           6.*Power(t0,5)*t1*xred + 21.*Power(t0,4)*Power(t1,2)*xred - 24.*Power(t0,3)*Power(t1,3)*xred + 6.*Power(t0,2)*Power(t1,4)*xred + 6.*t0*Power(t1,5)*xred - 
           3.*Power(t1,6)*xred + 6.*Power(t0,5)*te*xred - 12.*Power(t0,4)*t1*te*xred - 12.000000000000002*Power(t0,3)*Power(t1,2)*te*xred + 
           48.*Power(t0,2)*Power(t1,3)*te*xred - 42.*t0*Power(t1,4)*te*xred + 12.*Power(t1,5)*te*xred - 9.*Power(t0,4)*Power(te,2)*xred + 
           36.*Power(t0,3)*t1*Power(te,2)*xred - 54.*Power(t0,2)*Power(t1,2)*Power(te,2)*xred + 36.00000000000001*t0*Power(t1,3)*Power(te,2)*xred - 
           9.000000000000004*Power(t1,4)*Power(te,2)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.08333333333333333*(-t1 + te)*(4.*Power(t0,5)*Power(t1,3)*v0 - 20.*Power(t0,4)*Power(t1,4)*v0 + 40.*Power(t0,3)*Power(t1,5)*v0 - 
           40.*Power(t0,2)*Power(t1,6)*v0 + 20.*t0*Power(t1,7)*v0 - 4.*Power(t1,8)*v0 - 12.*Power(t0,5)*Power(t1,2)*te*v0 + 60.*Power(t0,4)*Power(t1,3)*te*v0 - 
           120.*Power(t0,3)*Power(t1,4)*te*v0 + 120.*Power(t0,2)*Power(t1,5)*te*v0 - 60.*t0*Power(t1,6)*te*v0 + 12.*Power(t1,7)*te*v0 + 
           12.*Power(t0,5)*t1*Power(te,2)*v0 - 60.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 120.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
           120.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 60.000000000000014*t0*Power(t1,5)*Power(te,2)*v0 - 12.000000000000004*Power(t1,6)*Power(te,2)*v0 - 
           4.*Power(t0,5)*Power(te,3)*v0 + 20.*Power(t0,4)*t1*Power(te,3)*v0 - 40.00000000000001*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
           40.00000000000001*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 20.000000000000007*t0*Power(t1,4)*Power(te,3)*v0 + 4.000000000000003*Power(t1,5)*Power(te,3)*v0 - 
           3.0000000000000004*Power(t0,6)*Power(t1,2)*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*ve - 5.000000000000001*Power(t0,4)*Power(t1,4)*ve - 
           20.*Power(t0,3)*Power(t1,5)*ve + 35.*Power(t0,2)*Power(t1,6)*ve - 22.*t0*Power(t1,7)*ve + 5.*Power(t1,8)*ve + 6.000000000000001*Power(t0,6)*t1*te*ve - 
           12.000000000000002*Power(t0,5)*Power(t1,2)*te*ve - 30.*Power(t0,4)*Power(t1,3)*te*ve + 120.*Power(t0,3)*Power(t1,4)*te*ve - 
           150.*Power(t0,2)*Power(t1,5)*te*ve + 84.*t0*Power(t1,6)*te*ve - 18.*Power(t1,7)*te*ve - 3.000000000000001*Power(t0,6)*Power(te,2)*ve - 
           5.999999999999999*Power(t0,5)*t1*Power(te,2)*ve + 75.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 180.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 
           195.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 102.*t0*Power(t1,5)*Power(te,2)*ve + 21.000000000000004*Power(t1,6)*Power(te,2)*ve + 
           8.*Power(t0,5)*Power(te,3)*ve - 40.*Power(t0,4)*t1*Power(te,3)*ve + 80.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 80.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           40.*t0*Power(t1,4)*Power(te,3)*ve - 8.*Power(t1,5)*Power(te,3)*ve - 12.*Power(t0,4)*Power(t1,3)*x0 + 48.*Power(t0,3)*Power(t1,4)*x0 - 
           72.*Power(t0,2)*Power(t1,5)*x0 + 48.*t0*Power(t1,6)*x0 - 12.*Power(t1,7)*x0 + 36.*Power(t0,4)*Power(t1,2)*te*x0 - 144.*Power(t0,3)*Power(t1,3)*te*x0 + 
           216.*Power(t0,2)*Power(t1,4)*te*x0 - 144.*t0*Power(t1,5)*te*x0 + 36.*Power(t1,6)*te*x0 - 36.00000000000001*Power(t0,4)*t1*Power(te,2)*x0 + 
           144.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 216.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 144.*t0*Power(t1,4)*Power(te,2)*x0 - 
           36.00000000000001*Power(t1,5)*Power(te,2)*x0 + 12.*Power(t0,4)*Power(te,3)*x0 - 48.*Power(t0,3)*t1*Power(te,3)*x0 + 
           72.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,3)*Power(te,3)*x0 + 12.000000000000004*Power(t1,4)*Power(te,3)*x0 - 
           2.0000000000000004*Power(t0,6)*t1*xe + 30.*Power(t0,4)*Power(t1,3)*xe - 80.*Power(t0,3)*Power(t1,4)*xe + 90.*Power(t0,2)*Power(t1,5)*xe - 
           48.*t0*Power(t1,6)*xe + 10.*Power(t1,7)*xe + 2.0000000000000004*Power(t0,6)*te*xe + 12.*Power(t0,5)*t1*te*xe - 90.*Power(t0,4)*Power(t1,2)*te*xe + 
           200.*Power(t0,3)*Power(t1,3)*te*xe - 210.*Power(t0,2)*Power(t1,4)*te*xe + 108.*t0*Power(t1,5)*te*xe - 22.*Power(t1,6)*te*xe - 
           12.*Power(t0,5)*Power(te,2)*xe + 60.*Power(t0,4)*t1*Power(te,2)*xe - 120.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 
           120.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 60.*t0*Power(t1,4)*Power(te,2)*xe + 12.*Power(t1,5)*Power(te,2)*xe + 2.0000000000000004*Power(t0,6)*t1*xred - 
           18.*Power(t0,4)*Power(t1,3)*xred + 32.*Power(t0,3)*Power(t1,4)*xred - 18.*Power(t0,2)*Power(t1,5)*xred + 2.*Power(t1,7)*xred - 
           2.0000000000000004*Power(t0,6)*te*xred - 12.*Power(t0,5)*t1*te*xred + 54.00000000000001*Power(t0,4)*Power(t1,2)*te*xred - 
           56.*Power(t0,3)*Power(t1,3)*te*xred - 5.9999999999999964*Power(t0,2)*Power(t1,4)*te*xred + 36.*t0*Power(t1,5)*te*xred - 14.*Power(t1,6)*te*xred + 
           12.*Power(t0,5)*Power(te,2)*xred - 24.*Power(t0,4)*t1*Power(te,2)*xred - 24.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
           96.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 84.00000000000001*t0*Power(t1,4)*Power(te,2)*xred + 24.*Power(t1,5)*Power(te,2)*xred - 
           12.*Power(t0,4)*Power(te,3)*xred + 48.*Power(t0,3)*t1*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 
           48.00000000000001*t0*Power(t1,3)*Power(te,3)*xred - 12.000000000000004*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.05555555555555555*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 0.16666666666666666*Power(t0,2)*Power(t1,4) - 
           0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 
           0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 0.08333333333333334*Power(t0,4)*Power(te,2) - 
           0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 0.25*Power(t1,4)*Power(te,2) - 
           0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 
           0.11111111111111113*Power(t1,3)*Power(te,3))*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 
           10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 
           30.*Power(t0,3)*Power(t1,4)*te*v0 + 30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) + 
      (0.08333333333333333*(-t1 + te)*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 
           0.16666666666666666*Power(t0,2)*Power(t1,4) - 0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 
           0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 
           0.08333333333333334*Power(t0,4)*Power(te,2) - 0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 
           0.25*Power(t1,4)*Power(te,2) - 0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 
           0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 0.11111111111111113*Power(t1,3)*Power(te,3))*
         (-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 
           1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
           40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
           30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
           30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
           40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
           4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
           10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) - 
      (0.08333333333333333*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 
           5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
           40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
           30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
           30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
           40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
           4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
           10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.16666666666666666*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,3)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))))*
    ((0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 10.*Power(t0,2)*Power(t1,6)*v0 + 
           5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 30.*Power(t0,3)*Power(t1,4)*te*v0 + 
           30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.08333333333333333*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)))) + 
   ((-0.05555555555555555*(-3.*Power(t0,5)*Power(t1,2)*v0 + 15.*Power(t0,4)*Power(t1,3)*v0 - 30.*Power(t0,3)*Power(t1,4)*v0 + 30.*Power(t0,2)*Power(t1,5)*v0 - 
           15.*t0*Power(t1,6)*v0 + 3.*Power(t1,7)*v0 + 6.*Power(t0,5)*t1*te*v0 - 30.*Power(t0,4)*Power(t1,2)*te*v0 + 60.*Power(t0,3)*Power(t1,3)*te*v0 - 
           60.*Power(t0,2)*Power(t1,4)*te*v0 + 30.000000000000004*t0*Power(t1,5)*te*v0 - 6.000000000000001*Power(t1,6)*te*v0 - 3.*Power(t0,5)*Power(te,2)*v0 + 
           15.*Power(t0,4)*t1*Power(te,2)*v0 - 30.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*v0 + 30.000000000000007*Power(t0,2)*Power(t1,3)*Power(te,2)*v0 - 
           15.000000000000005*t0*Power(t1,4)*Power(te,2)*v0 + 3.0000000000000013*Power(t1,5)*Power(te,2)*v0 + 3.*Power(t0,5)*Power(t1,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*ve + 30.*Power(t0,3)*Power(t1,4)*ve - 30.*Power(t0,2)*Power(t1,5)*ve + 15.*t0*Power(t1,6)*ve - 3.*Power(t1,7)*ve - 
           6.*Power(t0,5)*t1*te*ve + 30.*Power(t0,4)*Power(t1,2)*te*ve - 60.*Power(t0,3)*Power(t1,3)*te*ve + 60.*Power(t0,2)*Power(t1,4)*te*ve - 
           30.*t0*Power(t1,5)*te*ve + 6.*Power(t1,6)*te*ve + 3.*Power(t0,5)*Power(te,2)*ve - 15.*Power(t0,4)*t1*Power(te,2)*ve + 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*ve - 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*ve + 15.*t0*Power(t1,4)*Power(te,2)*ve - 3.*Power(t1,5)*Power(te,2)*ve + 
           9.*Power(t0,4)*Power(t1,2)*x0 - 36.*Power(t0,3)*Power(t1,3)*x0 + 54.*Power(t0,2)*Power(t1,4)*x0 - 36.*t0*Power(t1,5)*x0 + 9.*Power(t1,6)*x0 - 
           18.*Power(t0,4)*t1*te*x0 + 72.*Power(t0,3)*Power(t1,2)*te*x0 - 108.*Power(t0,2)*Power(t1,3)*te*x0 + 72.*t0*Power(t1,4)*te*x0 - 18.*Power(t1,5)*te*x0 + 
           9.*Power(t0,4)*Power(te,2)*x0 - 36.*Power(t0,3)*t1*Power(te,2)*x0 + 54.*Power(t0,2)*Power(t1,2)*Power(te,2)*x0 - 
           36.00000000000001*t0*Power(t1,3)*Power(te,2)*x0 + 9.000000000000004*Power(t1,4)*Power(te,2)*x0 + 6.*Power(t0,5)*t1*xe - 30.*Power(t0,4)*Power(t1,2)*xe + 
           60.*Power(t0,3)*Power(t1,3)*xe - 60.*Power(t0,2)*Power(t1,4)*xe + 30.*t0*Power(t1,5)*xe - 6.*Power(t1,6)*xe - 6.*Power(t0,5)*te*xe + 
           30.*Power(t0,4)*t1*te*xe - 60.*Power(t0,3)*Power(t1,2)*te*xe + 60.*Power(t0,2)*Power(t1,3)*te*xe - 30.*t0*Power(t1,4)*te*xe + 6.*Power(t1,5)*te*xe - 
           6.*Power(t0,5)*t1*xred + 21.*Power(t0,4)*Power(t1,2)*xred - 24.*Power(t0,3)*Power(t1,3)*xred + 6.*Power(t0,2)*Power(t1,4)*xred + 6.*t0*Power(t1,5)*xred - 
           3.*Power(t1,6)*xred + 6.*Power(t0,5)*te*xred - 12.*Power(t0,4)*t1*te*xred - 12.000000000000002*Power(t0,3)*Power(t1,2)*te*xred + 
           48.*Power(t0,2)*Power(t1,3)*te*xred - 42.*t0*Power(t1,4)*te*xred + 12.*Power(t1,5)*te*xred - 9.*Power(t0,4)*Power(te,2)*xred + 
           36.*Power(t0,3)*t1*Power(te,2)*xred - 54.*Power(t0,2)*Power(t1,2)*Power(te,2)*xred + 36.00000000000001*t0*Power(t1,3)*Power(te,2)*xred - 
           9.000000000000004*Power(t1,4)*Power(te,2)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
      (0.08333333333333333*(-t1 + te)*(4.*Power(t0,5)*Power(t1,3)*v0 - 20.*Power(t0,4)*Power(t1,4)*v0 + 40.*Power(t0,3)*Power(t1,5)*v0 - 
           40.*Power(t0,2)*Power(t1,6)*v0 + 20.*t0*Power(t1,7)*v0 - 4.*Power(t1,8)*v0 - 12.*Power(t0,5)*Power(t1,2)*te*v0 + 60.*Power(t0,4)*Power(t1,3)*te*v0 - 
           120.*Power(t0,3)*Power(t1,4)*te*v0 + 120.*Power(t0,2)*Power(t1,5)*te*v0 - 60.*t0*Power(t1,6)*te*v0 + 12.*Power(t1,7)*te*v0 + 
           12.*Power(t0,5)*t1*Power(te,2)*v0 - 60.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 120.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
           120.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 60.000000000000014*t0*Power(t1,5)*Power(te,2)*v0 - 12.000000000000004*Power(t1,6)*Power(te,2)*v0 - 
           4.*Power(t0,5)*Power(te,3)*v0 + 20.*Power(t0,4)*t1*Power(te,3)*v0 - 40.00000000000001*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
           40.00000000000001*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 20.000000000000007*t0*Power(t1,4)*Power(te,3)*v0 + 4.000000000000003*Power(t1,5)*Power(te,3)*v0 - 
           3.0000000000000004*Power(t0,6)*Power(t1,2)*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*ve - 5.000000000000001*Power(t0,4)*Power(t1,4)*ve - 
           20.*Power(t0,3)*Power(t1,5)*ve + 35.*Power(t0,2)*Power(t1,6)*ve - 22.*t0*Power(t1,7)*ve + 5.*Power(t1,8)*ve + 6.000000000000001*Power(t0,6)*t1*te*ve - 
           12.000000000000002*Power(t0,5)*Power(t1,2)*te*ve - 30.*Power(t0,4)*Power(t1,3)*te*ve + 120.*Power(t0,3)*Power(t1,4)*te*ve - 
           150.*Power(t0,2)*Power(t1,5)*te*ve + 84.*t0*Power(t1,6)*te*ve - 18.*Power(t1,7)*te*ve - 3.000000000000001*Power(t0,6)*Power(te,2)*ve - 
           5.999999999999999*Power(t0,5)*t1*Power(te,2)*ve + 75.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 180.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 
           195.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 102.*t0*Power(t1,5)*Power(te,2)*ve + 21.000000000000004*Power(t1,6)*Power(te,2)*ve + 
           8.*Power(t0,5)*Power(te,3)*ve - 40.*Power(t0,4)*t1*Power(te,3)*ve + 80.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 80.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           40.*t0*Power(t1,4)*Power(te,3)*ve - 8.*Power(t1,5)*Power(te,3)*ve - 12.*Power(t0,4)*Power(t1,3)*x0 + 48.*Power(t0,3)*Power(t1,4)*x0 - 
           72.*Power(t0,2)*Power(t1,5)*x0 + 48.*t0*Power(t1,6)*x0 - 12.*Power(t1,7)*x0 + 36.*Power(t0,4)*Power(t1,2)*te*x0 - 144.*Power(t0,3)*Power(t1,3)*te*x0 + 
           216.*Power(t0,2)*Power(t1,4)*te*x0 - 144.*t0*Power(t1,5)*te*x0 + 36.*Power(t1,6)*te*x0 - 36.00000000000001*Power(t0,4)*t1*Power(te,2)*x0 + 
           144.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 216.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 144.*t0*Power(t1,4)*Power(te,2)*x0 - 
           36.00000000000001*Power(t1,5)*Power(te,2)*x0 + 12.*Power(t0,4)*Power(te,3)*x0 - 48.*Power(t0,3)*t1*Power(te,3)*x0 + 
           72.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,3)*Power(te,3)*x0 + 12.000000000000004*Power(t1,4)*Power(te,3)*x0 - 
           2.0000000000000004*Power(t0,6)*t1*xe + 30.*Power(t0,4)*Power(t1,3)*xe - 80.*Power(t0,3)*Power(t1,4)*xe + 90.*Power(t0,2)*Power(t1,5)*xe - 
           48.*t0*Power(t1,6)*xe + 10.*Power(t1,7)*xe + 2.0000000000000004*Power(t0,6)*te*xe + 12.*Power(t0,5)*t1*te*xe - 90.*Power(t0,4)*Power(t1,2)*te*xe + 
           200.*Power(t0,3)*Power(t1,3)*te*xe - 210.*Power(t0,2)*Power(t1,4)*te*xe + 108.*t0*Power(t1,5)*te*xe - 22.*Power(t1,6)*te*xe - 
           12.*Power(t0,5)*Power(te,2)*xe + 60.*Power(t0,4)*t1*Power(te,2)*xe - 120.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 
           120.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 60.*t0*Power(t1,4)*Power(te,2)*xe + 12.*Power(t1,5)*Power(te,2)*xe + 2.0000000000000004*Power(t0,6)*t1*xred - 
           18.*Power(t0,4)*Power(t1,3)*xred + 32.*Power(t0,3)*Power(t1,4)*xred - 18.*Power(t0,2)*Power(t1,5)*xred + 2.*Power(t1,7)*xred - 
           2.0000000000000004*Power(t0,6)*te*xred - 12.*Power(t0,5)*t1*te*xred + 54.00000000000001*Power(t0,4)*Power(t1,2)*te*xred - 
           56.*Power(t0,3)*Power(t1,3)*te*xred - 5.9999999999999964*Power(t0,2)*Power(t1,4)*te*xred + 36.*t0*Power(t1,5)*te*xred - 14.*Power(t1,6)*te*xred + 
           12.*Power(t0,5)*Power(te,2)*xred - 24.*Power(t0,4)*t1*Power(te,2)*xred - 24.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
           96.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 84.00000000000001*t0*Power(t1,4)*Power(te,2)*xred + 24.*Power(t1,5)*Power(te,2)*xred - 
           12.*Power(t0,4)*Power(te,3)*xred + 48.*Power(t0,3)*t1*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 
           48.00000000000001*t0*Power(t1,3)*Power(te,3)*xred - 12.000000000000004*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
      (0.05555555555555555*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 0.16666666666666666*Power(t0,2)*Power(t1,4) - 
           0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 
           0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 0.08333333333333334*Power(t0,4)*Power(te,2) - 
           0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 0.25*Power(t1,4)*Power(te,2) - 
           0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 
           0.11111111111111113*Power(t1,3)*Power(te,3))*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 
           10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 
           30.*Power(t0,3)*Power(t1,4)*te*v0 + 30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) - 
      (0.08333333333333333*(-t1 + te)*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 
           0.16666666666666666*Power(t0,2)*Power(t1,4) - 0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 
           0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 
           0.08333333333333334*Power(t0,4)*Power(te,2) - 0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 
           0.25*Power(t1,4)*Power(te,2) - 0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 
           0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 0.11111111111111113*Power(t1,3)*Power(te,3))*
         (-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 
           1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
           40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
           30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
           30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
           40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
           4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
           10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) + 
      (0.08333333333333333*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 
           5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
           40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
           30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
           30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
           40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
           4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
           10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
      (0.16666666666666666*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,3)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))))*
    ((0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 10.*Power(t0,2)*Power(t1,6)*v0 + 
           5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 30.*Power(t0,3)*Power(t1,4)*te*v0 + 
           30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.08333333333333333*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)))) + 
   ((0.05555555555555555*(-3.*Power(t0,5)*Power(t1,2)*v0 + 15.*Power(t0,4)*Power(t1,3)*v0 - 30.*Power(t0,3)*Power(t1,4)*v0 + 30.*Power(t0,2)*Power(t1,5)*v0 - 
           15.*t0*Power(t1,6)*v0 + 3.*Power(t1,7)*v0 + 6.*Power(t0,5)*t1*te*v0 - 30.*Power(t0,4)*Power(t1,2)*te*v0 + 60.*Power(t0,3)*Power(t1,3)*te*v0 - 
           60.*Power(t0,2)*Power(t1,4)*te*v0 + 30.000000000000004*t0*Power(t1,5)*te*v0 - 6.000000000000001*Power(t1,6)*te*v0 - 3.*Power(t0,5)*Power(te,2)*v0 + 
           15.*Power(t0,4)*t1*Power(te,2)*v0 - 30.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*v0 + 30.000000000000007*Power(t0,2)*Power(t1,3)*Power(te,2)*v0 - 
           15.000000000000005*t0*Power(t1,4)*Power(te,2)*v0 + 3.0000000000000013*Power(t1,5)*Power(te,2)*v0 + 3.*Power(t0,5)*Power(t1,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*ve + 30.*Power(t0,3)*Power(t1,4)*ve - 30.*Power(t0,2)*Power(t1,5)*ve + 15.*t0*Power(t1,6)*ve - 3.*Power(t1,7)*ve - 
           6.*Power(t0,5)*t1*te*ve + 30.*Power(t0,4)*Power(t1,2)*te*ve - 60.*Power(t0,3)*Power(t1,3)*te*ve + 60.*Power(t0,2)*Power(t1,4)*te*ve - 
           30.*t0*Power(t1,5)*te*ve + 6.*Power(t1,6)*te*ve + 3.*Power(t0,5)*Power(te,2)*ve - 15.*Power(t0,4)*t1*Power(te,2)*ve + 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*ve - 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*ve + 15.*t0*Power(t1,4)*Power(te,2)*ve - 3.*Power(t1,5)*Power(te,2)*ve + 
           9.*Power(t0,4)*Power(t1,2)*x0 - 36.*Power(t0,3)*Power(t1,3)*x0 + 54.*Power(t0,2)*Power(t1,4)*x0 - 36.*t0*Power(t1,5)*x0 + 9.*Power(t1,6)*x0 - 
           18.*Power(t0,4)*t1*te*x0 + 72.*Power(t0,3)*Power(t1,2)*te*x0 - 108.*Power(t0,2)*Power(t1,3)*te*x0 + 72.*t0*Power(t1,4)*te*x0 - 18.*Power(t1,5)*te*x0 + 
           9.*Power(t0,4)*Power(te,2)*x0 - 36.*Power(t0,3)*t1*Power(te,2)*x0 + 54.*Power(t0,2)*Power(t1,2)*Power(te,2)*x0 - 
           36.00000000000001*t0*Power(t1,3)*Power(te,2)*x0 + 9.000000000000004*Power(t1,4)*Power(te,2)*x0 + 6.*Power(t0,5)*t1*xe - 30.*Power(t0,4)*Power(t1,2)*xe + 
           60.*Power(t0,3)*Power(t1,3)*xe - 60.*Power(t0,2)*Power(t1,4)*xe + 30.*t0*Power(t1,5)*xe - 6.*Power(t1,6)*xe - 6.*Power(t0,5)*te*xe + 
           30.*Power(t0,4)*t1*te*xe - 60.*Power(t0,3)*Power(t1,2)*te*xe + 60.*Power(t0,2)*Power(t1,3)*te*xe - 30.*t0*Power(t1,4)*te*xe + 6.*Power(t1,5)*te*xe - 
           6.*Power(t0,5)*t1*xred + 21.*Power(t0,4)*Power(t1,2)*xred - 24.*Power(t0,3)*Power(t1,3)*xred + 6.*Power(t0,2)*Power(t1,4)*xred + 6.*t0*Power(t1,5)*xred - 
           3.*Power(t1,6)*xred + 6.*Power(t0,5)*te*xred - 12.*Power(t0,4)*t1*te*xred - 12.000000000000002*Power(t0,3)*Power(t1,2)*te*xred + 
           48.*Power(t0,2)*Power(t1,3)*te*xred - 42.*t0*Power(t1,4)*te*xred + 12.*Power(t1,5)*te*xred - 9.*Power(t0,4)*Power(te,2)*xred + 
           36.*Power(t0,3)*t1*Power(te,2)*xred - 54.*Power(t0,2)*Power(t1,2)*Power(te,2)*xred + 36.00000000000001*t0*Power(t1,3)*Power(te,2)*xred - 
           9.000000000000004*Power(t1,4)*Power(te,2)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.08333333333333333*(-t1 + te)*(4.*Power(t0,5)*Power(t1,3)*v0 - 20.*Power(t0,4)*Power(t1,4)*v0 + 40.*Power(t0,3)*Power(t1,5)*v0 - 
           40.*Power(t0,2)*Power(t1,6)*v0 + 20.*t0*Power(t1,7)*v0 - 4.*Power(t1,8)*v0 - 12.*Power(t0,5)*Power(t1,2)*te*v0 + 60.*Power(t0,4)*Power(t1,3)*te*v0 - 
           120.*Power(t0,3)*Power(t1,4)*te*v0 + 120.*Power(t0,2)*Power(t1,5)*te*v0 - 60.*t0*Power(t1,6)*te*v0 + 12.*Power(t1,7)*te*v0 + 
           12.*Power(t0,5)*t1*Power(te,2)*v0 - 60.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 120.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
           120.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 60.000000000000014*t0*Power(t1,5)*Power(te,2)*v0 - 12.000000000000004*Power(t1,6)*Power(te,2)*v0 - 
           4.*Power(t0,5)*Power(te,3)*v0 + 20.*Power(t0,4)*t1*Power(te,3)*v0 - 40.00000000000001*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
           40.00000000000001*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 20.000000000000007*t0*Power(t1,4)*Power(te,3)*v0 + 4.000000000000003*Power(t1,5)*Power(te,3)*v0 - 
           3.0000000000000004*Power(t0,6)*Power(t1,2)*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*ve - 5.000000000000001*Power(t0,4)*Power(t1,4)*ve - 
           20.*Power(t0,3)*Power(t1,5)*ve + 35.*Power(t0,2)*Power(t1,6)*ve - 22.*t0*Power(t1,7)*ve + 5.*Power(t1,8)*ve + 6.000000000000001*Power(t0,6)*t1*te*ve - 
           12.000000000000002*Power(t0,5)*Power(t1,2)*te*ve - 30.*Power(t0,4)*Power(t1,3)*te*ve + 120.*Power(t0,3)*Power(t1,4)*te*ve - 
           150.*Power(t0,2)*Power(t1,5)*te*ve + 84.*t0*Power(t1,6)*te*ve - 18.*Power(t1,7)*te*ve - 3.000000000000001*Power(t0,6)*Power(te,2)*ve - 
           5.999999999999999*Power(t0,5)*t1*Power(te,2)*ve + 75.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 180.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 
           195.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 102.*t0*Power(t1,5)*Power(te,2)*ve + 21.000000000000004*Power(t1,6)*Power(te,2)*ve + 
           8.*Power(t0,5)*Power(te,3)*ve - 40.*Power(t0,4)*t1*Power(te,3)*ve + 80.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 80.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           40.*t0*Power(t1,4)*Power(te,3)*ve - 8.*Power(t1,5)*Power(te,3)*ve - 12.*Power(t0,4)*Power(t1,3)*x0 + 48.*Power(t0,3)*Power(t1,4)*x0 - 
           72.*Power(t0,2)*Power(t1,5)*x0 + 48.*t0*Power(t1,6)*x0 - 12.*Power(t1,7)*x0 + 36.*Power(t0,4)*Power(t1,2)*te*x0 - 144.*Power(t0,3)*Power(t1,3)*te*x0 + 
           216.*Power(t0,2)*Power(t1,4)*te*x0 - 144.*t0*Power(t1,5)*te*x0 + 36.*Power(t1,6)*te*x0 - 36.00000000000001*Power(t0,4)*t1*Power(te,2)*x0 + 
           144.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 216.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 144.*t0*Power(t1,4)*Power(te,2)*x0 - 
           36.00000000000001*Power(t1,5)*Power(te,2)*x0 + 12.*Power(t0,4)*Power(te,3)*x0 - 48.*Power(t0,3)*t1*Power(te,3)*x0 + 
           72.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,3)*Power(te,3)*x0 + 12.000000000000004*Power(t1,4)*Power(te,3)*x0 - 
           2.0000000000000004*Power(t0,6)*t1*xe + 30.*Power(t0,4)*Power(t1,3)*xe - 80.*Power(t0,3)*Power(t1,4)*xe + 90.*Power(t0,2)*Power(t1,5)*xe - 
           48.*t0*Power(t1,6)*xe + 10.*Power(t1,7)*xe + 2.0000000000000004*Power(t0,6)*te*xe + 12.*Power(t0,5)*t1*te*xe - 90.*Power(t0,4)*Power(t1,2)*te*xe + 
           200.*Power(t0,3)*Power(t1,3)*te*xe - 210.*Power(t0,2)*Power(t1,4)*te*xe + 108.*t0*Power(t1,5)*te*xe - 22.*Power(t1,6)*te*xe - 
           12.*Power(t0,5)*Power(te,2)*xe + 60.*Power(t0,4)*t1*Power(te,2)*xe - 120.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 
           120.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 60.*t0*Power(t1,4)*Power(te,2)*xe + 12.*Power(t1,5)*Power(te,2)*xe + 2.0000000000000004*Power(t0,6)*t1*xred - 
           18.*Power(t0,4)*Power(t1,3)*xred + 32.*Power(t0,3)*Power(t1,4)*xred - 18.*Power(t0,2)*Power(t1,5)*xred + 2.*Power(t1,7)*xred - 
           2.0000000000000004*Power(t0,6)*te*xred - 12.*Power(t0,5)*t1*te*xred + 54.00000000000001*Power(t0,4)*Power(t1,2)*te*xred - 
           56.*Power(t0,3)*Power(t1,3)*te*xred - 5.9999999999999964*Power(t0,2)*Power(t1,4)*te*xred + 36.*t0*Power(t1,5)*te*xred - 14.*Power(t1,6)*te*xred + 
           12.*Power(t0,5)*Power(te,2)*xred - 24.*Power(t0,4)*t1*Power(te,2)*xred - 24.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
           96.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 84.00000000000001*t0*Power(t1,4)*Power(te,2)*xred + 24.*Power(t1,5)*Power(te,2)*xred - 
           12.*Power(t0,4)*Power(te,3)*xred + 48.*Power(t0,3)*t1*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 
           48.00000000000001*t0*Power(t1,3)*Power(te,3)*xred - 12.000000000000004*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.05555555555555555*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 0.16666666666666666*Power(t0,2)*Power(t1,4) - 
           0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 
           0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 0.08333333333333334*Power(t0,4)*Power(te,2) - 
           0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 0.25*Power(t1,4)*Power(te,2) - 
           0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 
           0.11111111111111113*Power(t1,3)*Power(te,3))*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 
           10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 
           30.*Power(t0,3)*Power(t1,4)*te*v0 + 30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) + 
      (0.08333333333333333*(-t1 + te)*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 
           0.16666666666666666*Power(t0,2)*Power(t1,4) - 0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 
           0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 
           0.08333333333333334*Power(t0,4)*Power(te,2) - 0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 
           0.25*Power(t1,4)*Power(te,2) - 0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 
           0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 0.11111111111111113*Power(t1,3)*Power(te,3))*
         (-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 
           1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
           40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
           30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
           30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
           40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
           4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
           10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) - 
      (0.08333333333333333*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 
           5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
           40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
           30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
           30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
           40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
           4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
           10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
      (0.16666666666666666*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,3)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))))*
    ((-0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 10.*Power(t0,2)*Power(t1,6)*v0 + 
           5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 30.*Power(t0,3)*Power(t1,4)*te*v0 + 
           30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
           15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
           15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
           5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
           5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
           3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
           15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
           30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
           1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
           5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
           18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
           54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
           36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
           3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
           15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
           6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
           30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
           30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 
           3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 
           3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
           6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 
           6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
           21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 12.*Power(t0,3)*t1*Power(te,3)*xred - 
           18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
           0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
           0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
           0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
           0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
           0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
           0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
      (0.08333333333333333*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
           10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
           40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
           6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
           60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
           20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
           20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
           5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
           5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
           4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
           1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
           5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
           5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
           15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
           42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
           1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
           65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
           2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
           10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
           18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
           72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
           72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 
           18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 
           72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 
           3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 
           14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 
           90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 
           6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
           105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 
           20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
           20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
           2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
           2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
           18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
           27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
           2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
           4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
           32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
           3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
           12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
       ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
         (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
           0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
           0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
           0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
           0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
           0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
           0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)))) - 
   (0.08333333333333333*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 
        5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
        40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
        30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
        30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
        40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
        4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
        10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
        5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
        4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
        1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
        5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
        5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
        15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
        42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
        1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
        65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 2.*Power(t0,5)*Power(te,4)*ve - 
        10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
        10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
        18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
        72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
        72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 
        12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 
        48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
        18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
        1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 
        25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 
        80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 
        1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
        100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 
        4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
        20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
        2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
        2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
        18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
        27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
        2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
        4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
        32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
        3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
        12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred)*
      ((0.05555555555555555*(-t1 + te)*(-3.*Power(t0,5)*Power(t1,2)*v0 + 15.*Power(t0,4)*Power(t1,3)*v0 - 30.*Power(t0,3)*Power(t1,4)*v0 + 
             30.*Power(t0,2)*Power(t1,5)*v0 - 15.*t0*Power(t1,6)*v0 + 3.*Power(t1,7)*v0 + 6.*Power(t0,5)*t1*te*v0 - 30.*Power(t0,4)*Power(t1,2)*te*v0 + 
             60.*Power(t0,3)*Power(t1,3)*te*v0 - 60.*Power(t0,2)*Power(t1,4)*te*v0 + 30.000000000000004*t0*Power(t1,5)*te*v0 - 6.000000000000001*Power(t1,6)*te*v0 - 
             3.*Power(t0,5)*Power(te,2)*v0 + 15.*Power(t0,4)*t1*Power(te,2)*v0 - 30.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*v0 + 
             30.000000000000007*Power(t0,2)*Power(t1,3)*Power(te,2)*v0 - 15.000000000000005*t0*Power(t1,4)*Power(te,2)*v0 + 
             3.0000000000000013*Power(t1,5)*Power(te,2)*v0 + 3.*Power(t0,5)*Power(t1,2)*ve - 15.*Power(t0,4)*Power(t1,3)*ve + 30.*Power(t0,3)*Power(t1,4)*ve - 
             30.*Power(t0,2)*Power(t1,5)*ve + 15.*t0*Power(t1,6)*ve - 3.*Power(t1,7)*ve - 6.*Power(t0,5)*t1*te*ve + 30.*Power(t0,4)*Power(t1,2)*te*ve - 
             60.*Power(t0,3)*Power(t1,3)*te*ve + 60.*Power(t0,2)*Power(t1,4)*te*ve - 30.*t0*Power(t1,5)*te*ve + 6.*Power(t1,6)*te*ve + 3.*Power(t0,5)*Power(te,2)*ve - 
             15.*Power(t0,4)*t1*Power(te,2)*ve + 30.*Power(t0,3)*Power(t1,2)*Power(te,2)*ve - 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*ve + 
             15.*t0*Power(t1,4)*Power(te,2)*ve - 3.*Power(t1,5)*Power(te,2)*ve + 9.*Power(t0,4)*Power(t1,2)*x0 - 36.*Power(t0,3)*Power(t1,3)*x0 + 
             54.*Power(t0,2)*Power(t1,4)*x0 - 36.*t0*Power(t1,5)*x0 + 9.*Power(t1,6)*x0 - 18.*Power(t0,4)*t1*te*x0 + 72.*Power(t0,3)*Power(t1,2)*te*x0 - 
             108.*Power(t0,2)*Power(t1,3)*te*x0 + 72.*t0*Power(t1,4)*te*x0 - 18.*Power(t1,5)*te*x0 + 9.*Power(t0,4)*Power(te,2)*x0 - 36.*Power(t0,3)*t1*Power(te,2)*x0 + 
             54.*Power(t0,2)*Power(t1,2)*Power(te,2)*x0 - 36.00000000000001*t0*Power(t1,3)*Power(te,2)*x0 + 9.000000000000004*Power(t1,4)*Power(te,2)*x0 + 
             6.*Power(t0,5)*t1*xe - 30.*Power(t0,4)*Power(t1,2)*xe + 60.*Power(t0,3)*Power(t1,3)*xe - 60.*Power(t0,2)*Power(t1,4)*xe + 30.*t0*Power(t1,5)*xe - 
             6.*Power(t1,6)*xe - 6.*Power(t0,5)*te*xe + 30.*Power(t0,4)*t1*te*xe - 60.*Power(t0,3)*Power(t1,2)*te*xe + 60.*Power(t0,2)*Power(t1,3)*te*xe - 
             30.*t0*Power(t1,4)*te*xe + 6.*Power(t1,5)*te*xe - 6.*Power(t0,5)*t1*xred + 21.*Power(t0,4)*Power(t1,2)*xred - 24.*Power(t0,3)*Power(t1,3)*xred + 
             6.*Power(t0,2)*Power(t1,4)*xred + 6.*t0*Power(t1,5)*xred - 3.*Power(t1,6)*xred + 6.*Power(t0,5)*te*xred - 12.*Power(t0,4)*t1*te*xred - 
             12.000000000000002*Power(t0,3)*Power(t1,2)*te*xred + 48.*Power(t0,2)*Power(t1,3)*te*xred - 42.*t0*Power(t1,4)*te*xred + 12.*Power(t1,5)*te*xred - 
             9.*Power(t0,4)*Power(te,2)*xred + 36.*Power(t0,3)*t1*Power(te,2)*xred - 54.*Power(t0,2)*Power(t1,2)*Power(te,2)*xred + 
             36.00000000000001*t0*Power(t1,3)*Power(te,2)*xred - 9.000000000000004*Power(t1,4)*Power(te,2)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
        (0.041666666666666664*Power(-t1 + te,2)*(4.*Power(t0,5)*Power(t1,3)*v0 - 20.*Power(t0,4)*Power(t1,4)*v0 + 40.*Power(t0,3)*Power(t1,5)*v0 - 
             40.*Power(t0,2)*Power(t1,6)*v0 + 20.*t0*Power(t1,7)*v0 - 4.*Power(t1,8)*v0 - 12.*Power(t0,5)*Power(t1,2)*te*v0 + 60.*Power(t0,4)*Power(t1,3)*te*v0 - 
             120.*Power(t0,3)*Power(t1,4)*te*v0 + 120.*Power(t0,2)*Power(t1,5)*te*v0 - 60.*t0*Power(t1,6)*te*v0 + 12.*Power(t1,7)*te*v0 + 
             12.*Power(t0,5)*t1*Power(te,2)*v0 - 60.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 120.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
             120.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 60.000000000000014*t0*Power(t1,5)*Power(te,2)*v0 - 12.000000000000004*Power(t1,6)*Power(te,2)*v0 - 
             4.*Power(t0,5)*Power(te,3)*v0 + 20.*Power(t0,4)*t1*Power(te,3)*v0 - 40.00000000000001*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
             40.00000000000001*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 20.000000000000007*t0*Power(t1,4)*Power(te,3)*v0 + 
             4.000000000000003*Power(t1,5)*Power(te,3)*v0 - 3.0000000000000004*Power(t0,6)*Power(t1,2)*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*ve - 
             5.000000000000001*Power(t0,4)*Power(t1,4)*ve - 20.*Power(t0,3)*Power(t1,5)*ve + 35.*Power(t0,2)*Power(t1,6)*ve - 22.*t0*Power(t1,7)*ve + 
             5.*Power(t1,8)*ve + 6.000000000000001*Power(t0,6)*t1*te*ve - 12.000000000000002*Power(t0,5)*Power(t1,2)*te*ve - 30.*Power(t0,4)*Power(t1,3)*te*ve + 
             120.*Power(t0,3)*Power(t1,4)*te*ve - 150.*Power(t0,2)*Power(t1,5)*te*ve + 84.*t0*Power(t1,6)*te*ve - 18.*Power(t1,7)*te*ve - 
             3.000000000000001*Power(t0,6)*Power(te,2)*ve - 5.999999999999999*Power(t0,5)*t1*Power(te,2)*ve + 75.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
             180.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 195.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 102.*t0*Power(t1,5)*Power(te,2)*ve + 
             21.000000000000004*Power(t1,6)*Power(te,2)*ve + 8.*Power(t0,5)*Power(te,3)*ve - 40.*Power(t0,4)*t1*Power(te,3)*ve + 
             80.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 80.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 40.*t0*Power(t1,4)*Power(te,3)*ve - 
             8.*Power(t1,5)*Power(te,3)*ve - 12.*Power(t0,4)*Power(t1,3)*x0 + 48.*Power(t0,3)*Power(t1,4)*x0 - 72.*Power(t0,2)*Power(t1,5)*x0 + 48.*t0*Power(t1,6)*x0 - 
             12.*Power(t1,7)*x0 + 36.*Power(t0,4)*Power(t1,2)*te*x0 - 144.*Power(t0,3)*Power(t1,3)*te*x0 + 216.*Power(t0,2)*Power(t1,4)*te*x0 - 
             144.*t0*Power(t1,5)*te*x0 + 36.*Power(t1,6)*te*x0 - 36.00000000000001*Power(t0,4)*t1*Power(te,2)*x0 + 144.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 
             216.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 144.*t0*Power(t1,4)*Power(te,2)*x0 - 36.00000000000001*Power(t1,5)*Power(te,2)*x0 + 
             12.*Power(t0,4)*Power(te,3)*x0 - 48.*Power(t0,3)*t1*Power(te,3)*x0 + 72.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
             48.00000000000001*t0*Power(t1,3)*Power(te,3)*x0 + 12.000000000000004*Power(t1,4)*Power(te,3)*x0 - 2.0000000000000004*Power(t0,6)*t1*xe + 
             30.*Power(t0,4)*Power(t1,3)*xe - 80.*Power(t0,3)*Power(t1,4)*xe + 90.*Power(t0,2)*Power(t1,5)*xe - 48.*t0*Power(t1,6)*xe + 10.*Power(t1,7)*xe + 
             2.0000000000000004*Power(t0,6)*te*xe + 12.*Power(t0,5)*t1*te*xe - 90.*Power(t0,4)*Power(t1,2)*te*xe + 200.*Power(t0,3)*Power(t1,3)*te*xe - 
             210.*Power(t0,2)*Power(t1,4)*te*xe + 108.*t0*Power(t1,5)*te*xe - 22.*Power(t1,6)*te*xe - 12.*Power(t0,5)*Power(te,2)*xe + 
             60.*Power(t0,4)*t1*Power(te,2)*xe - 120.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 120.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 
             60.*t0*Power(t1,4)*Power(te,2)*xe + 12.*Power(t1,5)*Power(te,2)*xe + 2.0000000000000004*Power(t0,6)*t1*xred - 18.*Power(t0,4)*Power(t1,3)*xred + 
             32.*Power(t0,3)*Power(t1,4)*xred - 18.*Power(t0,2)*Power(t1,5)*xred + 2.*Power(t1,7)*xred - 2.0000000000000004*Power(t0,6)*te*xred - 
             12.*Power(t0,5)*t1*te*xred + 54.00000000000001*Power(t0,4)*Power(t1,2)*te*xred - 56.*Power(t0,3)*Power(t1,3)*te*xred - 
             5.9999999999999964*Power(t0,2)*Power(t1,4)*te*xred + 36.*t0*Power(t1,5)*te*xred - 14.*Power(t1,6)*te*xred + 12.*Power(t0,5)*Power(te,2)*xred - 
             24.*Power(t0,4)*t1*Power(te,2)*xred - 24.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 96.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
             84.00000000000001*t0*Power(t1,4)*Power(te,2)*xred + 24.*Power(t1,5)*Power(te,2)*xred - 12.*Power(t0,4)*Power(te,3)*xred + 
             48.*Power(t0,3)*t1*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 48.00000000000001*t0*Power(t1,3)*Power(te,3)*xred - 
             12.000000000000004*Power(t1,4)*Power(te,3)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
           (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
             0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
             0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
             0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
             0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
             0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
             0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
        (0.05555555555555555*(-t1 + te)*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 
             0.16666666666666666*Power(t0,2)*Power(t1,4) - 0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 
             0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 
             0.08333333333333334*Power(t0,4)*Power(te,2) - 0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 
             0.25*Power(t1,4)*Power(te,2) - 0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 
             0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 0.11111111111111113*Power(t1,3)*Power(te,3))*
           (1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 
             1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 30.*Power(t0,3)*Power(t1,4)*te*v0 + 
             30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
             15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
             15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
             5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
             5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
             3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
             15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
             30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 
             3.*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 
             10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 
             12.*Power(t0,3)*Power(t1,4)*x0 - 18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 
             36.*Power(t0,3)*Power(t1,3)*te*x0 + 54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 
             9.*Power(t0,4)*t1*Power(te,2)*x0 + 36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 
             36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 
             3.*Power(t0,5)*Power(t1,2)*xe + 15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 
             3.*Power(t1,7)*xe + 6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
             30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
             30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 
             3.*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 
             12.*Power(t0,2)*Power(t1,5)*xred + 3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 
             24.*Power(t0,3)*Power(t1,3)*te*xred + 6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 
             3.*Power(t0,5)*Power(te,2)*xred - 6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
             24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 
             12.*Power(t0,3)*t1*Power(te,3)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) + 
        (0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 10.*Power(t0,2)*Power(t1,6)*v0 + 
             5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 30.*Power(t0,3)*Power(t1,4)*te*v0 + 
             30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
             15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
             15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
             5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
             5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
             3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
             15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
             30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 
             3.*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 
             10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 
             12.*Power(t0,3)*Power(t1,4)*x0 - 18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 
             36.*Power(t0,3)*Power(t1,3)*te*x0 + 54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 
             9.*Power(t0,4)*t1*Power(te,2)*x0 + 36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 
             36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 
             3.*Power(t0,5)*Power(t1,2)*xe + 15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 
             3.*Power(t1,7)*xe + 6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
             30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
             30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 
             3.*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 
             12.*Power(t0,2)*Power(t1,5)*xred + 3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 
             24.*Power(t0,3)*Power(t1,3)*te*xred + 6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 
             3.*Power(t0,5)*Power(te,2)*xred - 6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
             24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 
             12.*Power(t0,3)*t1*Power(te,3)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
        (0.01388888888888889*(-3.9999999999999996*Power(t0,5)*Power(t1,3)*v0 + 20.*Power(t0,4)*Power(t1,4)*v0 - 40.*Power(t0,3)*Power(t1,5)*v0 + 
             40.*Power(t0,2)*Power(t1,6)*v0 - 20.*t0*Power(t1,7)*v0 + 3.999999999999999*Power(t1,8)*v0 + 12.*Power(t0,5)*Power(t1,2)*te*v0 - 
             60.*Power(t0,4)*Power(t1,3)*te*v0 + 119.99999999999999*Power(t0,3)*Power(t1,4)*te*v0 - 120.00000000000001*Power(t0,2)*Power(t1,5)*te*v0 + 
             60.00000000000001*t0*Power(t1,6)*te*v0 - 12.000000000000002*Power(t1,7)*te*v0 - 11.999999999999998*Power(t0,5)*t1*Power(te,2)*v0 + 
             60.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 - 120.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 + 120.00000000000003*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 - 
             60.000000000000014*t0*Power(t1,5)*Power(te,2)*v0 + 12.000000000000004*Power(t1,6)*Power(te,2)*v0 + 4.*Power(t0,5)*Power(te,3)*v0 - 
             20.*Power(t0,4)*t1*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 + 
             20.*t0*Power(t1,4)*Power(te,3)*v0 - 4.*Power(t1,5)*Power(te,3)*v0 - 3.*Power(t0,6)*Power(t1,2)*ve + 18.*Power(t0,5)*Power(t1,3)*ve - 
             44.99999999999999*Power(t0,4)*Power(t1,4)*ve + 59.99999999999999*Power(t0,3)*Power(t1,5)*ve - 44.99999999999999*Power(t0,2)*Power(t1,6)*ve + 
             18.*t0*Power(t1,7)*ve - 3.*Power(t1,8)*ve + 6.*Power(t0,6)*t1*te*ve - 36.*Power(t0,5)*Power(t1,2)*te*ve + 89.99999999999999*Power(t0,4)*Power(t1,3)*te*ve - 
             119.99999999999999*Power(t0,3)*Power(t1,4)*te*ve + 89.99999999999999*Power(t0,2)*Power(t1,5)*te*ve - 36.*t0*Power(t1,6)*te*ve + 6.*Power(t1,7)*te*ve - 
             3.*Power(t0,6)*Power(te,2)*ve + 18.*Power(t0,5)*t1*Power(te,2)*ve - 44.99999999999999*Power(t0,4)*Power(t1,2)*Power(te,2)*ve + 
             59.999999999999986*Power(t0,3)*Power(t1,3)*Power(te,2)*ve - 44.99999999999999*Power(t0,2)*Power(t1,4)*Power(te,2)*ve + 18.*t0*Power(t1,5)*Power(te,2)*ve - 
             3.*Power(t1,6)*Power(te,2)*ve + 12.*Power(t0,4)*Power(t1,3)*x0 - 48.*Power(t0,3)*Power(t1,4)*x0 + 72.*Power(t0,2)*Power(t1,5)*x0 - 48.*t0*Power(t1,6)*x0 + 
             12.*Power(t1,7)*x0 - 36.*Power(t0,4)*Power(t1,2)*te*x0 + 144.*Power(t0,3)*Power(t1,3)*te*x0 - 216.*Power(t0,2)*Power(t1,4)*te*x0 + 
             144.*t0*Power(t1,5)*te*x0 - 36.*Power(t1,6)*te*x0 + 36.*Power(t0,4)*t1*Power(te,2)*x0 - 144.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 + 
             216.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 - 144.00000000000003*t0*Power(t1,4)*Power(te,2)*x0 + 36.00000000000001*Power(t1,5)*Power(te,2)*x0 - 
             12.*Power(t0,4)*Power(te,3)*x0 + 48.*Power(t0,3)*t1*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 + 48.*t0*Power(t1,3)*Power(te,3)*x0 - 
             12.*Power(t1,4)*Power(te,3)*x0 - 6.*Power(t0,6)*t1*xe + 36.*Power(t0,5)*Power(t1,2)*xe - 90.*Power(t0,4)*Power(t1,3)*xe + 120.*Power(t0,3)*Power(t1,4)*xe - 
             90.*Power(t0,2)*Power(t1,5)*xe + 36.*t0*Power(t1,6)*xe - 6.*Power(t1,7)*xe + 6.*Power(t0,6)*te*xe - 36.*Power(t0,5)*t1*te*xe + 
             90.*Power(t0,4)*Power(t1,2)*te*xe - 120.*Power(t0,3)*Power(t1,3)*te*xe + 90.*Power(t0,2)*Power(t1,4)*te*xe - 36.*t0*Power(t1,5)*te*xe + 
             6.*Power(t1,6)*te*xe + 6.*Power(t0,6)*t1*xred - 36.*Power(t0,5)*Power(t1,2)*xred + 78.*Power(t0,4)*Power(t1,3)*xred - 72.*Power(t0,3)*Power(t1,4)*xred + 
             18.*Power(t0,2)*Power(t1,5)*xred + 11.999999999999998*t0*Power(t1,6)*xred - 5.999999999999999*Power(t1,7)*xred - 6.*Power(t0,6)*te*xred + 
             36.*Power(t0,5)*t1*te*xred - 53.99999999999999*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 
             126.*Power(t0,2)*Power(t1,4)*te*xred - 108.00000000000001*t0*Power(t1,5)*te*xred + 30.*Power(t1,6)*te*xred - 36.*Power(t0,4)*t1*Power(te,2)*xred + 
             144.*Power(t0,3)*Power(t1,2)*Power(te,2)*xred - 216.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred + 144.00000000000003*t0*Power(t1,4)*Power(te,2)*xred - 
             36.00000000000001*Power(t1,5)*Power(te,2)*xred + 12.*Power(t0,4)*Power(te,3)*xred - 48.*Power(t0,3)*t1*Power(te,3)*xred + 
             72.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred - 48.*t0*Power(t1,3)*Power(te,3)*xred + 12.*Power(t1,4)*Power(te,3)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
        (0.041666666666666664*Power(-t1 + te,2)*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 
             0.16666666666666666*Power(t0,2)*Power(t1,4) - 0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 
             0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 
             0.08333333333333334*Power(t0,4)*Power(te,2) - 0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 
             0.25*Power(t1,4)*Power(te,2) - 0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 
             0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 0.11111111111111113*Power(t1,3)*Power(te,3))*
           (-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 
             1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
             40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
             30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
             30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
             40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
             4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
             10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
             4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
             1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
             5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
             5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
             15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
             42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
             1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
             65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
             2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 
             20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 
             12.*Power(t0,3)*Power(t1,5)*x0 + 18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 
             48.*Power(t0,3)*Power(t1,4)*te*x0 - 72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 
             18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 
             72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 
             20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 
             30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 
             10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
             100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 
             11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 
             40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 
             8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 
             18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 
             1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 
             28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 
             7.*Power(t1,6)*Power(te,2)*xred + 4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 
             8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 
             28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 3.*Power(t0,4)*Power(te,4)*xred + 
             12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
           Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
             0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
             0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
             0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
             0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
             0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
             0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) - 
        (0.08333333333333333*(-t1 + te)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
             10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
             40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
             60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
             20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
             20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
             5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
             4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
             1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
             5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
             5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
             15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
             42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
             1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
             65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
             2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 
             20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 
             12.*Power(t0,3)*Power(t1,5)*x0 + 18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 
             48.*Power(t0,3)*Power(t1,4)*te*x0 - 72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 
             18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 
             72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 
             20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 
             30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 
             10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
             100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 
             11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 
             40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 
             8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 
             18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 
             1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 
             28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 
             7.*Power(t1,6)*Power(te,2)*xred + 4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 
             8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 
             28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 3.*Power(t0,4)*Power(te,4)*xred + 
             12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
           (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
             0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
             0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
             0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
             0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
             0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
             0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
        (0.08333333333333333*Power(-t1 + te,2)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
             10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
             40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
             60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
             20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
             20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
             5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
             4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
             1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
             5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
             5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
             15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
             42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
             1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
             65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
             2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 
             20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 
             12.*Power(t0,3)*Power(t1,5)*x0 + 18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 
             48.*Power(t0,3)*Power(t1,4)*te*x0 - 72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 
             18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 
             72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 
             20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 
             30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 
             10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
             100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 
             11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 
             40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 
             8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 
             18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 
             1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 
             28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 
             7.*Power(t1,6)*Power(te,2)*xred + 4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 
             8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 
             28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 3.*Power(t0,4)*Power(te,4)*xred + 
             12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,3)*
           (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
             0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
             0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
             0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
             0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
             0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
             0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
        (0.01388888888888889*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 0.16666666666666666*Power(t0,2)*Power(t1,4) - 
             0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 
             0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 0.08333333333333334*Power(t0,4)*Power(te,2) - 
             0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 0.25*Power(t1,4)*Power(te,2) - 
             0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 
             0.11111111111111113*Power(t1,3)*Power(te,3))*(1.*Power(t0,5)*Power(t1,4)*v0 - 5.*Power(t0,4)*Power(t1,5)*v0 + 10.*Power(t0,3)*Power(t1,6)*v0 - 
             10.*Power(t0,2)*Power(t1,7)*v0 + 5.*t0*Power(t1,8)*v0 - 1.*Power(t1,9)*v0 - 3.9999999999999996*Power(t0,5)*Power(t1,3)*te*v0 + 
             20.*Power(t0,4)*Power(t1,4)*te*v0 - 40.*Power(t0,3)*Power(t1,5)*te*v0 + 40.*Power(t0,2)*Power(t1,6)*te*v0 - 20.*t0*Power(t1,7)*te*v0 + 
             3.999999999999999*Power(t1,8)*te*v0 + 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 - 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 + 
             59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 - 60.00000000000001*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 + 
             30.000000000000004*t0*Power(t1,6)*Power(te,2)*v0 - 6.000000000000001*Power(t1,7)*Power(te,2)*v0 - 3.9999999999999996*Power(t0,5)*t1*Power(te,3)*v0 + 
             20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 - 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 + 40.00000000000001*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 - 
             20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 + 4.000000000000001*Power(t1,6)*Power(te,3)*v0 + 1.*Power(t0,5)*Power(te,4)*v0 - 
             5.*Power(t0,4)*t1*Power(te,4)*v0 + 10.*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 - 10.*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 + 
             5.*t0*Power(t1,4)*Power(te,4)*v0 - 1.*Power(t1,5)*Power(te,4)*v0 + 1.*Power(t0,6)*Power(t1,3)*ve - 6.*Power(t0,5)*Power(t1,4)*ve + 
             14.999999999999998*Power(t0,4)*Power(t1,5)*ve - 19.999999999999996*Power(t0,3)*Power(t1,6)*ve + 14.999999999999998*Power(t0,2)*Power(t1,7)*ve - 
             6.*t0*Power(t1,8)*ve + 1.*Power(t1,9)*ve - 3.*Power(t0,6)*Power(t1,2)*te*ve + 18.*Power(t0,5)*Power(t1,3)*te*ve - 
             44.99999999999999*Power(t0,4)*Power(t1,4)*te*ve + 59.99999999999999*Power(t0,3)*Power(t1,5)*te*ve - 44.99999999999999*Power(t0,2)*Power(t1,6)*te*ve + 
             18.*t0*Power(t1,7)*te*ve - 3.*Power(t1,8)*te*ve + 3.*Power(t0,6)*t1*Power(te,2)*ve - 18.*Power(t0,5)*Power(t1,2)*Power(te,2)*ve + 
             44.99999999999999*Power(t0,4)*Power(t1,3)*Power(te,2)*ve - 59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*ve + 
             44.99999999999999*Power(t0,2)*Power(t1,5)*Power(te,2)*ve - 18.*t0*Power(t1,6)*Power(te,2)*ve + 3.*Power(t1,7)*Power(te,2)*ve - 
             1.*Power(t0,6)*Power(te,3)*ve + 6.*Power(t0,5)*t1*Power(te,3)*ve - 14.999999999999998*Power(t0,4)*Power(t1,2)*Power(te,3)*ve + 
             19.999999999999996*Power(t0,3)*Power(t1,3)*Power(te,3)*ve - 14.999999999999998*Power(t0,2)*Power(t1,4)*Power(te,3)*ve + 6.*t0*Power(t1,5)*Power(te,3)*ve - 
             1.*Power(t1,6)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,4)*x0 + 12.*Power(t0,3)*Power(t1,5)*x0 - 18.*Power(t0,2)*Power(t1,6)*x0 + 12.*t0*Power(t1,7)*x0 - 
             3.*Power(t1,8)*x0 + 12.*Power(t0,4)*Power(t1,3)*te*x0 - 48.*Power(t0,3)*Power(t1,4)*te*x0 + 72.*Power(t0,2)*Power(t1,5)*te*x0 - 48.*t0*Power(t1,6)*te*x0 + 
             12.*Power(t1,7)*te*x0 - 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 + 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 - 
             108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 + 72.*t0*Power(t1,5)*Power(te,2)*x0 - 18.*Power(t1,6)*Power(te,2)*x0 + 12.*Power(t0,4)*t1*Power(te,3)*x0 - 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 + 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,4)*Power(te,3)*x0 + 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 - 3.*Power(t0,4)*Power(te,4)*x0 + 12.*Power(t0,3)*t1*Power(te,4)*x0 - 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 + 12.*t0*Power(t1,3)*Power(te,4)*x0 - 3.*Power(t1,4)*Power(te,4)*x0 + 3.*Power(t0,6)*Power(t1,2)*xe - 
             18.*Power(t0,5)*Power(t1,3)*xe + 45.*Power(t0,4)*Power(t1,4)*xe - 60.*Power(t0,3)*Power(t1,5)*xe + 45.*Power(t0,2)*Power(t1,6)*xe - 18.*t0*Power(t1,7)*xe + 
             3.*Power(t1,8)*xe - 6.*Power(t0,6)*t1*te*xe + 36.*Power(t0,5)*Power(t1,2)*te*xe - 90.*Power(t0,4)*Power(t1,3)*te*xe + 120.*Power(t0,3)*Power(t1,4)*te*xe - 
             90.*Power(t0,2)*Power(t1,5)*te*xe + 36.*t0*Power(t1,6)*te*xe - 6.*Power(t1,7)*te*xe + 3.*Power(t0,6)*Power(te,2)*xe - 18.*Power(t0,5)*t1*Power(te,2)*xe + 
             45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe - 60.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe + 45.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe - 
             18.*t0*Power(t1,5)*Power(te,2)*xe + 3.*Power(t1,6)*Power(te,2)*xe - 3.*Power(t0,6)*Power(t1,2)*xred + 18.*Power(t0,5)*Power(t1,3)*xred - 
             42.*Power(t0,4)*Power(t1,4)*xred + 48.*Power(t0,3)*Power(t1,5)*xred - 26.999999999999996*Power(t0,2)*Power(t1,6)*xred + 
             5.999999999999999*t0*Power(t1,7)*xred + 6.*Power(t0,6)*t1*te*xred - 36.*Power(t0,5)*Power(t1,2)*te*xred + 78.*Power(t0,4)*Power(t1,3)*te*xred - 
             72.*Power(t0,3)*Power(t1,4)*te*xred + 18.*Power(t0,2)*Power(t1,5)*te*xred + 11.999999999999998*t0*Power(t1,6)*te*xred - 
             5.999999999999999*Power(t1,7)*te*xred - 3.*Power(t0,6)*Power(te,2)*xred + 18.*Power(t0,5)*t1*Power(te,2)*xred - 
             26.999999999999996*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 12.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred + 63.*Power(t0,2)*Power(t1,4)*Power(te,2)*xred - 
             54.00000000000001*t0*Power(t1,5)*Power(te,2)*xred + 15.*Power(t1,6)*Power(te,2)*xred - 12.*Power(t0,4)*t1*Power(te,3)*xred + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred + 48.00000000000001*t0*Power(t1,4)*Power(te,3)*xred - 
             12.000000000000002*Power(t1,5)*Power(te,3)*xred + 3.*Power(t0,4)*Power(te,4)*xred - 12.*Power(t0,3)*t1*Power(te,4)*xred + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred - 12.*t0*Power(t1,3)*Power(te,4)*xred + 3.*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2))))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
      (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
        0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
        0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
        0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
        0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
        0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
        0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
   (0.08333333333333333*(4.*Power(t0,5)*Power(t1,3)*v0 - 20.*Power(t0,4)*Power(t1,4)*v0 + 40.*Power(t0,3)*Power(t1,5)*v0 - 40.*Power(t0,2)*Power(t1,6)*v0 + 
        20.*t0*Power(t1,7)*v0 - 4.*Power(t1,8)*v0 - 12.*Power(t0,5)*Power(t1,2)*te*v0 + 60.*Power(t0,4)*Power(t1,3)*te*v0 - 120.*Power(t0,3)*Power(t1,4)*te*v0 + 
        120.*Power(t0,2)*Power(t1,5)*te*v0 - 60.*t0*Power(t1,6)*te*v0 + 12.*Power(t1,7)*te*v0 + 12.*Power(t0,5)*t1*Power(te,2)*v0 - 
        60.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 120.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 120.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
        60.000000000000014*t0*Power(t1,5)*Power(te,2)*v0 - 12.000000000000004*Power(t1,6)*Power(te,2)*v0 - 4.*Power(t0,5)*Power(te,3)*v0 + 
        20.*Power(t0,4)*t1*Power(te,3)*v0 - 40.00000000000001*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 40.00000000000001*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
        20.000000000000007*t0*Power(t1,4)*Power(te,3)*v0 + 4.000000000000003*Power(t1,5)*Power(te,3)*v0 - 3.0000000000000004*Power(t0,6)*Power(t1,2)*ve + 
        10.000000000000002*Power(t0,5)*Power(t1,3)*ve - 5.000000000000001*Power(t0,4)*Power(t1,4)*ve - 20.*Power(t0,3)*Power(t1,5)*ve + 35.*Power(t0,2)*Power(t1,6)*ve - 
        22.*t0*Power(t1,7)*ve + 5.*Power(t1,8)*ve + 6.000000000000001*Power(t0,6)*t1*te*ve - 12.000000000000002*Power(t0,5)*Power(t1,2)*te*ve - 
        30.*Power(t0,4)*Power(t1,3)*te*ve + 120.*Power(t0,3)*Power(t1,4)*te*ve - 150.*Power(t0,2)*Power(t1,5)*te*ve + 84.*t0*Power(t1,6)*te*ve - 18.*Power(t1,7)*te*ve - 
        3.000000000000001*Power(t0,6)*Power(te,2)*ve - 5.999999999999999*Power(t0,5)*t1*Power(te,2)*ve + 75.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
        180.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 195.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 102.*t0*Power(t1,5)*Power(te,2)*ve + 
        21.000000000000004*Power(t1,6)*Power(te,2)*ve + 8.*Power(t0,5)*Power(te,3)*ve - 40.*Power(t0,4)*t1*Power(te,3)*ve + 80.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 
        80.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 40.*t0*Power(t1,4)*Power(te,3)*ve - 8.*Power(t1,5)*Power(te,3)*ve - 12.*Power(t0,4)*Power(t1,3)*x0 + 
        48.*Power(t0,3)*Power(t1,4)*x0 - 72.*Power(t0,2)*Power(t1,5)*x0 + 48.*t0*Power(t1,6)*x0 - 12.*Power(t1,7)*x0 + 36.*Power(t0,4)*Power(t1,2)*te*x0 - 
        144.*Power(t0,3)*Power(t1,3)*te*x0 + 216.*Power(t0,2)*Power(t1,4)*te*x0 - 144.*t0*Power(t1,5)*te*x0 + 36.*Power(t1,6)*te*x0 - 
        36.00000000000001*Power(t0,4)*t1*Power(te,2)*x0 + 144.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 216.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 
        144.*t0*Power(t1,4)*Power(te,2)*x0 - 36.00000000000001*Power(t1,5)*Power(te,2)*x0 + 12.*Power(t0,4)*Power(te,3)*x0 - 48.*Power(t0,3)*t1*Power(te,3)*x0 + 
        72.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,3)*Power(te,3)*x0 + 12.000000000000004*Power(t1,4)*Power(te,3)*x0 - 
        2.0000000000000004*Power(t0,6)*t1*xe + 30.*Power(t0,4)*Power(t1,3)*xe - 80.*Power(t0,3)*Power(t1,4)*xe + 90.*Power(t0,2)*Power(t1,5)*xe - 
        48.*t0*Power(t1,6)*xe + 10.*Power(t1,7)*xe + 2.0000000000000004*Power(t0,6)*te*xe + 12.*Power(t0,5)*t1*te*xe - 90.*Power(t0,4)*Power(t1,2)*te*xe + 
        200.*Power(t0,3)*Power(t1,3)*te*xe - 210.*Power(t0,2)*Power(t1,4)*te*xe + 108.*t0*Power(t1,5)*te*xe - 22.*Power(t1,6)*te*xe - 12.*Power(t0,5)*Power(te,2)*xe + 
        60.*Power(t0,4)*t1*Power(te,2)*xe - 120.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 120.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 
        60.*t0*Power(t1,4)*Power(te,2)*xe + 12.*Power(t1,5)*Power(te,2)*xe + 2.0000000000000004*Power(t0,6)*t1*xred - 18.*Power(t0,4)*Power(t1,3)*xred + 
        32.*Power(t0,3)*Power(t1,4)*xred - 18.*Power(t0,2)*Power(t1,5)*xred + 2.*Power(t1,7)*xred - 2.0000000000000004*Power(t0,6)*te*xred - 
        12.*Power(t0,5)*t1*te*xred + 54.00000000000001*Power(t0,4)*Power(t1,2)*te*xred - 56.*Power(t0,3)*Power(t1,3)*te*xred - 
        5.9999999999999964*Power(t0,2)*Power(t1,4)*te*xred + 36.*t0*Power(t1,5)*te*xred - 14.*Power(t1,6)*te*xred + 12.*Power(t0,5)*Power(te,2)*xred - 
        24.*Power(t0,4)*t1*Power(te,2)*xred - 24.000000000000007*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 96.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 
        84.00000000000001*t0*Power(t1,4)*Power(te,2)*xred + 24.*Power(t1,5)*Power(te,2)*xred - 12.*Power(t0,4)*Power(te,3)*xred + 48.*Power(t0,3)*t1*Power(te,3)*xred - 
        72.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 48.00000000000001*t0*Power(t1,3)*Power(te,3)*xred - 12.000000000000004*Power(t1,4)*Power(te,3)*xred)*
      ((0.05555555555555555*(-t1 + te)*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 
             10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 
             30.*Power(t0,3)*Power(t1,4)*te*v0 + 30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 
             3.*Power(t0,5)*t1*Power(te,2)*v0 - 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
             30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 
             1.*Power(t0,5)*Power(te,3)*v0 + 5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
             10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 
             1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 
             10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 
             30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 
             3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 
             30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 
             5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
             5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
             18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
             54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
             36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 
             9.*Power(t1,5)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
             12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
             15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
             6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
             30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
             30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 
             3.*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 
             12.*Power(t0,2)*Power(t1,5)*xred + 3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 
             24.*Power(t0,3)*Power(t1,3)*te*xred + 6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 
             3.*Power(t0,5)*Power(te,2)*xred - 6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
             24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 
             12.*Power(t0,3)*t1*Power(te,3)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
        (0.041666666666666664*Power(-t1 + te,2)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
             10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
             40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
             60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
             20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
             20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
             5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
             4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
             1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
             5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
             5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
             15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
             42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
             1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
             65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
             2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 
             20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 
             12.*Power(t0,3)*Power(t1,5)*x0 + 18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 
             48.*Power(t0,3)*Power(t1,4)*te*x0 - 72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 
             18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 
             72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 
             20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 
             30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 
             10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
             100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 
             11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 
             40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 
             8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 
             18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 
             1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 
             28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 
             7.*Power(t1,6)*Power(te,2)*xred + 4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 
             8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 
             28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 3.*Power(t0,4)*Power(te,4)*xred + 
             12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
           (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
             0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
             0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
             0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
             0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
             0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
             0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
        (0.01388888888888889*(1.*Power(t0,5)*Power(t1,4)*v0 - 5.*Power(t0,4)*Power(t1,5)*v0 + 10.*Power(t0,3)*Power(t1,6)*v0 - 10.*Power(t0,2)*Power(t1,7)*v0 + 
             5.*t0*Power(t1,8)*v0 - 1.*Power(t1,9)*v0 - 3.9999999999999996*Power(t0,5)*Power(t1,3)*te*v0 + 20.*Power(t0,4)*Power(t1,4)*te*v0 - 
             40.*Power(t0,3)*Power(t1,5)*te*v0 + 40.*Power(t0,2)*Power(t1,6)*te*v0 - 20.*t0*Power(t1,7)*te*v0 + 3.999999999999999*Power(t1,8)*te*v0 + 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 - 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 + 59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 - 
             60.00000000000001*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 + 30.000000000000004*t0*Power(t1,6)*Power(te,2)*v0 - 
             6.000000000000001*Power(t1,7)*Power(te,2)*v0 - 3.9999999999999996*Power(t0,5)*t1*Power(te,3)*v0 + 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 - 
             40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 + 40.00000000000001*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 - 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 + 
             4.000000000000001*Power(t1,6)*Power(te,3)*v0 + 1.*Power(t0,5)*Power(te,4)*v0 - 5.*Power(t0,4)*t1*Power(te,4)*v0 + 
             10.*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 - 10.*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 + 5.*t0*Power(t1,4)*Power(te,4)*v0 - 
             1.*Power(t1,5)*Power(te,4)*v0 + 1.*Power(t0,6)*Power(t1,3)*ve - 6.*Power(t0,5)*Power(t1,4)*ve + 14.999999999999998*Power(t0,4)*Power(t1,5)*ve - 
             19.999999999999996*Power(t0,3)*Power(t1,6)*ve + 14.999999999999998*Power(t0,2)*Power(t1,7)*ve - 6.*t0*Power(t1,8)*ve + 1.*Power(t1,9)*ve - 
             3.*Power(t0,6)*Power(t1,2)*te*ve + 18.*Power(t0,5)*Power(t1,3)*te*ve - 44.99999999999999*Power(t0,4)*Power(t1,4)*te*ve + 
             59.99999999999999*Power(t0,3)*Power(t1,5)*te*ve - 44.99999999999999*Power(t0,2)*Power(t1,6)*te*ve + 18.*t0*Power(t1,7)*te*ve - 3.*Power(t1,8)*te*ve + 
             3.*Power(t0,6)*t1*Power(te,2)*ve - 18.*Power(t0,5)*Power(t1,2)*Power(te,2)*ve + 44.99999999999999*Power(t0,4)*Power(t1,3)*Power(te,2)*ve - 
             59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*ve + 44.99999999999999*Power(t0,2)*Power(t1,5)*Power(te,2)*ve - 18.*t0*Power(t1,6)*Power(te,2)*ve + 
             3.*Power(t1,7)*Power(te,2)*ve - 1.*Power(t0,6)*Power(te,3)*ve + 6.*Power(t0,5)*t1*Power(te,3)*ve - 
             14.999999999999998*Power(t0,4)*Power(t1,2)*Power(te,3)*ve + 19.999999999999996*Power(t0,3)*Power(t1,3)*Power(te,3)*ve - 
             14.999999999999998*Power(t0,2)*Power(t1,4)*Power(te,3)*ve + 6.*t0*Power(t1,5)*Power(te,3)*ve - 1.*Power(t1,6)*Power(te,3)*ve - 
             3.*Power(t0,4)*Power(t1,4)*x0 + 12.*Power(t0,3)*Power(t1,5)*x0 - 18.*Power(t0,2)*Power(t1,6)*x0 + 12.*t0*Power(t1,7)*x0 - 3.*Power(t1,8)*x0 + 
             12.*Power(t0,4)*Power(t1,3)*te*x0 - 48.*Power(t0,3)*Power(t1,4)*te*x0 + 72.*Power(t0,2)*Power(t1,5)*te*x0 - 48.*t0*Power(t1,6)*te*x0 + 
             12.*Power(t1,7)*te*x0 - 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 + 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 - 
             108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 + 72.*t0*Power(t1,5)*Power(te,2)*x0 - 18.*Power(t1,6)*Power(te,2)*x0 + 12.*Power(t0,4)*t1*Power(te,3)*x0 - 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 + 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,4)*Power(te,3)*x0 + 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 - 3.*Power(t0,4)*Power(te,4)*x0 + 12.*Power(t0,3)*t1*Power(te,4)*x0 - 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 + 12.*t0*Power(t1,3)*Power(te,4)*x0 - 3.*Power(t1,4)*Power(te,4)*x0 + 3.*Power(t0,6)*Power(t1,2)*xe - 
             18.*Power(t0,5)*Power(t1,3)*xe + 45.*Power(t0,4)*Power(t1,4)*xe - 60.*Power(t0,3)*Power(t1,5)*xe + 45.*Power(t0,2)*Power(t1,6)*xe - 18.*t0*Power(t1,7)*xe + 
             3.*Power(t1,8)*xe - 6.*Power(t0,6)*t1*te*xe + 36.*Power(t0,5)*Power(t1,2)*te*xe - 90.*Power(t0,4)*Power(t1,3)*te*xe + 120.*Power(t0,3)*Power(t1,4)*te*xe - 
             90.*Power(t0,2)*Power(t1,5)*te*xe + 36.*t0*Power(t1,6)*te*xe - 6.*Power(t1,7)*te*xe + 3.*Power(t0,6)*Power(te,2)*xe - 18.*Power(t0,5)*t1*Power(te,2)*xe + 
             45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe - 60.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe + 45.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe - 
             18.*t0*Power(t1,5)*Power(te,2)*xe + 3.*Power(t1,6)*Power(te,2)*xe - 3.*Power(t0,6)*Power(t1,2)*xred + 18.*Power(t0,5)*Power(t1,3)*xred - 
             42.*Power(t0,4)*Power(t1,4)*xred + 48.*Power(t0,3)*Power(t1,5)*xred - 26.999999999999996*Power(t0,2)*Power(t1,6)*xred + 
             5.999999999999999*t0*Power(t1,7)*xred + 6.*Power(t0,6)*t1*te*xred - 36.*Power(t0,5)*Power(t1,2)*te*xred + 78.*Power(t0,4)*Power(t1,3)*te*xred - 
             72.*Power(t0,3)*Power(t1,4)*te*xred + 18.*Power(t0,2)*Power(t1,5)*te*xred + 11.999999999999998*t0*Power(t1,6)*te*xred - 
             5.999999999999999*Power(t1,7)*te*xred - 3.*Power(t0,6)*Power(te,2)*xred + 18.*Power(t0,5)*t1*Power(te,2)*xred - 
             26.999999999999996*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 12.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred + 63.*Power(t0,2)*Power(t1,4)*Power(te,2)*xred - 
             54.00000000000001*t0*Power(t1,5)*Power(te,2)*xred + 15.*Power(t1,6)*Power(te,2)*xred - 12.*Power(t0,4)*t1*Power(te,3)*xred + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred + 48.00000000000001*t0*Power(t1,4)*Power(te,3)*xred - 
             12.000000000000002*Power(t1,5)*Power(te,3)*xred + 3.*Power(t0,4)*Power(te,4)*xred - 12.*Power(t0,3)*t1*Power(te,4)*xred + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred - 12.*t0*Power(t1,3)*Power(te,4)*xred + 3.*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)))))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
      (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
        0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
        0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
        0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
        0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
        0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
        0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
   (0.08333333333333333*(0.08333333333333334*Power(t0,4)*Power(t1,2) - 0.22222222222222224*Power(t0,3)*Power(t1,3) + 0.16666666666666666*Power(t0,2)*Power(t1,4) - 
        0.027777777777777776*Power(t1,6) - 0.16666666666666669*Power(t0,4)*t1*te + 0.3333333333333333*Power(t0,3)*Power(t1,2)*te - 
        0.3333333333333334*t0*Power(t1,4)*te + 0.16666666666666669*Power(t1,5)*te + 0.08333333333333334*Power(t0,4)*Power(te,2) - 
        0.5*Power(t0,2)*Power(t1,2)*Power(te,2) + 0.6666666666666667*t0*Power(t1,3)*Power(te,2) - 0.25*Power(t1,4)*Power(te,2) - 
        0.11111111111111113*Power(t0,3)*Power(te,3) + 0.33333333333333337*Power(t0,2)*t1*Power(te,3) - 0.33333333333333337*t0*Power(t1,2)*Power(te,3) + 
        0.11111111111111113*Power(t1,3)*Power(te,3))*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
        10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
        40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
        6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
        60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
        20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
        20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
        5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
        5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
        4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
        1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
        5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
        5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
        15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
        42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
        1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
        65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 2.*Power(t0,5)*Power(te,4)*ve - 
        10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
        10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
        18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
        72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
        72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 
        12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 
        48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
        18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
        1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 
        25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 
        80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 
        1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
        100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 
        4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
        20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
        2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
        2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
        18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
        27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
        2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
        4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
        32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
        3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
        12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred)*
      ((0.05555555555555555*(-t1 + te)*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 
             10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 
             30.*Power(t0,3)*Power(t1,4)*te*v0 + 30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 
             3.*Power(t0,5)*t1*Power(te,2)*v0 - 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
             30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 
             1.*Power(t0,5)*Power(te,3)*v0 + 5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
             10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 
             1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 
             10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 
             30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 
             3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 
             30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 
             5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
             5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
             18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
             54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
             36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 
             9.*Power(t1,5)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
             12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
             15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
             6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
             30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
             30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 
             3.*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 
             12.*Power(t0,2)*Power(t1,5)*xred + 3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 
             24.*Power(t0,3)*Power(t1,3)*te*xred + 6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 
             3.*Power(t0,5)*Power(te,2)*xred - 6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
             24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 
             12.*Power(t0,3)*t1*Power(te,3)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
        (0.041666666666666664*Power(-t1 + te,2)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
             10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
             40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
             60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
             20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
             20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
             5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
             4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
             1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
             5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
             5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
             15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
             42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
             1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
             65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
             2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 
             20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 
             12.*Power(t0,3)*Power(t1,5)*x0 + 18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 
             48.*Power(t0,3)*Power(t1,4)*te*x0 - 72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 
             18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 
             72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 
             20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 
             30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 
             10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
             100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 
             11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 
             40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 
             8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 
             18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 
             1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 
             28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 
             7.*Power(t1,6)*Power(te,2)*xred + 4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 
             8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 
             28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 3.*Power(t0,4)*Power(te,4)*xred + 
             12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
           (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
             0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
             0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
             0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
             0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
             0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
             0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
        (0.01388888888888889*(1.*Power(t0,5)*Power(t1,4)*v0 - 5.*Power(t0,4)*Power(t1,5)*v0 + 10.*Power(t0,3)*Power(t1,6)*v0 - 10.*Power(t0,2)*Power(t1,7)*v0 + 
             5.*t0*Power(t1,8)*v0 - 1.*Power(t1,9)*v0 - 3.9999999999999996*Power(t0,5)*Power(t1,3)*te*v0 + 20.*Power(t0,4)*Power(t1,4)*te*v0 - 
             40.*Power(t0,3)*Power(t1,5)*te*v0 + 40.*Power(t0,2)*Power(t1,6)*te*v0 - 20.*t0*Power(t1,7)*te*v0 + 3.999999999999999*Power(t1,8)*te*v0 + 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 - 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 + 59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 - 
             60.00000000000001*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 + 30.000000000000004*t0*Power(t1,6)*Power(te,2)*v0 - 
             6.000000000000001*Power(t1,7)*Power(te,2)*v0 - 3.9999999999999996*Power(t0,5)*t1*Power(te,3)*v0 + 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 - 
             40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 + 40.00000000000001*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 - 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 + 
             4.000000000000001*Power(t1,6)*Power(te,3)*v0 + 1.*Power(t0,5)*Power(te,4)*v0 - 5.*Power(t0,4)*t1*Power(te,4)*v0 + 
             10.*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 - 10.*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 + 5.*t0*Power(t1,4)*Power(te,4)*v0 - 
             1.*Power(t1,5)*Power(te,4)*v0 + 1.*Power(t0,6)*Power(t1,3)*ve - 6.*Power(t0,5)*Power(t1,4)*ve + 14.999999999999998*Power(t0,4)*Power(t1,5)*ve - 
             19.999999999999996*Power(t0,3)*Power(t1,6)*ve + 14.999999999999998*Power(t0,2)*Power(t1,7)*ve - 6.*t0*Power(t1,8)*ve + 1.*Power(t1,9)*ve - 
             3.*Power(t0,6)*Power(t1,2)*te*ve + 18.*Power(t0,5)*Power(t1,3)*te*ve - 44.99999999999999*Power(t0,4)*Power(t1,4)*te*ve + 
             59.99999999999999*Power(t0,3)*Power(t1,5)*te*ve - 44.99999999999999*Power(t0,2)*Power(t1,6)*te*ve + 18.*t0*Power(t1,7)*te*ve - 3.*Power(t1,8)*te*ve + 
             3.*Power(t0,6)*t1*Power(te,2)*ve - 18.*Power(t0,5)*Power(t1,2)*Power(te,2)*ve + 44.99999999999999*Power(t0,4)*Power(t1,3)*Power(te,2)*ve - 
             59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*ve + 44.99999999999999*Power(t0,2)*Power(t1,5)*Power(te,2)*ve - 18.*t0*Power(t1,6)*Power(te,2)*ve + 
             3.*Power(t1,7)*Power(te,2)*ve - 1.*Power(t0,6)*Power(te,3)*ve + 6.*Power(t0,5)*t1*Power(te,3)*ve - 
             14.999999999999998*Power(t0,4)*Power(t1,2)*Power(te,3)*ve + 19.999999999999996*Power(t0,3)*Power(t1,3)*Power(te,3)*ve - 
             14.999999999999998*Power(t0,2)*Power(t1,4)*Power(te,3)*ve + 6.*t0*Power(t1,5)*Power(te,3)*ve - 1.*Power(t1,6)*Power(te,3)*ve - 
             3.*Power(t0,4)*Power(t1,4)*x0 + 12.*Power(t0,3)*Power(t1,5)*x0 - 18.*Power(t0,2)*Power(t1,6)*x0 + 12.*t0*Power(t1,7)*x0 - 3.*Power(t1,8)*x0 + 
             12.*Power(t0,4)*Power(t1,3)*te*x0 - 48.*Power(t0,3)*Power(t1,4)*te*x0 + 72.*Power(t0,2)*Power(t1,5)*te*x0 - 48.*t0*Power(t1,6)*te*x0 + 
             12.*Power(t1,7)*te*x0 - 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 + 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 - 
             108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 + 72.*t0*Power(t1,5)*Power(te,2)*x0 - 18.*Power(t1,6)*Power(te,2)*x0 + 12.*Power(t0,4)*t1*Power(te,3)*x0 - 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 + 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,4)*Power(te,3)*x0 + 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 - 3.*Power(t0,4)*Power(te,4)*x0 + 12.*Power(t0,3)*t1*Power(te,4)*x0 - 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 + 12.*t0*Power(t1,3)*Power(te,4)*x0 - 3.*Power(t1,4)*Power(te,4)*x0 + 3.*Power(t0,6)*Power(t1,2)*xe - 
             18.*Power(t0,5)*Power(t1,3)*xe + 45.*Power(t0,4)*Power(t1,4)*xe - 60.*Power(t0,3)*Power(t1,5)*xe + 45.*Power(t0,2)*Power(t1,6)*xe - 18.*t0*Power(t1,7)*xe + 
             3.*Power(t1,8)*xe - 6.*Power(t0,6)*t1*te*xe + 36.*Power(t0,5)*Power(t1,2)*te*xe - 90.*Power(t0,4)*Power(t1,3)*te*xe + 120.*Power(t0,3)*Power(t1,4)*te*xe - 
             90.*Power(t0,2)*Power(t1,5)*te*xe + 36.*t0*Power(t1,6)*te*xe - 6.*Power(t1,7)*te*xe + 3.*Power(t0,6)*Power(te,2)*xe - 18.*Power(t0,5)*t1*Power(te,2)*xe + 
             45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe - 60.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe + 45.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe - 
             18.*t0*Power(t1,5)*Power(te,2)*xe + 3.*Power(t1,6)*Power(te,2)*xe - 3.*Power(t0,6)*Power(t1,2)*xred + 18.*Power(t0,5)*Power(t1,3)*xred - 
             42.*Power(t0,4)*Power(t1,4)*xred + 48.*Power(t0,3)*Power(t1,5)*xred - 26.999999999999996*Power(t0,2)*Power(t1,6)*xred + 
             5.999999999999999*t0*Power(t1,7)*xred + 6.*Power(t0,6)*t1*te*xred - 36.*Power(t0,5)*Power(t1,2)*te*xred + 78.*Power(t0,4)*Power(t1,3)*te*xred - 
             72.*Power(t0,3)*Power(t1,4)*te*xred + 18.*Power(t0,2)*Power(t1,5)*te*xred + 11.999999999999998*t0*Power(t1,6)*te*xred - 
             5.999999999999999*Power(t1,7)*te*xred - 3.*Power(t0,6)*Power(te,2)*xred + 18.*Power(t0,5)*t1*Power(te,2)*xred - 
             26.999999999999996*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 12.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred + 63.*Power(t0,2)*Power(t1,4)*Power(te,2)*xred - 
             54.00000000000001*t0*Power(t1,5)*Power(te,2)*xred + 15.*Power(t1,6)*Power(te,2)*xred - 12.*Power(t0,4)*t1*Power(te,3)*xred + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred + 48.00000000000001*t0*Power(t1,4)*Power(te,3)*xred - 
             12.000000000000002*Power(t1,5)*Power(te,3)*xred + 3.*Power(t0,4)*Power(te,4)*xred - 12.*Power(t0,3)*t1*Power(te,4)*xred + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred - 12.*t0*Power(t1,3)*Power(te,4)*xred + 3.*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)))))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
      Power(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
        0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
        0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
        0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
        0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
        0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
        0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4),2)) - 
   (0.16666666666666666*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 
        5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
        40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
        30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
        30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
        40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
        4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
        10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
        5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
        4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
        1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
        5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
        5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
        15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
        42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
        1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
        65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 2.*Power(t0,5)*Power(te,4)*ve - 
        10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
        10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
        18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
        72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
        72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 
        12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 
        48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
        18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
        1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 
        25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 
        80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 
        1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
        100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 
        4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
        20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
        2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
        2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
        18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
        27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
        2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
        4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
        32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
        3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
        12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred)*
      ((0.05555555555555555*(-t1 + te)*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 
             10.*Power(t0,2)*Power(t1,6)*v0 + 5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 
             30.*Power(t0,3)*Power(t1,4)*te*v0 + 30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 
             3.*Power(t0,5)*t1*Power(te,2)*v0 - 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
             30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 
             1.*Power(t0,5)*Power(te,3)*v0 + 5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
             10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 
             1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 
             10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 
             30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 
             3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 
             30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 
             5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
             5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
             18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
             54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
             36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 
             9.*Power(t1,5)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
             12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
             15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
             6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 
             30.*t0*Power(t1,5)*te*xe - 6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 
             30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 
             3.*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 
             12.*Power(t0,2)*Power(t1,5)*xred + 3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 21.*Power(t0,4)*Power(t1,2)*te*xred - 
             24.*Power(t0,3)*Power(t1,3)*te*xred + 6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 3.*Power(t1,6)*te*xred + 
             3.*Power(t0,5)*Power(te,2)*xred - 6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
             24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 
             12.*Power(t0,3)*t1*Power(te,3)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
        (0.041666666666666664*Power(-t1 + te,2)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 
             10.*Power(t0,2)*Power(t1,7)*v0 - 5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 
             40.*Power(t0,3)*Power(t1,5)*te*v0 - 40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
             60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 
             20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 
             20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 
             5.*Power(t0,4)*t1*Power(te,4)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
             5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
             4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
             1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
             5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
             5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
             15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
             42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
             1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
             65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 
             2.*Power(t0,5)*Power(te,4)*ve - 10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 
             20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 
             12.*Power(t0,3)*Power(t1,5)*x0 + 18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 
             48.*Power(t0,3)*Power(t1,4)*te*x0 - 72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 
             18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 
             72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 48.*t0*Power(t1,4)*Power(te,3)*x0 - 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 
             20.*Power(t0,3)*Power(t1,5)*xe - 25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 
             30.*Power(t0,4)*Power(t1,3)*te*xe - 80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 
             10.*Power(t1,7)*te*xe + 1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
             100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 
             11.*Power(t1,6)*Power(te,2)*xe - 4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 
             40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 
             1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 
             8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 
             18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 
             1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 
             28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 
             7.*Power(t1,6)*Power(te,2)*xred + 4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 
             8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 
             28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 3.*Power(t0,4)*Power(te,4)*xred + 
             12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 
             3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
           (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
             0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
             0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
             0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
             0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
             0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
             0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
        (0.01388888888888889*(1.*Power(t0,5)*Power(t1,4)*v0 - 5.*Power(t0,4)*Power(t1,5)*v0 + 10.*Power(t0,3)*Power(t1,6)*v0 - 10.*Power(t0,2)*Power(t1,7)*v0 + 
             5.*t0*Power(t1,8)*v0 - 1.*Power(t1,9)*v0 - 3.9999999999999996*Power(t0,5)*Power(t1,3)*te*v0 + 20.*Power(t0,4)*Power(t1,4)*te*v0 - 
             40.*Power(t0,3)*Power(t1,5)*te*v0 + 40.*Power(t0,2)*Power(t1,6)*te*v0 - 20.*t0*Power(t1,7)*te*v0 + 3.999999999999999*Power(t1,8)*te*v0 + 
             6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 - 30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 + 59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 - 
             60.00000000000001*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 + 30.000000000000004*t0*Power(t1,6)*Power(te,2)*v0 - 
             6.000000000000001*Power(t1,7)*Power(te,2)*v0 - 3.9999999999999996*Power(t0,5)*t1*Power(te,3)*v0 + 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 - 
             40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 + 40.00000000000001*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 - 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 + 
             4.000000000000001*Power(t1,6)*Power(te,3)*v0 + 1.*Power(t0,5)*Power(te,4)*v0 - 5.*Power(t0,4)*t1*Power(te,4)*v0 + 
             10.*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 - 10.*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 + 5.*t0*Power(t1,4)*Power(te,4)*v0 - 
             1.*Power(t1,5)*Power(te,4)*v0 + 1.*Power(t0,6)*Power(t1,3)*ve - 6.*Power(t0,5)*Power(t1,4)*ve + 14.999999999999998*Power(t0,4)*Power(t1,5)*ve - 
             19.999999999999996*Power(t0,3)*Power(t1,6)*ve + 14.999999999999998*Power(t0,2)*Power(t1,7)*ve - 6.*t0*Power(t1,8)*ve + 1.*Power(t1,9)*ve - 
             3.*Power(t0,6)*Power(t1,2)*te*ve + 18.*Power(t0,5)*Power(t1,3)*te*ve - 44.99999999999999*Power(t0,4)*Power(t1,4)*te*ve + 
             59.99999999999999*Power(t0,3)*Power(t1,5)*te*ve - 44.99999999999999*Power(t0,2)*Power(t1,6)*te*ve + 18.*t0*Power(t1,7)*te*ve - 3.*Power(t1,8)*te*ve + 
             3.*Power(t0,6)*t1*Power(te,2)*ve - 18.*Power(t0,5)*Power(t1,2)*Power(te,2)*ve + 44.99999999999999*Power(t0,4)*Power(t1,3)*Power(te,2)*ve - 
             59.99999999999999*Power(t0,3)*Power(t1,4)*Power(te,2)*ve + 44.99999999999999*Power(t0,2)*Power(t1,5)*Power(te,2)*ve - 18.*t0*Power(t1,6)*Power(te,2)*ve + 
             3.*Power(t1,7)*Power(te,2)*ve - 1.*Power(t0,6)*Power(te,3)*ve + 6.*Power(t0,5)*t1*Power(te,3)*ve - 
             14.999999999999998*Power(t0,4)*Power(t1,2)*Power(te,3)*ve + 19.999999999999996*Power(t0,3)*Power(t1,3)*Power(te,3)*ve - 
             14.999999999999998*Power(t0,2)*Power(t1,4)*Power(te,3)*ve + 6.*t0*Power(t1,5)*Power(te,3)*ve - 1.*Power(t1,6)*Power(te,3)*ve - 
             3.*Power(t0,4)*Power(t1,4)*x0 + 12.*Power(t0,3)*Power(t1,5)*x0 - 18.*Power(t0,2)*Power(t1,6)*x0 + 12.*t0*Power(t1,7)*x0 - 3.*Power(t1,8)*x0 + 
             12.*Power(t0,4)*Power(t1,3)*te*x0 - 48.*Power(t0,3)*Power(t1,4)*te*x0 + 72.*Power(t0,2)*Power(t1,5)*te*x0 - 48.*t0*Power(t1,6)*te*x0 + 
             12.*Power(t1,7)*te*x0 - 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 + 72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 - 
             108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 + 72.*t0*Power(t1,5)*Power(te,2)*x0 - 18.*Power(t1,6)*Power(te,2)*x0 + 12.*Power(t0,4)*t1*Power(te,3)*x0 - 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 + 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 - 48.00000000000001*t0*Power(t1,4)*Power(te,3)*x0 + 
             12.000000000000002*Power(t1,5)*Power(te,3)*x0 - 3.*Power(t0,4)*Power(te,4)*x0 + 12.*Power(t0,3)*t1*Power(te,4)*x0 - 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 + 12.*t0*Power(t1,3)*Power(te,4)*x0 - 3.*Power(t1,4)*Power(te,4)*x0 + 3.*Power(t0,6)*Power(t1,2)*xe - 
             18.*Power(t0,5)*Power(t1,3)*xe + 45.*Power(t0,4)*Power(t1,4)*xe - 60.*Power(t0,3)*Power(t1,5)*xe + 45.*Power(t0,2)*Power(t1,6)*xe - 18.*t0*Power(t1,7)*xe + 
             3.*Power(t1,8)*xe - 6.*Power(t0,6)*t1*te*xe + 36.*Power(t0,5)*Power(t1,2)*te*xe - 90.*Power(t0,4)*Power(t1,3)*te*xe + 120.*Power(t0,3)*Power(t1,4)*te*xe - 
             90.*Power(t0,2)*Power(t1,5)*te*xe + 36.*t0*Power(t1,6)*te*xe - 6.*Power(t1,7)*te*xe + 3.*Power(t0,6)*Power(te,2)*xe - 18.*Power(t0,5)*t1*Power(te,2)*xe + 
             45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe - 60.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe + 45.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe - 
             18.*t0*Power(t1,5)*Power(te,2)*xe + 3.*Power(t1,6)*Power(te,2)*xe - 3.*Power(t0,6)*Power(t1,2)*xred + 18.*Power(t0,5)*Power(t1,3)*xred - 
             42.*Power(t0,4)*Power(t1,4)*xred + 48.*Power(t0,3)*Power(t1,5)*xred - 26.999999999999996*Power(t0,2)*Power(t1,6)*xred + 
             5.999999999999999*t0*Power(t1,7)*xred + 6.*Power(t0,6)*t1*te*xred - 36.*Power(t0,5)*Power(t1,2)*te*xred + 78.*Power(t0,4)*Power(t1,3)*te*xred - 
             72.*Power(t0,3)*Power(t1,4)*te*xred + 18.*Power(t0,2)*Power(t1,5)*te*xred + 11.999999999999998*t0*Power(t1,6)*te*xred - 
             5.999999999999999*Power(t1,7)*te*xred - 3.*Power(t0,6)*Power(te,2)*xred + 18.*Power(t0,5)*t1*Power(te,2)*xred - 
             26.999999999999996*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 12.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred + 63.*Power(t0,2)*Power(t1,4)*Power(te,2)*xred - 
             54.00000000000001*t0*Power(t1,5)*Power(te,2)*xred + 15.*Power(t1,6)*Power(te,2)*xred - 12.*Power(t0,4)*t1*Power(te,3)*xred + 
             48.*Power(t0,3)*Power(t1,2)*Power(te,3)*xred - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred + 48.00000000000001*t0*Power(t1,4)*Power(te,3)*xred - 
             12.000000000000002*Power(t1,5)*Power(te,3)*xred + 3.*Power(t0,4)*Power(te,4)*xred - 12.*Power(t0,3)*t1*Power(te,4)*xred + 
             18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred - 12.*t0*Power(t1,3)*Power(te,4)*xred + 3.*Power(t1,4)*Power(te,4)*xred))/
         ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
             0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
             0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
             0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
             0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
             0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
             0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)))))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,3)*
      (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
        0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
        0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
        0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
        0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
        0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
        0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)));

	return res;
}

static double newton_raphson_cp(double t0, double x0, double v0, double te, double t1, double xred) {
	/* Newton - Raphson Method */
	int iter = 0;
	double h = te_func_cp(t0, x0, v0, te, t1, xred) / dte_func_cp(t0, x0, v0, te, t1, xred);
	while (fabs(h) >= 0.001 && iter < 100)
	{
		h = te_func_cp(t0, x0, v0, te, t1, xred) / dte_func_cp(t0, x0, v0, te, t1, xred);
		  
		te = te - h;
		iter += 1;
	}

	return te;
}

static void UP_problem(double v0, double x0, double t0) {
	int i;
	double step, temp_te, check_te, temp_cost = DBL_MAX;
	double te, t1 = T_max, xred = XMAX;

	/* ~~~~~~~~~~~~~ calculate UP final time ~~~~~~~~~~~~~ */
	double te_init = t0 + 1;
	for (int i = 0; i < newtonIntervals; i++) {
		double te_cost = 0.0;

		check_te = newton_raphson(t0, x0, v0, te_init);
		double checkSpeed = v0;
		int negativeSpeed = 0;
		//if (temp_te > 2.0*max_te || isnan(temp_te) || temp_te < 0.0)
		if (isnan(check_te) || check_te < 0.0)
			continue;

		for (int i = 0; i < P.numsteps; i++) {
			//double step = (te - here.k) / (P.numsteps);
			//double t = (i*(te - here.k) / (P.numsteps - 1) + here.k);
			double step = (check_te - t0) / (P.numsteps);
			double t = (i*(check_te - t0) / (P.numsteps - 1) + t0);

			te_cost += 0.5*pow(UP_control(t, t0, x0, v0, check_te), 2);
			/* skip optimal te that leads to negative speed trajectories */
			checkSpeed += UP_control(t, t0, x0, v0, check_te)*P.T;
			if (checkSpeed < 0)
				negativeSpeed = 1;
		}
		if ((check_te > t0) && (te_cost < temp_cost) && (negativeSpeed != 1)) {
			temp_cost = te_cost;
			temp_te = check_te;
		}
		te_init += 5.0;
	}
	te_up = te = temp_te;
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

	for (i = 0; i < P.numsteps; i++) {
		step = (te - t0) / (P.numsteps);
		double t = (i*(te - t0) / (P.numsteps - 1) + t0);

		initial_u[i] = (6.0*t*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
			(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4)) -
			(2.0*(1.*Power(t0, 4)*v0 - 1.0*Power(t0, 3)*te*v0 - 3.0*Power(t0, 2)*Power(te, 2)*v0 + 5.0*t0*Power(te, 3)*v0 - 2.0*Power(te, 4)*v0 + 2.0*Power(t0, 4)*ve -
				5.0*Power(t0, 3)*te*ve + 3.0*Power(t0, 2)*Power(te, 2)*ve + 1.0*t0*Power(te, 3)*ve - 1.0*Power(te, 4)*ve - 3.0*Power(t0, 3)*x0 + 3.0*Power(t0, 2)*te*x0 +
				3.0*t0*Power(te, 2)*x0 - 3.0*Power(te, 3)*x0 + 3.0*Power(t0, 3)*xe - 3.0*Power(t0, 2)*te*xe - 3.0*t0*Power(te, 2)*xe + 3.0*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4)));

		initial_v[i] = (3.0*Power(t, 2)*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
			(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4)) -
			(2.0*t*(1.*Power(t0, 4)*v0 - 1.0*Power(t0, 3)*te*v0 - 3.0*Power(t0, 2)*Power(te, 2)*v0 + 5.0*t0*Power(te, 3)*v0 - 2.0*Power(te, 4)*v0 + 2.0*Power(t0, 4)*ve -
				5.0*Power(t0, 3)*te*ve + 3.0*Power(t0, 2)*Power(te, 2)*ve + 1.0*t0*Power(te, 3)*ve - 1.0*Power(te, 4)*ve - 3.0*Power(t0, 3)*x0 + 3.0*Power(t0, 2)*te*x0 +
				3.0*t0*Power(te, 2)*x0 - 3.0*Power(te, 3)*x0 + 3.0*Power(t0, 3)*xe - 3.0*Power(t0, 2)*te*xe - 3.0*t0*Power(te, 2)*xe + 3.0*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4))) +
			(1.*(2.0*Power(t0, 4)*te*v0 - 5.0*Power(t0, 3)*Power(te, 2)*v0 + 3.0*Power(t0, 2)*Power(te, 3)*v0 + 1.0*t0*Power(te, 4)*v0 - 1.0*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve -
				1.0*Power(t0, 4)*te*ve - 3.0*Power(t0, 3)*Power(te, 2)*ve + 5.0*Power(t0, 2)*Power(te, 3)*ve - 2.0*t0*Power(te, 4)*ve - 6.0*Power(t0, 3)*te*x0 +
				12.0*Power(t0, 2)*Power(te, 2)*x0 - 6.0*t0*Power(te, 3)*x0 + 6.0*Power(t0, 3)*te*xe - 12.0*Power(t0, 2)*Power(te, 2)*xe + 6.0*t0*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4)));

		initial_x[i] = (1.0*Power(t, 3)*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
			(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4)) -
			(1.0*Power(t, 2)*(1.*Power(t0, 4)*v0 - 1.0*Power(t0, 3)*te*v0 - 3.0*Power(t0, 2)*Power(te, 2)*v0 + 5.0*t0*Power(te, 3)*v0 - 2.0*Power(te, 4)*v0 + 2.0*Power(t0, 4)*ve -
				5.0*Power(t0, 3)*te*ve + 3.0*Power(t0, 2)*Power(te, 2)*ve + 1.0*t0*Power(te, 3)*ve - 1.0*Power(te, 4)*ve - 3.0*Power(t0, 3)*x0 + 3.0*Power(t0, 2)*te*x0 +
				3.0*t0*Power(te, 2)*x0 - 3.0*Power(te, 3)*x0 + 3.0*Power(t0, 3)*xe - 3.0*Power(t0, 2)*te*xe - 3.0*t0*Power(te, 2)*xe + 3.0*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4))) +
			(1.*t*(2.0*Power(t0, 4)*te*v0 - 5.0*Power(t0, 3)*Power(te, 2)*v0 + 3.0*Power(t0, 2)*Power(te, 3)*v0 + 1.0*t0*Power(te, 4)*v0 - 1.0*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve -
				1.0*Power(t0, 4)*te*ve - 3.0*Power(t0, 3)*Power(te, 2)*ve + 5.0*Power(t0, 2)*Power(te, 3)*ve - 2.0*t0*Power(te, 4)*ve - 6.0*Power(t0, 3)*te*x0 +
				12.0*Power(t0, 2)*Power(te, 2)*x0 - 6.0*t0*Power(te, 3)*x0 + 6.0*Power(t0, 3)*te*xe - 12.0*Power(t0, 2)*Power(te, 2)*xe + 6.0*t0*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4))) -
			(1.*(1.*Power(t0, 4)*Power(te, 2)*v0 - 3.0*Power(t0, 3)*Power(te, 3)*v0 + 3.*Power(t0, 2)*Power(te, 4)*v0 - 1.0*t0*Power(te, 5)*v0 + 1.*Power(t0, 5)*te*ve - 3.0*Power(t0, 4)*Power(te, 2)*ve +
				3.*Power(t0, 3)*Power(te, 3)*ve - 1.0*Power(t0, 2)*Power(te, 4)*ve - 3.*Power(t0, 3)*Power(te, 2)*x0 + 7.*Power(t0, 2)*Power(te, 3)*x0 - 5.*t0*Power(te, 4)*x0 + 1.0*Power(te, 5)*x0 - 1.*Power(t0, 5)*xe + 5.*Power(t0, 4)*te*xe -
				7.*Power(t0, 3)*Power(te, 2)*xe + 3.0*Power(t0, 2)*Power(te, 3)*xe)) / ((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 6.0*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0*Power(te, 4)));
	}
}

static double a1(double t, double t0, double x0, double v0, double te, double t1, double xred) {
	double res = -((t*((0.5*Power(t0,2) - 1.*t0*t1 + 0.5*Power(t1,2))*(0. - 1.*(-0.16666666666666666*Power(-1.*t1 + te,3)*(-1.*v0 + 1.*(0. + 1.*ve)) - 
                0.5*Power(-1.*t1 + te,2)*(1.*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te))*v0 + 1.*(0. - 1.*(-1.*x0 + 1.*(0. + 1.*xe)))))) - 
          1.*(0.08333333333333334*Power(-1.*t1 + te,4) - 1.*(-0.16666666666666666*(1.*t0 - 1.*t1)*Power(-1.*t1 + te,3) - 
                0.5*Power(-1.*t1 + te,2)*(-1.*t0*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 1.*(-1.*(0.5*Power(t0,2) - 0.5*Power(t1,2)) + 1.*t1*(-1.*t1 + te)))))*
           (1.*(-1.*t0 + 1.*t1)*v0 + 1.*(1.*x0 - 1.*xred))))/
      (-1.*(-0.5*Power(t0,2)*(-1.*t0 + 1.*t1) + 1.*(-0.16666666666666666*Power(t0,3) + 0.16666666666666666*Power(t1,3)))*
         (0.08333333333333334*Power(-1.*t1 + te,4) - 1.*(-0.16666666666666666*(1.*t0 - 1.*t1)*Power(-1.*t1 + te,3) - 
              0.5*Power(-1.*t1 + te,2)*(-1.*t0*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 1.*(-1.*(0.5*Power(t0,2) - 0.5*Power(t1,2)) + 1.*t1*(-1.*t1 + te))))) + 
        (0.5*Power(t0,2) - 1.*t0*t1 + 0.5*Power(t1,2))*(0.08333333333333334*t1*Power(-1.*t1 + te,4) - 
           1.*(-0.16666666666666666*(0.5*Power(t0,2) - 0.5*Power(t1,2))*Power(-1.*t1 + te,3) - 
              0.5*Power(-1.*t1 + te,2)*(-0.5*Power(t0,2)*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 
                 1.*(-1.*(0.16666666666666666*Power(t0,3) - 0.16666666666666666*Power(t1,3)) + 0.5*Power(t1,2)*(-1.*t1 + te))))))) + 
   (0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 2.*Power(t0,4)*Power(t1,4)*v0 - 0.49999999999999994*Power(t0,3)*Power(t1,5)*v0 + 
        3.5*Power(t0,2)*Power(t1,6)*v0 - 2.5*t0*Power(t1,7)*v0 + 0.49999999999999994*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 
        6.*Power(t0,4)*Power(t1,3)*te*v0 - 6.*Power(t0,2)*Power(t1,5)*te*v0 + 3.*t0*Power(t1,6)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
        6.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 2.9999999999999996*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
        3.0000000000000013*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 6.000000000000003*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 
        1.*Power(t0,5)*Power(te,3)*v0 + 2.*Power(t0,4)*t1*Power(te,3)*v0 - 4.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
        10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 11.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 4.000000000000001*Power(t1,5)*Power(te,3)*v0 + 
        1.5000000000000002*Power(t0,3)*t1*Power(te,4)*v0 - 4.500000000000001*Power(t0,2)*Power(t1,2)*Power(te,4)*v0 + 4.500000000000001*t0*Power(t1,3)*Power(te,4)*v0 - 
        1.5000000000000002*Power(t1,4)*Power(te,4)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 3.5*Power(t0,4)*Power(t1,4)*ve - 4.*Power(t0,3)*Power(t1,5)*ve + 
        1.*Power(t0,2)*Power(t1,6)*ve + 0.9999999999999999*t0*Power(t1,7)*ve - 0.49999999999999994*Power(t1,8)*ve + 3.*Power(t0,5)*Power(t1,2)*te*ve - 
        10.5*Power(t0,4)*Power(t1,3)*te*ve + 12.*Power(t0,3)*Power(t1,4)*te*ve - 3.*Power(t0,2)*Power(t1,5)*te*ve - 2.9999999999999996*t0*Power(t1,6)*te*ve + 
        1.4999999999999998*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 10.5*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
        12.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 3.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve + 2.9999999999999996*t0*Power(t1,5)*Power(te,2)*ve - 
        1.4999999999999998*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 3.5*Power(t0,4)*t1*Power(te,3)*ve + 4.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 
        1.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve - 0.9999999999999999*t0*Power(t1,4)*Power(te,3)*ve + 0.49999999999999994*Power(t1,5)*Power(te,3)*ve - 
        3.*Power(t0,4)*Power(t1,3)*x0 + 6.*Power(t0,3)*Power(t1,4)*x0 - 1.5000000000000002*Power(t0,2)*Power(t1,5)*x0 - 3.*t0*Power(t1,6)*x0 + 1.5*Power(t1,7)*x0 + 
        9.*Power(t0,4)*Power(t1,2)*te*x0 - 18.*Power(t0,3)*Power(t1,3)*te*x0 + 6.*Power(t0,2)*Power(t1,4)*te*x0 + 6.*t0*Power(t1,5)*te*x0 - 3.*Power(t1,6)*te*x0 - 
        9.*Power(t0,4)*t1*Power(te,2)*x0 + 18.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 9.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 
        6.*Power(t0,3)*t1*Power(te,3)*x0 + 6.000000000000001*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 6.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 
        3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 1.5000000000000002*Power(t0,2)*t1*Power(te,4)*x0 + 3.0000000000000004*t0*Power(t1,2)*Power(te,4)*x0 - 
        1.5000000000000002*Power(t1,3)*Power(te,4)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 10.5*Power(t0,4)*Power(t1,3)*xe - 12.*Power(t0,3)*Power(t1,4)*xe + 
        3.*Power(t0,2)*Power(t1,5)*xe + 3.*t0*Power(t1,6)*xe - 1.5*Power(t1,7)*xe + 6.*Power(t0,5)*t1*te*xe - 21.*Power(t0,4)*Power(t1,2)*te*xe + 
        24.*Power(t0,3)*Power(t1,3)*te*xe - 6.*Power(t0,2)*Power(t1,4)*te*xe - 6.*t0*Power(t1,5)*te*xe + 3.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 
        10.5*Power(t0,4)*t1*Power(te,2)*xe - 12.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 3.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe + 3.*t0*Power(t1,4)*Power(te,2)*xe - 
        1.5*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 7.5*Power(t0,4)*Power(t1,3)*xred + 6.*Power(t0,3)*Power(t1,4)*xred - 
        1.5*Power(t0,2)*Power(t1,5)*xred - 6.*Power(t0,5)*t1*te*xred + 12.*Power(t0,4)*Power(t1,2)*te*xred - 6.*Power(t0,3)*Power(t1,3)*te*xred + 
        3.*Power(t0,5)*Power(te,2)*xred - 1.5*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
        6.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 3.0000000000000013*t0*Power(t1,4)*Power(te,2)*xred + 1.5000000000000007*Power(t1,5)*Power(te,2)*xred - 
        3.*Power(t0,4)*Power(te,3)*xred + 6.*Power(t0,3)*t1*Power(te,3)*xred - 6.000000000000001*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 
        6.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred + 1.5000000000000002*Power(t0,2)*t1*Power(te,4)*xred - 
        3.0000000000000004*t0*Power(t1,2)*Power(te,4)*xred + 1.5000000000000002*Power(t1,3)*Power(te,4)*xred))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
        0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
        0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
        0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
        0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
        0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
        0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)));

	return res;
}

static double a2(double t, double t0, double x0, double v0, double te, double t1, double xred) {
	double res = (0.05555555555555555*(1.*Power(t0,5)*Power(t1,3)*v0 - 5.*Power(t0,4)*Power(t1,4)*v0 + 10.*Power(t0,3)*Power(t1,5)*v0 - 10.*Power(t0,2)*Power(t1,6)*v0 + 
        5.*t0*Power(t1,7)*v0 - 1.*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 15.*Power(t0,4)*Power(t1,3)*te*v0 - 30.*Power(t0,3)*Power(t1,4)*te*v0 + 
        30.*Power(t0,2)*Power(t1,5)*te*v0 - 15.*t0*Power(t1,6)*te*v0 + 3.*Power(t1,7)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
        15.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 30.*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 
        15.000000000000002*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 1.*Power(t0,5)*Power(te,3)*v0 + 
        5.*Power(t0,4)*t1*Power(te,3)*v0 - 10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 
        5.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 1.0000000000000004*Power(t1,5)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 
        5.*Power(t0,4)*Power(t1,4)*ve - 10.*Power(t0,3)*Power(t1,5)*ve + 10.*Power(t0,2)*Power(t1,6)*ve - 5.*t0*Power(t1,7)*ve + 1.*Power(t1,8)*ve + 
        3.*Power(t0,5)*Power(t1,2)*te*ve - 15.*Power(t0,4)*Power(t1,3)*te*ve + 30.*Power(t0,3)*Power(t1,4)*te*ve - 30.*Power(t0,2)*Power(t1,5)*te*ve + 
        15.*t0*Power(t1,6)*te*ve - 3.*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 15.*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
        30.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 30.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve - 15.*t0*Power(t1,5)*Power(te,2)*ve + 3.*Power(t1,6)*Power(te,2)*ve + 
        1.*Power(t0,5)*Power(te,3)*ve - 5.*Power(t0,4)*t1*Power(te,3)*ve + 10.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 10.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve + 
        5.*t0*Power(t1,4)*Power(te,3)*ve - 1.*Power(t1,5)*Power(te,3)*ve - 3.*Power(t0,4)*Power(t1,3)*x0 + 12.*Power(t0,3)*Power(t1,4)*x0 - 
        18.*Power(t0,2)*Power(t1,5)*x0 + 12.*t0*Power(t1,6)*x0 - 3.*Power(t1,7)*x0 + 9.*Power(t0,4)*Power(t1,2)*te*x0 - 36.*Power(t0,3)*Power(t1,3)*te*x0 + 
        54.*Power(t0,2)*Power(t1,4)*te*x0 - 36.*t0*Power(t1,5)*te*x0 + 9.*Power(t1,6)*te*x0 - 9.*Power(t0,4)*t1*Power(te,2)*x0 + 
        36.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 54.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 36.*t0*Power(t1,4)*Power(te,2)*x0 - 9.*Power(t1,5)*Power(te,2)*x0 + 
        3.*Power(t0,4)*Power(te,3)*x0 - 12.*Power(t0,3)*t1*Power(te,3)*x0 + 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 
        12.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 
        15.*Power(t0,4)*Power(t1,3)*xe - 30.*Power(t0,3)*Power(t1,4)*xe + 30.*Power(t0,2)*Power(t1,5)*xe - 15.*t0*Power(t1,6)*xe + 3.*Power(t1,7)*xe + 
        6.*Power(t0,5)*t1*te*xe - 30.*Power(t0,4)*Power(t1,2)*te*xe + 60.*Power(t0,3)*Power(t1,3)*te*xe - 60.*Power(t0,2)*Power(t1,4)*te*xe + 30.*t0*Power(t1,5)*te*xe - 
        6.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 15.*Power(t0,4)*t1*Power(te,2)*xe - 30.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 
        30.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe - 15.*t0*Power(t1,4)*Power(te,2)*xe + 3.*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 
        12.*Power(t0,4)*Power(t1,3)*xred + 18.*Power(t0,3)*Power(t1,4)*xred - 12.*Power(t0,2)*Power(t1,5)*xred + 3.*t0*Power(t1,6)*xred - 6.*Power(t0,5)*t1*te*xred + 
        21.*Power(t0,4)*Power(t1,2)*te*xred - 24.*Power(t0,3)*Power(t1,3)*te*xred + 6.*Power(t0,2)*Power(t1,4)*te*xred + 6.*t0*Power(t1,5)*te*xred - 
        3.*Power(t1,6)*te*xred + 3.*Power(t0,5)*Power(te,2)*xred - 6.*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
        24.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 21.*t0*Power(t1,4)*Power(te,2)*xred + 6.*Power(t1,5)*Power(te,2)*xred - 3.*Power(t0,4)*Power(te,3)*xred + 
        12.*Power(t0,3)*t1*Power(te,3)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 12.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 
        3.000000000000001*Power(t1,4)*Power(te,3)*xred))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
        0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
        0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
        0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
        0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
        0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
        0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
   (0.08333333333333333*(t - t1)*(-1.*Power(t0,5)*Power(t1,4)*v0 + 5.*Power(t0,4)*Power(t1,5)*v0 - 10.*Power(t0,3)*Power(t1,6)*v0 + 10.*Power(t0,2)*Power(t1,7)*v0 - 
        5.*t0*Power(t1,8)*v0 + 1.*Power(t1,9)*v0 + 4.*Power(t0,5)*Power(t1,3)*te*v0 - 20.*Power(t0,4)*Power(t1,4)*te*v0 + 40.*Power(t0,3)*Power(t1,5)*te*v0 - 
        40.*Power(t0,2)*Power(t1,6)*te*v0 + 20.*t0*Power(t1,7)*te*v0 - 4.*Power(t1,8)*te*v0 - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
        30.*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 60.*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
        30.*t0*Power(t1,6)*Power(te,2)*v0 + 6.*Power(t1,7)*Power(te,2)*v0 + 4.*Power(t0,5)*t1*Power(te,3)*v0 - 20.*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 + 
        40.*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 - 40.*Power(t0,2)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*t0*Power(t1,5)*Power(te,3)*v0 - 
        4.000000000000001*Power(t1,6)*Power(te,3)*v0 - 1.*Power(t0,5)*Power(te,4)*v0 + 5.*Power(t0,4)*t1*Power(te,4)*v0 - 
        10.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 + 10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
        5.000000000000002*t0*Power(t1,4)*Power(te,4)*v0 + 1.0000000000000007*Power(t1,5)*Power(te,4)*v0 + 1.0000000000000002*Power(t0,6)*Power(t1,3)*ve - 
        4.000000000000001*Power(t0,5)*Power(t1,4)*ve + 5.000000000000001*Power(t0,4)*Power(t1,5)*ve - 5.*Power(t0,2)*Power(t1,7)*ve + 4.*t0*Power(t1,8)*ve - 
        1.*Power(t1,9)*ve - 3.0000000000000004*Power(t0,6)*Power(t1,2)*te*ve + 10.000000000000002*Power(t0,5)*Power(t1,3)*te*ve - 
        5.000000000000001*Power(t0,4)*Power(t1,4)*te*ve - 20.*Power(t0,3)*Power(t1,5)*te*ve + 35.*Power(t0,2)*Power(t1,6)*te*ve - 22.*t0*Power(t1,7)*te*ve + 
        5.*Power(t1,8)*te*ve + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*ve - 6.000000000000001*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 
        15.*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 60.*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 75.*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
        42.*t0*Power(t1,6)*Power(te,2)*ve - 9.*Power(t1,7)*Power(te,2)*ve - 1.0000000000000002*Power(t0,6)*Power(te,3)*ve - 
        1.9999999999999998*Power(t0,5)*t1*Power(te,3)*ve + 25.*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 60.*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
        65.*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 34.*t0*Power(t1,5)*Power(te,3)*ve + 7.000000000000001*Power(t1,6)*Power(te,3)*ve + 2.*Power(t0,5)*Power(te,4)*ve - 
        10.*Power(t0,4)*t1*Power(te,4)*ve + 20.*Power(t0,3)*Power(t1,2)*Power(te,4)*ve - 20.*Power(t0,2)*Power(t1,3)*Power(te,4)*ve + 
        10.*t0*Power(t1,4)*Power(te,4)*ve - 2.*Power(t1,5)*Power(te,4)*ve + 3.*Power(t0,4)*Power(t1,4)*x0 - 12.*Power(t0,3)*Power(t1,5)*x0 + 
        18.*Power(t0,2)*Power(t1,6)*x0 - 12.*t0*Power(t1,7)*x0 + 3.*Power(t1,8)*x0 - 12.*Power(t0,4)*Power(t1,3)*te*x0 + 48.*Power(t0,3)*Power(t1,4)*te*x0 - 
        72.*Power(t0,2)*Power(t1,5)*te*x0 + 48.*t0*Power(t1,6)*te*x0 - 12.*Power(t1,7)*te*x0 + 18.*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 
        72.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 108.*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 72.*t0*Power(t1,5)*Power(te,2)*x0 + 18.*Power(t1,6)*Power(te,2)*x0 - 
        12.000000000000002*Power(t0,4)*t1*Power(te,3)*x0 + 48.*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 - 72.*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 + 
        48.*t0*Power(t1,4)*Power(te,3)*x0 - 12.000000000000002*Power(t1,5)*Power(te,3)*x0 + 3.*Power(t0,4)*Power(te,4)*x0 - 12.*Power(t0,3)*t1*Power(te,4)*x0 + 
        18.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 - 12.000000000000002*t0*Power(t1,3)*Power(te,4)*x0 + 3.000000000000001*Power(t1,4)*Power(te,4)*x0 + 
        1.0000000000000002*Power(t0,6)*Power(t1,2)*xe - 2.0000000000000004*Power(t0,5)*Power(t1,3)*xe - 5.*Power(t0,4)*Power(t1,4)*xe + 20.*Power(t0,3)*Power(t1,5)*xe - 
        25.*Power(t0,2)*Power(t1,6)*xe + 14.*t0*Power(t1,7)*xe - 3.*Power(t1,8)*xe - 2.0000000000000004*Power(t0,6)*t1*te*xe + 30.*Power(t0,4)*Power(t1,3)*te*xe - 
        80.*Power(t0,3)*Power(t1,4)*te*xe + 90.*Power(t0,2)*Power(t1,5)*te*xe - 48.*t0*Power(t1,6)*te*xe + 10.*Power(t1,7)*te*xe + 
        1.0000000000000002*Power(t0,6)*Power(te,2)*xe + 6.*Power(t0,5)*t1*Power(te,2)*xe - 45.*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 
        100.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 105.*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 54.*t0*Power(t1,5)*Power(te,2)*xe - 11.*Power(t1,6)*Power(te,2)*xe - 
        4.*Power(t0,5)*Power(te,3)*xe + 20.*Power(t0,4)*t1*Power(te,3)*xe - 40.*Power(t0,3)*Power(t1,2)*Power(te,3)*xe + 40.*Power(t0,2)*Power(t1,3)*Power(te,3)*xe - 
        20.*t0*Power(t1,4)*Power(te,3)*xe + 4.*Power(t1,5)*Power(te,3)*xe - 1.0000000000000002*Power(t0,6)*Power(t1,2)*xred + 
        2.0000000000000004*Power(t0,5)*Power(t1,3)*xred + 2.*Power(t0,4)*Power(t1,4)*xred - 8.*Power(t0,3)*Power(t1,5)*xred + 7.*Power(t0,2)*Power(t1,6)*xred - 
        2.*t0*Power(t1,7)*xred + 2.0000000000000004*Power(t0,6)*t1*te*xred - 18.*Power(t0,4)*Power(t1,3)*te*xred + 32.*Power(t0,3)*Power(t1,4)*te*xred - 
        18.*Power(t0,2)*Power(t1,5)*te*xred + 2.*Power(t1,7)*te*xred - 1.0000000000000002*Power(t0,6)*Power(te,2)*xred - 6.*Power(t0,5)*t1*Power(te,2)*xred + 
        27.000000000000004*Power(t0,4)*Power(t1,2)*Power(te,2)*xred - 28.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred - 
        2.9999999999999982*Power(t0,2)*Power(t1,4)*Power(te,2)*xred + 18.*t0*Power(t1,5)*Power(te,2)*xred - 7.*Power(t1,6)*Power(te,2)*xred + 
        4.*Power(t0,5)*Power(te,3)*xred - 8.*Power(t0,4)*t1*Power(te,3)*xred - 8.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred + 
        32.*Power(t0,2)*Power(t1,3)*Power(te,3)*xred - 28.000000000000004*t0*Power(t1,4)*Power(te,3)*xred + 8.*Power(t1,5)*Power(te,3)*xred - 
        3.*Power(t0,4)*Power(te,4)*xred + 12.*Power(t0,3)*t1*Power(te,4)*xred - 18.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred + 
        12.000000000000002*t0*Power(t1,3)*Power(te,4)*xred - 3.000000000000001*Power(t1,4)*Power(te,4)*xred))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*Power(1.*t1 - 1.*te,2)*
      (0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 0.08333333333333333*Power(t0,2)*Power(t1,5) + 
        0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 
        0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 
        0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 0.08333333333333334*Power(t1,5)*Power(te,2) + 
        0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 
        0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 
        0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)));

	return res;
}

static double v1(double t, double t0, double x0, double v0, double te, double t1, double xred) {

	double res = (-0.5*Power(t,2)*((0.5*Power(t0,2) - 1.*t0*t1 + 0.5*Power(t1,2))*(0. - 
           1.*(-0.16666666666666666*Power(-1.*t1 + te,3)*(-1.*v0 + 1.*(0. + 1.*ve)) - 
              0.5*Power(-1.*t1 + te,2)*(1.*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te))*v0 + 1.*(0. - 1.*(-1.*x0 + 1.*(0. + 1.*xe)))))) - 
        1.*(0.08333333333333334*Power(-1.*t1 + te,4) - 1.*(-0.16666666666666666*(1.*t0 - 1.*t1)*Power(-1.*t1 + te,3) - 
              0.5*Power(-1.*t1 + te,2)*(-1.*t0*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 1.*(-1.*(0.5*Power(t0,2) - 0.5*Power(t1,2)) + 1.*t1*(-1.*t1 + te)))))*
         (1.*(-1.*t0 + 1.*t1)*v0 + 1.*(1.*x0 - 1.*xred))))/
    (-1.*(-0.5*Power(t0,2)*(-1.*t0 + 1.*t1) + 1.*(-0.16666666666666666*Power(t0,3) + 0.16666666666666666*Power(t1,3)))*
       (0.08333333333333334*Power(-1.*t1 + te,4) - 1.*(-0.16666666666666666*(1.*t0 - 1.*t1)*Power(-1.*t1 + te,3) - 
            0.5*Power(-1.*t1 + te,2)*(-1.*t0*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 1.*(-1.*(0.5*Power(t0,2) - 0.5*Power(t1,2)) + 1.*t1*(-1.*t1 + te))))) + 
      (0.5*Power(t0,2) - 1.*t0*t1 + 0.5*Power(t1,2))*(0.08333333333333334*t1*Power(-1.*t1 + te,4) - 
         1.*(-0.16666666666666666*(0.5*Power(t0,2) - 0.5*Power(t1,2))*Power(-1.*t1 + te,3) - 
            0.5*Power(-1.*t1 + te,2)*(-0.5*Power(t0,2)*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 
               1.*(-1.*(0.16666666666666666*Power(t0,3) - 0.16666666666666666*Power(t1,3)) + 0.5*Power(t1,2)*(-1.*t1 + te)))))) + 
   (0.05555555555555555*t*(1.*Power(t0,5)*Power(t1,3)*v0 - 2.*Power(t0,4)*Power(t1,4)*v0 - 0.49999999999999994*Power(t0,3)*Power(t1,5)*v0 + 
        3.5*Power(t0,2)*Power(t1,6)*v0 - 2.5*t0*Power(t1,7)*v0 + 0.49999999999999994*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 
        6.*Power(t0,4)*Power(t1,3)*te*v0 - 6.*Power(t0,2)*Power(t1,5)*te*v0 + 3.*t0*Power(t1,6)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
        6.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 2.9999999999999996*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
        3.0000000000000013*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 6.000000000000003*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 
        1.*Power(t0,5)*Power(te,3)*v0 + 2.*Power(t0,4)*t1*Power(te,3)*v0 - 4.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
        10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 11.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 4.000000000000001*Power(t1,5)*Power(te,3)*v0 + 
        1.5000000000000002*Power(t0,3)*t1*Power(te,4)*v0 - 4.500000000000001*Power(t0,2)*Power(t1,2)*Power(te,4)*v0 + 4.500000000000001*t0*Power(t1,3)*Power(te,4)*v0 - 
        1.5000000000000002*Power(t1,4)*Power(te,4)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 3.5*Power(t0,4)*Power(t1,4)*ve - 4.*Power(t0,3)*Power(t1,5)*ve + 
        1.*Power(t0,2)*Power(t1,6)*ve + 0.9999999999999999*t0*Power(t1,7)*ve - 0.49999999999999994*Power(t1,8)*ve + 3.*Power(t0,5)*Power(t1,2)*te*ve - 
        10.5*Power(t0,4)*Power(t1,3)*te*ve + 12.*Power(t0,3)*Power(t1,4)*te*ve - 3.*Power(t0,2)*Power(t1,5)*te*ve - 2.9999999999999996*t0*Power(t1,6)*te*ve + 
        1.4999999999999998*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 10.5*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
        12.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 3.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve + 2.9999999999999996*t0*Power(t1,5)*Power(te,2)*ve - 
        1.4999999999999998*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 3.5*Power(t0,4)*t1*Power(te,3)*ve + 4.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 
        1.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve - 0.9999999999999999*t0*Power(t1,4)*Power(te,3)*ve + 0.49999999999999994*Power(t1,5)*Power(te,3)*ve - 
        3.*Power(t0,4)*Power(t1,3)*x0 + 6.*Power(t0,3)*Power(t1,4)*x0 - 1.5000000000000002*Power(t0,2)*Power(t1,5)*x0 - 3.*t0*Power(t1,6)*x0 + 1.5*Power(t1,7)*x0 + 
        9.*Power(t0,4)*Power(t1,2)*te*x0 - 18.*Power(t0,3)*Power(t1,3)*te*x0 + 6.*Power(t0,2)*Power(t1,4)*te*x0 + 6.*t0*Power(t1,5)*te*x0 - 3.*Power(t1,6)*te*x0 - 
        9.*Power(t0,4)*t1*Power(te,2)*x0 + 18.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 9.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 
        6.*Power(t0,3)*t1*Power(te,3)*x0 + 6.000000000000001*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 6.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 
        3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 1.5000000000000002*Power(t0,2)*t1*Power(te,4)*x0 + 3.0000000000000004*t0*Power(t1,2)*Power(te,4)*x0 - 
        1.5000000000000002*Power(t1,3)*Power(te,4)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 10.5*Power(t0,4)*Power(t1,3)*xe - 12.*Power(t0,3)*Power(t1,4)*xe + 
        3.*Power(t0,2)*Power(t1,5)*xe + 3.*t0*Power(t1,6)*xe - 1.5*Power(t1,7)*xe + 6.*Power(t0,5)*t1*te*xe - 21.*Power(t0,4)*Power(t1,2)*te*xe + 
        24.*Power(t0,3)*Power(t1,3)*te*xe - 6.*Power(t0,2)*Power(t1,4)*te*xe - 6.*t0*Power(t1,5)*te*xe + 3.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 
        10.5*Power(t0,4)*t1*Power(te,2)*xe - 12.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 3.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe + 3.*t0*Power(t1,4)*Power(te,2)*xe - 
        1.5*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 7.5*Power(t0,4)*Power(t1,3)*xred + 6.*Power(t0,3)*Power(t1,4)*xred - 
        1.5*Power(t0,2)*Power(t1,5)*xred - 6.*Power(t0,5)*t1*te*xred + 12.*Power(t0,4)*Power(t1,2)*te*xred - 6.*Power(t0,3)*Power(t1,3)*te*xred + 
        3.*Power(t0,5)*Power(te,2)*xred - 1.5*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
        6.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 3.0000000000000013*t0*Power(t1,4)*Power(te,2)*xred + 1.5000000000000007*Power(t1,5)*Power(te,2)*xred - 
        3.*Power(t0,4)*Power(te,3)*xred + 6.*Power(t0,3)*t1*Power(te,3)*xred - 6.000000000000001*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 
        6.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred + 1.5000000000000002*Power(t0,2)*t1*Power(te,4)*xred - 
        3.0000000000000004*t0*Power(t1,2)*Power(te,4)*xred + 1.5000000000000002*Power(t1,3)*Power(te,4)*xred))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
        0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
        0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
        0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
        0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
        0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
        0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
   (0.01388888888888889*(-2.9999999999999996*Power(t0,5)*Power(t1,4)*v0 + 8.999999999999998*Power(t0,4)*Power(t1,5)*v0 - 8.999999999999998*Power(t0,3)*Power(t1,6)*v0 + 
        2.9999999999999996*Power(t0,2)*Power(t1,7)*v0 + 7.999999999999999*Power(t0,5)*Power(t1,3)*te*v0 - 21.999999999999996*Power(t0,4)*Power(t1,4)*te*v0 + 
        19.999999999999996*Power(t0,3)*Power(t1,5)*te*v0 - 7.999999999999999*Power(t0,2)*Power(t1,6)*te*v0 + 3.9999999999999996*t0*Power(t1,7)*te*v0 - 
        1.9999999999999996*Power(t1,8)*te*v0 - 5.999999999999999*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 11.999999999999998*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 
        5.9999999999999964*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 5.999999999999997*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 
        12.000000000000002*t0*Power(t1,6)*Power(te,2)*v0 + 6.000000000000001*Power(t1,7)*Power(te,2)*v0 + 6.000000000000001*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 - 
        12.000000000000002*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 + 12.000000000000002*t0*Power(t1,5)*Power(te,3)*v0 - 6.000000000000001*Power(t1,6)*Power(te,3)*v0 + 
        1.*Power(t0,5)*Power(te,4)*v0 - 5.*Power(t0,4)*t1*Power(te,4)*v0 + 7.*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 - 1.*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 
        4.*t0*Power(t1,4)*Power(te,4)*v0 + 2.*Power(t1,5)*Power(te,4)*v0 + 1.*Power(t0,6)*Power(t1,3)*ve - 2.*Power(t0,5)*Power(t1,4)*ve - 
        1.9999999999999998*Power(t0,4)*Power(t1,5)*ve + 7.999999999999999*Power(t0,3)*Power(t1,6)*ve - 6.999999999999999*Power(t0,2)*Power(t1,7)*ve + 
        1.9999999999999998*t0*Power(t1,8)*ve - 3.*Power(t0,6)*Power(t1,2)*te*ve + 6.*Power(t0,5)*Power(t1,3)*te*ve + 5.999999999999999*Power(t0,4)*Power(t1,4)*te*ve - 
        23.999999999999996*Power(t0,3)*Power(t1,5)*te*ve + 20.999999999999996*Power(t0,2)*Power(t1,6)*te*ve - 5.999999999999999*t0*Power(t1,7)*te*ve + 
        3.*Power(t0,6)*t1*Power(te,2)*ve - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 5.999999999999999*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 
        23.999999999999996*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 20.999999999999996*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
        5.999999999999999*t0*Power(t1,6)*Power(te,2)*ve - 1.*Power(t0,6)*Power(te,3)*ve + 2.*Power(t0,5)*t1*Power(te,3)*ve + 
        1.9999999999999998*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 7.999999999999999*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
        6.999999999999999*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 1.9999999999999998*t0*Power(t1,5)*Power(te,3)*ve + 9.*Power(t0,4)*Power(t1,4)*x0 - 
        24.*Power(t0,3)*Power(t1,5)*x0 + 20.999999999999996*Power(t0,2)*Power(t1,6)*x0 - 5.999999999999999*t0*Power(t1,7)*x0 - 24.*Power(t0,4)*Power(t1,3)*te*x0 + 
        59.99999999999999*Power(t0,3)*Power(t1,4)*te*x0 - 47.99999999999999*Power(t0,2)*Power(t1,5)*te*x0 + 11.999999999999998*t0*Power(t1,6)*te*x0 + 
        17.999999999999996*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 36.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 
        17.999999999999993*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 12.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 + 
        24.000000000000004*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 - 12.000000000000002*t0*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,4)*Power(te,4)*x0 + 
        12.*Power(t0,3)*t1*Power(te,4)*x0 - 15.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 + 6.*t0*Power(t1,3)*Power(te,4)*x0 + 3.*Power(t0,6)*Power(t1,2)*xe - 
        6.*Power(t0,5)*Power(t1,3)*xe - 5.999999999999999*Power(t0,4)*Power(t1,4)*xe + 24.*Power(t0,3)*Power(t1,5)*xe - 20.999999999999996*Power(t0,2)*Power(t1,6)*xe + 
        5.999999999999999*t0*Power(t1,7)*xe - 6.*Power(t0,6)*t1*te*xe + 12.*Power(t0,5)*Power(t1,2)*te*xe + 11.999999999999998*Power(t0,4)*Power(t1,3)*te*xe - 
        48.*Power(t0,3)*Power(t1,4)*te*xe + 41.99999999999999*Power(t0,2)*Power(t1,5)*te*xe - 11.999999999999998*t0*Power(t1,6)*te*xe + 3.*Power(t0,6)*Power(te,2)*xe - 
        6.*Power(t0,5)*t1*Power(te,2)*xe - 5.999999999999999*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 24.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
        20.999999999999996*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 5.999999999999999*t0*Power(t1,5)*Power(te,2)*xe - 3.*Power(t0,6)*Power(t1,2)*xred + 
        6.*Power(t0,5)*Power(t1,3)*xred - 3.*Power(t0,4)*Power(t1,4)*xred + 6.*Power(t0,6)*t1*te*xred - 12.*Power(t0,5)*Power(t1,2)*te*xred + 
        12.*Power(t0,4)*Power(t1,3)*te*xred - 11.999999999999998*Power(t0,3)*Power(t1,4)*te*xred + 5.999999999999999*Power(t0,2)*Power(t1,5)*te*xred - 
        3.*Power(t0,6)*Power(te,2)*xred + 6.*Power(t0,5)*t1*Power(te,2)*xred - 11.999999999999998*Power(t0,4)*Power(t1,2)*Power(te,2)*xred + 
        12.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred + 3.000000000000006*Power(t0,2)*Power(t1,4)*Power(te,2)*xred - 6.000000000000003*t0*Power(t1,5)*Power(te,2)*xred + 
        12.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred - 24.000000000000004*Power(t0,2)*Power(t1,3)*Power(te,3)*xred + 
        12.000000000000002*t0*Power(t1,4)*Power(te,3)*xred + 3.*Power(t0,4)*Power(te,4)*xred - 12.*Power(t0,3)*t1*Power(te,4)*xred + 
        15.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred - 6.*t0*Power(t1,3)*Power(te,4)*xred))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
        0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
        0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
        0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
        0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
        0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
        0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)));

	return res;
}

static double x1(double t, double t0, double x0, double v0, double te, double t1, double xred) {

	double res = (-0.16666666666666666*Power(t,3)*((0.5*Power(t0,2) - 1.*t0*t1 + 0.5*Power(t1,2))*
         (0. - 1.*(-0.16666666666666666*Power(-1.*t1 + te,3)*(-1.*v0 + 1.*(0. + 1.*ve)) - 
              0.5*Power(-1.*t1 + te,2)*(1.*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te))*v0 + 1.*(0. - 1.*(-1.*x0 + 1.*(0. + 1.*xe)))))) - 
        1.*(0.08333333333333334*Power(-1.*t1 + te,4) - 1.*(-0.16666666666666666*(1.*t0 - 1.*t1)*Power(-1.*t1 + te,3) - 
              0.5*Power(-1.*t1 + te,2)*(-1.*t0*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 1.*(-1.*(0.5*Power(t0,2) - 0.5*Power(t1,2)) + 1.*t1*(-1.*t1 + te)))))*
         (1.*(-1.*t0 + 1.*t1)*v0 + 1.*(1.*x0 - 1.*xred))))/
    (-1.*(-0.5*Power(t0,2)*(-1.*t0 + 1.*t1) + 1.*(-0.16666666666666666*Power(t0,3) + 0.16666666666666666*Power(t1,3)))*
       (0.08333333333333334*Power(-1.*t1 + te,4) - 1.*(-0.16666666666666666*(1.*t0 - 1.*t1)*Power(-1.*t1 + te,3) - 
            0.5*Power(-1.*t1 + te,2)*(-1.*t0*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 1.*(-1.*(0.5*Power(t0,2) - 0.5*Power(t1,2)) + 1.*t1*(-1.*t1 + te))))) + 
      (0.5*Power(t0,2) - 1.*t0*t1 + 0.5*Power(t1,2))*(0.08333333333333334*t1*Power(-1.*t1 + te,4) - 
         1.*(-0.16666666666666666*(0.5*Power(t0,2) - 0.5*Power(t1,2))*Power(-1.*t1 + te,3) - 
            0.5*Power(-1.*t1 + te,2)*(-0.5*Power(t0,2)*(-1.*(1.*t0 - 1.*t1) + 1.*(-1.*t1 + te)) + 
               1.*(-1.*(0.16666666666666666*Power(t0,3) - 0.16666666666666666*Power(t1,3)) + 0.5*Power(t1,2)*(-1.*t1 + te)))))) + 
   (0.027777777777777776*Power(t,2)*(1.*Power(t0,5)*Power(t1,3)*v0 - 2.*Power(t0,4)*Power(t1,4)*v0 - 0.49999999999999994*Power(t0,3)*Power(t1,5)*v0 + 
        3.5*Power(t0,2)*Power(t1,6)*v0 - 2.5*t0*Power(t1,7)*v0 + 0.49999999999999994*Power(t1,8)*v0 - 3.*Power(t0,5)*Power(t1,2)*te*v0 + 
        6.*Power(t0,4)*Power(t1,3)*te*v0 - 6.*Power(t0,2)*Power(t1,5)*te*v0 + 3.*t0*Power(t1,6)*te*v0 + 3.*Power(t0,5)*t1*Power(te,2)*v0 - 
        6.*Power(t0,4)*Power(t1,2)*Power(te,2)*v0 + 2.9999999999999996*Power(t0,3)*Power(t1,3)*Power(te,2)*v0 - 
        3.0000000000000013*Power(t0,2)*Power(t1,4)*Power(te,2)*v0 + 6.000000000000003*t0*Power(t1,5)*Power(te,2)*v0 - 3.0000000000000004*Power(t1,6)*Power(te,2)*v0 - 
        1.*Power(t0,5)*Power(te,3)*v0 + 2.*Power(t0,4)*t1*Power(te,3)*v0 - 4.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,3)*v0 + 
        10.000000000000002*Power(t0,2)*Power(t1,3)*Power(te,3)*v0 - 11.000000000000002*t0*Power(t1,4)*Power(te,3)*v0 + 4.000000000000001*Power(t1,5)*Power(te,3)*v0 + 
        1.5000000000000002*Power(t0,3)*t1*Power(te,4)*v0 - 4.500000000000001*Power(t0,2)*Power(t1,2)*Power(te,4)*v0 + 4.500000000000001*t0*Power(t1,3)*Power(te,4)*v0 - 
        1.5000000000000002*Power(t1,4)*Power(te,4)*v0 - 1.*Power(t0,5)*Power(t1,3)*ve + 3.5*Power(t0,4)*Power(t1,4)*ve - 4.*Power(t0,3)*Power(t1,5)*ve + 
        1.*Power(t0,2)*Power(t1,6)*ve + 0.9999999999999999*t0*Power(t1,7)*ve - 0.49999999999999994*Power(t1,8)*ve + 3.*Power(t0,5)*Power(t1,2)*te*ve - 
        10.5*Power(t0,4)*Power(t1,3)*te*ve + 12.*Power(t0,3)*Power(t1,4)*te*ve - 3.*Power(t0,2)*Power(t1,5)*te*ve - 2.9999999999999996*t0*Power(t1,6)*te*ve + 
        1.4999999999999998*Power(t1,7)*te*ve - 3.*Power(t0,5)*t1*Power(te,2)*ve + 10.5*Power(t0,4)*Power(t1,2)*Power(te,2)*ve - 
        12.*Power(t0,3)*Power(t1,3)*Power(te,2)*ve + 3.*Power(t0,2)*Power(t1,4)*Power(te,2)*ve + 2.9999999999999996*t0*Power(t1,5)*Power(te,2)*ve - 
        1.4999999999999998*Power(t1,6)*Power(te,2)*ve + 1.*Power(t0,5)*Power(te,3)*ve - 3.5*Power(t0,4)*t1*Power(te,3)*ve + 4.*Power(t0,3)*Power(t1,2)*Power(te,3)*ve - 
        1.*Power(t0,2)*Power(t1,3)*Power(te,3)*ve - 0.9999999999999999*t0*Power(t1,4)*Power(te,3)*ve + 0.49999999999999994*Power(t1,5)*Power(te,3)*ve - 
        3.*Power(t0,4)*Power(t1,3)*x0 + 6.*Power(t0,3)*Power(t1,4)*x0 - 1.5000000000000002*Power(t0,2)*Power(t1,5)*x0 - 3.*t0*Power(t1,6)*x0 + 1.5*Power(t1,7)*x0 + 
        9.*Power(t0,4)*Power(t1,2)*te*x0 - 18.*Power(t0,3)*Power(t1,3)*te*x0 + 6.*Power(t0,2)*Power(t1,4)*te*x0 + 6.*t0*Power(t1,5)*te*x0 - 3.*Power(t1,6)*te*x0 - 
        9.*Power(t0,4)*t1*Power(te,2)*x0 + 18.*Power(t0,3)*Power(t1,2)*Power(te,2)*x0 - 9.*Power(t0,2)*Power(t1,3)*Power(te,2)*x0 + 3.*Power(t0,4)*Power(te,3)*x0 - 
        6.*Power(t0,3)*t1*Power(te,3)*x0 + 6.000000000000001*Power(t0,2)*Power(t1,2)*Power(te,3)*x0 - 6.000000000000002*t0*Power(t1,3)*Power(te,3)*x0 + 
        3.000000000000001*Power(t1,4)*Power(te,3)*x0 - 1.5000000000000002*Power(t0,2)*t1*Power(te,4)*x0 + 3.0000000000000004*t0*Power(t1,2)*Power(te,4)*x0 - 
        1.5000000000000002*Power(t1,3)*Power(te,4)*x0 - 3.*Power(t0,5)*Power(t1,2)*xe + 10.5*Power(t0,4)*Power(t1,3)*xe - 12.*Power(t0,3)*Power(t1,4)*xe + 
        3.*Power(t0,2)*Power(t1,5)*xe + 3.*t0*Power(t1,6)*xe - 1.5*Power(t1,7)*xe + 6.*Power(t0,5)*t1*te*xe - 21.*Power(t0,4)*Power(t1,2)*te*xe + 
        24.*Power(t0,3)*Power(t1,3)*te*xe - 6.*Power(t0,2)*Power(t1,4)*te*xe - 6.*t0*Power(t1,5)*te*xe + 3.*Power(t1,6)*te*xe - 3.*Power(t0,5)*Power(te,2)*xe + 
        10.5*Power(t0,4)*t1*Power(te,2)*xe - 12.*Power(t0,3)*Power(t1,2)*Power(te,2)*xe + 3.*Power(t0,2)*Power(t1,3)*Power(te,2)*xe + 3.*t0*Power(t1,4)*Power(te,2)*xe - 
        1.5*Power(t1,5)*Power(te,2)*xe + 3.*Power(t0,5)*Power(t1,2)*xred - 7.5*Power(t0,4)*Power(t1,3)*xred + 6.*Power(t0,3)*Power(t1,4)*xred - 
        1.5*Power(t0,2)*Power(t1,5)*xred - 6.*Power(t0,5)*t1*te*xred + 12.*Power(t0,4)*Power(t1,2)*te*xred - 6.*Power(t0,3)*Power(t1,3)*te*xred + 
        3.*Power(t0,5)*Power(te,2)*xred - 1.5*Power(t0,4)*t1*Power(te,2)*xred - 6.000000000000001*Power(t0,3)*Power(t1,2)*Power(te,2)*xred + 
        6.*Power(t0,2)*Power(t1,3)*Power(te,2)*xred - 3.0000000000000013*t0*Power(t1,4)*Power(te,2)*xred + 1.5000000000000007*Power(t1,5)*Power(te,2)*xred - 
        3.*Power(t0,4)*Power(te,3)*xred + 6.*Power(t0,3)*t1*Power(te,3)*xred - 6.000000000000001*Power(t0,2)*Power(t1,2)*Power(te,3)*xred + 
        6.000000000000002*t0*Power(t1,3)*Power(te,3)*xred - 3.000000000000001*Power(t1,4)*Power(te,3)*xred + 1.5000000000000002*Power(t0,2)*t1*Power(te,4)*xred - 
        3.0000000000000004*t0*Power(t1,2)*Power(te,4)*xred + 1.5000000000000002*Power(t1,3)*Power(te,4)*xred))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
        0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
        0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
        0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
        0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
        0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
        0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) + 
   (0.01388888888888889*t*(-2.9999999999999996*Power(t0,5)*Power(t1,4)*v0 + 8.999999999999998*Power(t0,4)*Power(t1,5)*v0 - 
        8.999999999999998*Power(t0,3)*Power(t1,6)*v0 + 2.9999999999999996*Power(t0,2)*Power(t1,7)*v0 + 7.999999999999999*Power(t0,5)*Power(t1,3)*te*v0 - 
        21.999999999999996*Power(t0,4)*Power(t1,4)*te*v0 + 19.999999999999996*Power(t0,3)*Power(t1,5)*te*v0 - 7.999999999999999*Power(t0,2)*Power(t1,6)*te*v0 + 
        3.9999999999999996*t0*Power(t1,7)*te*v0 - 1.9999999999999996*Power(t1,8)*te*v0 - 5.999999999999999*Power(t0,5)*Power(t1,2)*Power(te,2)*v0 + 
        11.999999999999998*Power(t0,4)*Power(t1,3)*Power(te,2)*v0 - 5.9999999999999964*Power(t0,3)*Power(t1,4)*Power(te,2)*v0 + 
        5.999999999999997*Power(t0,2)*Power(t1,5)*Power(te,2)*v0 - 12.000000000000002*t0*Power(t1,6)*Power(te,2)*v0 + 6.000000000000001*Power(t1,7)*Power(te,2)*v0 + 
        6.000000000000001*Power(t0,4)*Power(t1,2)*Power(te,3)*v0 - 12.000000000000002*Power(t0,3)*Power(t1,3)*Power(te,3)*v0 + 
        12.000000000000002*t0*Power(t1,5)*Power(te,3)*v0 - 6.000000000000001*Power(t1,6)*Power(te,3)*v0 + 1.*Power(t0,5)*Power(te,4)*v0 - 
        5.*Power(t0,4)*t1*Power(te,4)*v0 + 7.*Power(t0,3)*Power(t1,2)*Power(te,4)*v0 - 1.*Power(t0,2)*Power(t1,3)*Power(te,4)*v0 - 4.*t0*Power(t1,4)*Power(te,4)*v0 + 
        2.*Power(t1,5)*Power(te,4)*v0 + 1.*Power(t0,6)*Power(t1,3)*ve - 2.*Power(t0,5)*Power(t1,4)*ve - 1.9999999999999998*Power(t0,4)*Power(t1,5)*ve + 
        7.999999999999999*Power(t0,3)*Power(t1,6)*ve - 6.999999999999999*Power(t0,2)*Power(t1,7)*ve + 1.9999999999999998*t0*Power(t1,8)*ve - 
        3.*Power(t0,6)*Power(t1,2)*te*ve + 6.*Power(t0,5)*Power(t1,3)*te*ve + 5.999999999999999*Power(t0,4)*Power(t1,4)*te*ve - 
        23.999999999999996*Power(t0,3)*Power(t1,5)*te*ve + 20.999999999999996*Power(t0,2)*Power(t1,6)*te*ve - 5.999999999999999*t0*Power(t1,7)*te*ve + 
        3.*Power(t0,6)*t1*Power(te,2)*ve - 6.*Power(t0,5)*Power(t1,2)*Power(te,2)*ve - 5.999999999999999*Power(t0,4)*Power(t1,3)*Power(te,2)*ve + 
        23.999999999999996*Power(t0,3)*Power(t1,4)*Power(te,2)*ve - 20.999999999999996*Power(t0,2)*Power(t1,5)*Power(te,2)*ve + 
        5.999999999999999*t0*Power(t1,6)*Power(te,2)*ve - 1.*Power(t0,6)*Power(te,3)*ve + 2.*Power(t0,5)*t1*Power(te,3)*ve + 
        1.9999999999999998*Power(t0,4)*Power(t1,2)*Power(te,3)*ve - 7.999999999999999*Power(t0,3)*Power(t1,3)*Power(te,3)*ve + 
        6.999999999999999*Power(t0,2)*Power(t1,4)*Power(te,3)*ve - 1.9999999999999998*t0*Power(t1,5)*Power(te,3)*ve + 9.*Power(t0,4)*Power(t1,4)*x0 - 
        24.*Power(t0,3)*Power(t1,5)*x0 + 20.999999999999996*Power(t0,2)*Power(t1,6)*x0 - 5.999999999999999*t0*Power(t1,7)*x0 - 24.*Power(t0,4)*Power(t1,3)*te*x0 + 
        59.99999999999999*Power(t0,3)*Power(t1,4)*te*x0 - 47.99999999999999*Power(t0,2)*Power(t1,5)*te*x0 + 11.999999999999998*t0*Power(t1,6)*te*x0 + 
        17.999999999999996*Power(t0,4)*Power(t1,2)*Power(te,2)*x0 - 36.*Power(t0,3)*Power(t1,3)*Power(te,2)*x0 + 
        17.999999999999993*Power(t0,2)*Power(t1,4)*Power(te,2)*x0 - 12.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*x0 + 
        24.000000000000004*Power(t0,2)*Power(t1,3)*Power(te,3)*x0 - 12.000000000000002*t0*Power(t1,4)*Power(te,3)*x0 - 3.*Power(t0,4)*Power(te,4)*x0 + 
        12.*Power(t0,3)*t1*Power(te,4)*x0 - 15.*Power(t0,2)*Power(t1,2)*Power(te,4)*x0 + 6.*t0*Power(t1,3)*Power(te,4)*x0 + 3.*Power(t0,6)*Power(t1,2)*xe - 
        6.*Power(t0,5)*Power(t1,3)*xe - 5.999999999999999*Power(t0,4)*Power(t1,4)*xe + 24.*Power(t0,3)*Power(t1,5)*xe - 20.999999999999996*Power(t0,2)*Power(t1,6)*xe + 
        5.999999999999999*t0*Power(t1,7)*xe - 6.*Power(t0,6)*t1*te*xe + 12.*Power(t0,5)*Power(t1,2)*te*xe + 11.999999999999998*Power(t0,4)*Power(t1,3)*te*xe - 
        48.*Power(t0,3)*Power(t1,4)*te*xe + 41.99999999999999*Power(t0,2)*Power(t1,5)*te*xe - 11.999999999999998*t0*Power(t1,6)*te*xe + 3.*Power(t0,6)*Power(te,2)*xe - 
        6.*Power(t0,5)*t1*Power(te,2)*xe - 5.999999999999999*Power(t0,4)*Power(t1,2)*Power(te,2)*xe + 24.*Power(t0,3)*Power(t1,3)*Power(te,2)*xe - 
        20.999999999999996*Power(t0,2)*Power(t1,4)*Power(te,2)*xe + 5.999999999999999*t0*Power(t1,5)*Power(te,2)*xe - 3.*Power(t0,6)*Power(t1,2)*xred + 
        6.*Power(t0,5)*Power(t1,3)*xred - 3.*Power(t0,4)*Power(t1,4)*xred + 6.*Power(t0,6)*t1*te*xred - 12.*Power(t0,5)*Power(t1,2)*te*xred + 
        12.*Power(t0,4)*Power(t1,3)*te*xred - 11.999999999999998*Power(t0,3)*Power(t1,4)*te*xred + 5.999999999999999*Power(t0,2)*Power(t1,5)*te*xred - 
        3.*Power(t0,6)*Power(te,2)*xred + 6.*Power(t0,5)*t1*Power(te,2)*xred - 11.999999999999998*Power(t0,4)*Power(t1,2)*Power(te,2)*xred + 
        12.*Power(t0,3)*Power(t1,3)*Power(te,2)*xred + 3.000000000000006*Power(t0,2)*Power(t1,4)*Power(te,2)*xred - 6.000000000000003*t0*Power(t1,5)*Power(te,2)*xred + 
        12.000000000000002*Power(t0,3)*Power(t1,2)*Power(te,3)*xred - 24.000000000000004*Power(t0,2)*Power(t1,3)*Power(te,3)*xred + 
        12.000000000000002*t0*Power(t1,4)*Power(te,3)*xred + 3.*Power(t0,4)*Power(te,4)*xred - 12.*Power(t0,3)*t1*Power(te,4)*xred + 
        15.*Power(t0,2)*Power(t1,2)*Power(te,4)*xred - 6.*t0*Power(t1,3)*Power(te,4)*xred))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
        0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
        0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
        0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
        0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
        0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
        0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4))) - 
   (0.01388888888888889*(-0.9999999999999999*Power(t0,5)*Power(t1,5)*v0 + 2.9999999999999996*Power(t0,4)*Power(t1,6)*v0 - 
        2.9999999999999996*Power(t0,3)*Power(t1,7)*v0 + 0.9999999999999998*Power(t0,2)*Power(t1,8)*v0 + 1.9999999999999996*Power(t0,5)*Power(t1,4)*te*v0 - 
        3.999999999999999*Power(t0,4)*Power(t1,5)*te*v0 + 3.999999999999999*Power(t0,2)*Power(t1,7)*te*v0 - 1.9999999999999996*t0*Power(t1,8)*te*v0 - 
        5.9999999999999964*Power(t0,4)*Power(t1,4)*Power(te,2)*v0 + 18.*Power(t0,3)*Power(t1,5)*Power(te,2)*v0 - 
        18.000000000000004*Power(t0,2)*Power(t1,6)*Power(te,2)*v0 + 6.000000000000001*t0*Power(t1,7)*Power(te,2)*v0 - 
        1.9999999999999998*Power(t0,5)*Power(t1,2)*Power(te,3)*v0 + 12.*Power(t0,4)*Power(t1,3)*Power(te,3)*v0 - 
        24.000000000000004*Power(t0,3)*Power(t1,4)*Power(te,3)*v0 + 20.000000000000004*Power(t0,2)*Power(t1,5)*Power(te,3)*v0 - 
        6.000000000000001*t0*Power(t1,6)*Power(te,3)*v0 + 1.*Power(t0,5)*t1*Power(te,4)*v0 - 5.*Power(t0,4)*Power(t1,2)*Power(te,4)*v0 + 
        9.*Power(t0,3)*Power(t1,3)*Power(te,4)*v0 - 7.000000000000001*Power(t0,2)*Power(t1,4)*Power(te,4)*v0 + 2.*t0*Power(t1,5)*Power(te,4)*v0 + 
        1.*Power(t0,6)*Power(t1,4)*ve - 4.*Power(t0,5)*Power(t1,5)*ve + 5.999999999999999*Power(t0,4)*Power(t1,6)*ve - 3.9999999999999996*Power(t0,3)*Power(t1,7)*ve + 
        0.9999999999999999*Power(t0,2)*Power(t1,8)*ve - 3.*Power(t0,6)*Power(t1,3)*te*ve + 12.*Power(t0,5)*Power(t1,4)*te*ve - 18.*Power(t0,4)*Power(t1,5)*te*ve + 
        11.999999999999998*Power(t0,3)*Power(t1,6)*te*ve - 2.9999999999999996*Power(t0,2)*Power(t1,7)*te*ve + 3.*Power(t0,6)*Power(t1,2)*Power(te,2)*ve - 
        12.*Power(t0,5)*Power(t1,3)*Power(te,2)*ve + 18.*Power(t0,4)*Power(t1,4)*Power(te,2)*ve - 11.999999999999998*Power(t0,3)*Power(t1,5)*Power(te,2)*ve + 
        2.9999999999999996*Power(t0,2)*Power(t1,6)*Power(te,2)*ve - 1.*Power(t0,6)*t1*Power(te,3)*ve + 4.*Power(t0,5)*Power(t1,2)*Power(te,3)*ve - 
        5.999999999999999*Power(t0,4)*Power(t1,3)*Power(te,3)*ve + 3.9999999999999996*Power(t0,3)*Power(t1,4)*Power(te,3)*ve - 
        0.9999999999999999*Power(t0,2)*Power(t1,5)*Power(te,3)*ve + 2.9999999999999996*Power(t0,4)*Power(t1,5)*x0 - 7.999999999999999*Power(t0,3)*Power(t1,6)*x0 + 
        6.999999999999999*Power(t0,2)*Power(t1,7)*x0 - 1.9999999999999996*t0*Power(t1,8)*x0 - 5.999999999999999*Power(t0,4)*Power(t1,4)*te*x0 + 
        11.999999999999998*Power(t0,3)*Power(t1,5)*te*x0 - 3.9999999999999996*Power(t0,2)*Power(t1,6)*te*x0 - 3.999999999999999*t0*Power(t1,7)*te*x0 + 
        1.9999999999999996*Power(t1,8)*te*x0 + 11.999999999999996*Power(t0,3)*Power(t1,4)*Power(te,2)*x0 - 30.000000000000004*Power(t0,2)*Power(t1,5)*Power(te,2)*x0 + 
        24.000000000000004*t0*Power(t1,6)*Power(te,2)*x0 - 6.000000000000001*Power(t1,7)*Power(te,2)*x0 + 6.*Power(t0,4)*Power(t1,2)*Power(te,3)*x0 - 
        28.*Power(t0,3)*Power(t1,3)*Power(te,3)*x0 + 44.*Power(t0,2)*Power(t1,4)*Power(te,3)*x0 - 28.000000000000004*t0*Power(t1,5)*Power(te,3)*x0 + 
        6.000000000000001*Power(t1,6)*Power(te,3)*x0 - 3.*Power(t0,4)*t1*Power(te,4)*x0 + 12.*Power(t0,3)*Power(t1,2)*Power(te,4)*x0 - 
        17.*Power(t0,2)*Power(t1,3)*Power(te,4)*x0 + 10.*t0*Power(t1,4)*Power(te,4)*x0 - 2.*Power(t1,5)*Power(te,4)*x0 + 3.*Power(t0,6)*Power(t1,3)*xe - 
        12.*Power(t0,5)*Power(t1,4)*xe + 18.*Power(t0,4)*Power(t1,5)*xe - 12.*Power(t0,3)*Power(t1,6)*xe + 2.9999999999999996*Power(t0,2)*Power(t1,7)*xe - 
        6.*Power(t0,6)*Power(t1,2)*te*xe + 24.*Power(t0,5)*Power(t1,3)*te*xe - 36.*Power(t0,4)*Power(t1,4)*te*xe + 24.*Power(t0,3)*Power(t1,5)*te*xe - 
        5.999999999999999*Power(t0,2)*Power(t1,6)*te*xe + 3.*Power(t0,6)*t1*Power(te,2)*xe - 12.*Power(t0,5)*Power(t1,2)*Power(te,2)*xe + 
        18.*Power(t0,4)*Power(t1,3)*Power(te,2)*xe - 12.*Power(t0,3)*Power(t1,4)*Power(te,2)*xe + 2.9999999999999996*Power(t0,2)*Power(t1,5)*Power(te,2)*xe - 
        1.*Power(t0,6)*Power(t1,3)*xred + 2.*Power(t0,5)*Power(t1,4)*xred - 1.*Power(t0,4)*Power(t1,5)*xred + 3.9999999999999996*Power(t0,5)*Power(t1,3)*te*xred - 
        7.999999999999999*Power(t0,4)*Power(t1,4)*te*xred + 3.9999999999999996*Power(t0,3)*Power(t1,5)*te*xred + 3.0000000000000004*Power(t0,6)*t1*Power(te,2)*xred - 
        11.999999999999998*Power(t0,5)*Power(t1,2)*Power(te,2)*xred + 12.*Power(t0,4)*Power(t1,3)*Power(te,2)*xred - 
        3.0000000000000013*Power(t0,2)*Power(t1,5)*Power(te,2)*xred - 2.*Power(t0,6)*Power(te,3)*xred + 4.*Power(t0,5)*t1*Power(te,3)*xred + 
        4.000000000000001*Power(t0,4)*Power(t1,2)*Power(te,3)*xred - 12.000000000000002*Power(t0,3)*Power(t1,3)*Power(te,3)*xred + 
        6.000000000000001*Power(t0,2)*Power(t1,4)*Power(te,3)*xred + 2.*Power(t0,5)*Power(te,4)*xred - 7.000000000000001*Power(t0,4)*t1*Power(te,4)*xred + 
        8.*Power(t0,3)*Power(t1,2)*Power(te,4)*xred - 3.*Power(t0,2)*Power(t1,3)*Power(te,4)*xred))/
    ((1.*Power(t0,2) - 2.*t0*t1 + 1.*Power(t1,2))*(0. - 0.027777777777777783*Power(t0,4)*Power(t1,3) + 0.08333333333333334*Power(t0,3)*Power(t1,4) - 
        0.08333333333333333*Power(t0,2)*Power(t1,5) + 0.027777777777777776*t0*Power(t1,6) + 0.08333333333333334*Power(t0,4)*Power(t1,2)*te - 
        0.22222222222222224*Power(t0,3)*Power(t1,3)*te + 0.16666666666666666*Power(t0,2)*Power(t1,4)*te - 0.027777777777777776*Power(t1,6)*te - 
        0.08333333333333334*Power(t0,4)*t1*Power(te,2) + 0.16666666666666666*Power(t0,3)*Power(t1,2)*Power(te,2) - 0.1666666666666667*t0*Power(t1,4)*Power(te,2) + 
        0.08333333333333334*Power(t1,5)*Power(te,2) + 0.027777777777777783*Power(t0,4)*Power(te,3) - 0.16666666666666669*Power(t0,2)*Power(t1,2)*Power(te,3) + 
        0.22222222222222227*t0*Power(t1,3)*Power(te,3) - 0.08333333333333334*Power(t1,4)*Power(te,3) - 0.027777777777777783*Power(t0,3)*Power(te,4) + 
        0.08333333333333334*Power(t0,2)*t1*Power(te,4) - 0.08333333333333334*t0*Power(t1,2)*Power(te,4) + 0.027777777777777783*Power(t1,3)*Power(te,4)));

	return res;
}

static void CP_problem(double v0, double x0, double t0) {

	/* UP Problem */
	double te_init, temp_te, check_te, temp_cost = DBL_MAX;
	double checkSpeed = v0;
	int negativeSpeed = 0;
	
	// calculate max te assuming the min speed is 1
	double max_te = (xe - x0) / 1.0;
	te_init = t0 + 1;

	for (int i = 0; i < newtonIntervals; i++) {
		double te_cost = 0.0;

		check_te = newton_raphson(t0, x0, v0, te_init);

		if (check_te > 2.0*max_te || isnan(check_te) || check_te < 0.0)
			continue;

		for (int i = 0; i < P.numsteps; i++) {
			double step = (check_te - t0) / (P.numsteps);
			double t = (i*(check_te - t0) / (P.numsteps - 1) + t0);
			
			te_cost += 0.5*pow(UP_control(t, t0, x0, v0, check_te), 2);
			/* skip optimal te that leads to negative speed trajectories */
			checkSpeed = (3*Power(t,2)*(t0*v0 - check_te*v0 + t0*ve - check_te*ve - 2*x0 + 2*xe))/Power(t0 - check_te,3) - (2*t*(Power(t0,2)*v0 + t0*check_te*v0 - 2*Power(check_te,2)*v0 + 2*Power(t0,2)*ve - t0*check_te*ve - Power(check_te,2)*ve - 3*t0*x0 - 3*check_te*x0 + 3*t0*xe + 3*check_te*xe))/Power(t0 - check_te,3) - (-2*Power(t0,2)*check_te*v0 + t0*Power(check_te,2)*v0 + Power(check_te,3)*v0 - Power(t0,3)*ve - Power(t0,2)*check_te*ve + 2*t0*Power(check_te,2)*ve + 6*t0*check_te*x0 - 6*t0*check_te*xe)/Power(t0 - check_te,3);
			if (checkSpeed < 0.0)
				negativeSpeed = 1;
		}

		if ((check_te > t0) && (te_cost < temp_cost) && (negativeSpeed != 1)) {
			temp_cost = te_cost;
			temp_te = check_te;
		}
		te_init += newton_add_term;
	}
	double te = temp_te;

	// fprintf(stderr, "te: %.4f \n", te);
	double a_up, v_up, x_up;
	double step = (te - t0) / (P.numsteps - 1);
	int violation = 0;
	for (double t = 0; t <= te; t = t + step) {
		a_up = (6*t*(t0*v0 - te*v0 + t0*ve - te*ve - 2*x0 + 2*xe))/Power((t0 - te),3) - (2*(Power(t0,2)*v0 + t0*te*v0 - 2*Power(te,2)*v0 + 2*Power(t0,2)*ve - t0*te*ve - Power(te,2)*ve - 3*t0*x0 - 3*te*x0 + 3*t0*xe + 3*te*xe))/Power((t0 - te),3);
		v_up = (3*Power(t,2)*(t0*v0 - te*v0 + t0*ve - te*ve - 2*x0 + 2*xe))/Power(t0 - te,3) - (2*t*(Power(t0,2)*v0 + t0*te*v0 - 2*Power(te,2)*v0 + 2*Power(t0,2)*ve - t0*te*ve - Power(te,2)*ve - 3*t0*x0 - 3*te*x0 + 3*t0*xe + 3*te*xe))/Power(t0 - te,3) - (-2*Power(t0,2)*te*v0 + t0*Power(te,2)*v0 + Power(te,3)*v0 - Power(t0,3)*ve - Power(t0,2)*te*ve + 2*t0*Power(te,2)*ve + 6*t0*te*x0 - 6*t0*te*xe)/Power(t0 - te,3);
		x_up = (Power(t,3)*(t0*v0 - te*v0 + t0*ve - te*ve - 2*x0 + 2*xe))/Power(t0 - te,3) - (Power(t,2)*(Power(t0,2)*v0 + t0*te*v0 - 2*Power(te,2)*v0 + 2*Power(t0,2)*ve - t0*te*ve - Power(te,2)*ve - 3*t0*x0 - 3*te*x0 + 3*t0*xe + 3*te*xe))/Power(t0 - te,3) - (t*(-2*Power(t0,2)*te*v0 + t0*Power(te,2)*v0 + Power(te,3)*v0 - Power(t0,3)*ve - Power(t0,2)*te*ve + 2*t0*Power(te,2)*ve + 6*t0*te*x0 - 6*t0*te*xe))/Power(t0 - te,3) - (Power(t0,2)*Power(te,2)*v0 - t0*Power(te,3)*v0 + Power(t0,3)*te*ve - Power(t0,2)*Power(te,2)*ve - 3*t0*Power(te,2)*x0 + Power(te,3)*x0 - Power(t0,3)*xe + 3*Power(t0,2)*te*xe)/Power(t0 - te,3);
		// fprintf(stderr, "(%.4f) \t %.4f \t %.4f \t %.4f \n", t, a_up, v_up, x_up);
		if (x_up > 150.0 && t < T_max)
			violation = 1;
	}

	if (violation == 1) {
		// fprintf(stderr, "UP violates the traffic signal \n");

		/* CP Problem */
		// double a_up, v_up, x_up;
		temp_cost = DBL_MAX;
    temp_te = -1.0;
		step = (te - t0) / (P.numsteps - 1);
		double t1 = 30.0, xred = 150.0;
    
		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
		for (int i = 0; i < 10; i++) {
			double te_cost = 0.0;

			check_te = newton_raphson_cp(t0, x0, v0, te_init, t1, xred);
			// fprintf(stderr, "check: %.4f -- init_te: %.4f\n", check_te, te_init);
			if (check_te > 2.0*max_te || isnan(check_te) || check_te < 0.0 || check_te <= T_max) {
				te_init += newton_add_term;
				continue;
			}
			
			for (int i = 0; i < P.numsteps; i++) {
				double step = (check_te - t0) / (P.numsteps);
				double t = (i*(check_te - t0) / (P.numsteps - 1) + t0);
				
				te_cost += 0.5*pow(a2(t, t0, x0, v0, check_te, t1, xred), 2);
				/* skip optimal te that leads to negative speed trajectories */
				checkSpeed = v1(t, t0, x0, v0, check_te, t1, xred);
				if (checkSpeed < 0.0) {
					negativeSpeed = 1;
          // fprintf(stderr, "negative speed \n");
        }
			}

			if ((check_te > t0) && (te_cost < temp_cost) && (negativeSpeed != 1)) {
				temp_cost = te_cost;
				temp_te = check_te;
        // fprintf(stderr, "found te: %.4f \n\n", temp_te);
			}
			te_init += newton_add_term;
		}
    if (temp_te == -1.0) {
      fprintf(stderr, "no optimal te results in positive speed trajectory \n");
      negSpeedInitCon = 1;
      // exit(1);
    }
		te = temp_te;
		// fprintf(stderr, "te: %.4f \n\n", te);
		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    for (int t = 0; t <= 30; t++) {
			initial_u[t] = a1(t, t0, x0, v0, te, t1, xred);
			initial_v[t] = v1(t, t0, x0, v0, te, t1, xred);
			initial_x[t] = x1(t, t0, x0, v0, te, t1, xred);
			// fprintf(stderr, "(%d) \t %.4f \t %.4f \t %.4f \n", t, initial_u[t], initial_v[t], initial_x[t]);
		}
	}

}
/* ------------------------------------------------------------------------------ */

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

double idm_accel(double x, double xl, double v, double vl, double v0, double ts, int leaderID) {
    
    double s = x - xl;
    double sStar = idm.s0 + fmax(v*ts + ((v * (v - vl)) / (2.0 * sqrt(idm.aMax * idm.b))), 0.0);
    double aFree = idm.aMax * (1.0 - pow(v/v0, idm.d));

    double accIDM, accIDM_x, accIDM_v;
    if (leaderID == -1)
        accIDM = aFree;
    else if (leaderID >= 0)
        accIDM = aFree - (idm.aMax * pow(sStar/s, 2));
    else if (leaderID == -2) {
		double a1 = 1.0 - pow(x/xl, 4);
		double a2 = 1.0 - pow(v/v0, 2);
		
		if (a1 == 0.0 && a2 == 0.0)
			accIDM = 0.0;
		else if (a1 == 0.0 && a2 != 0.0)
			accIDM = idm.aMax * a2;
		else if (a1 != 0.0 && a2 == 0.0)
			accIDM = idm.aMax * a1;
		else
        	accIDM = idm.aMax * fmin(a1, a2);

		/* second try */
		// sStar = idm.s0 + fmax(v*ts + ((v * (v - vl)) / (2.0 * sqrt(idm.aMax * idm.b))), 0.0);
		// accIDM = aFree - (idm.aMax * pow(sStar/s, 2));
    }
    else
        fprintf(stderr, "Unknown Leader Type. \n");
    
    
    return accIDM;
}

static void printsol(node_t initial, int it) {
	
	node_t *nextptr;
	int k = 0;
	node_t *state;

	assert((state = (node_t*)calloc(sizeof(node_t), P.numsteps)));
	/* trace sequence of states */
	state[0] = initial;
	ctrl_t c;

	FILE *f;
	if (DISPLAY) {
		fprintf(stderr, "\nk \t Position \t Speed    \t Control  \t Cost \t \t Hybrid  \t J_opt  \t te \t Probability \t Real Cost (1/2u^2) \t Real+escape \t Escape \n");
		fprintf(stderr, "-----------------------------------------------------------------------------------------------------------------\n");
	}
	if (OUTPUT) {
		char buf[0x100];
		snprintf(buf, sizeof(buf), "outputs/(%02d) x0=%.1f - v0=%.1f - delta=%.1f - du=%.3f.txt", it + 1, P.x0, P.v0, DELTA_V, du);
		f = fopen(buf,"w");
		fprintf(f, "k  Position  Speed  Control  Cost  Hybrid  J_opt  te  Probability  Real_Cost_(1/2u^2)  Real+escape   Escape \n");
	}
	
	double hybrid_cost = 0.0;
	double real_cost = 0.0;
	double real_escape_cost = 0.0; // 0.5u^2 + the escape cost in each stage
	for (k = 0; k < P.numsteps; k++) {
		c = ACCESS(C, state[k]);

		/* real cost*/
		real_cost += 0.5 * pow(c.u, 2);

		/* hybrid cost */
		if (k < T_min){
			hybrid_cost += 0.5 * pow(c.u, 2);
			real_escape_cost += 0.5 * pow(c.u, 2);
		}
		else{
			hybrid_cost += (1.0 - p[k]) * 0.5 * pow(c.u, 2) + p[k] * ACCESS(J_opt, state[k]);
			if (STOCHASTIC)
				real_escape_cost = real_cost + ACCESS(J_opt, state[k]);
			else {
				if (k == P.numsteps - 1)
					real_escape_cost = real_cost + ACCESS(J_opt, state[k]);	
				else
					real_escape_cost = real_cost + ACCESS(escape, state[k]);
			}
		}
		
		if (DISPLAY) {
			// fprintf(stderr, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n",
			fprintf(stderr, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n",
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
				real_escape_cost
				// ACCESS(escape, state[k]),
				// D.vmin[k],
				// D.vmax[k],
				// D.xmin[k],
				// D.xmax[k],
				// init_v[k],
				// init_x[k]
				);
		}
		if (OUTPUT) {
			fprintf(f, "%d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n",
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
				ACCESS(escape, state[k]),
				D.vmin[k],
				D.vmax[k],
				D.xmin[k],
				D.xmax[k],
				init_v[k],
				init_x[k],
				init_u[k]
				);
		}

		if (k < P.numsteps - 1) {
			if ((nextptr = getnext(state[k], c)) == NULL) {
				fprintf(stderr, "invalid \n");
				exit(1);
			}
			state[k + 1] = *nextptr;
		}

		opt_x[k] = state[k].x;
		opt_v[k] = state[k].v;
		opt_u[k] = c.u;
	}

	#if (ARRB)
		/* calculate fuel consumption for SDP solution without escape trajectories */
		double fuel = arrb(opt_v, opt_u, P.numsteps, P.T);
		fprintf(stderr, "fuel: %.4f \n", fuel);
		
		/* calculate deterministic GLOSA trajectories */
		double t, te, x0, v0, dt;
		int t0, steps, count;

		int RD = (T_max - T_min);
		// double a_oc[RD], x_oc[RD], v_oc[RD];
		double **a_oc, **x_oc, **v_oc;
		double **a_full, **x_full, **v_full;

		a_oc = (double**)calloc(RD + 1, sizeof(double*)); x_oc = (double**)calloc(RD + 1, sizeof(double*)); v_oc = (double**)calloc(RD + 1, sizeof(double*));
		a_full = (double**)calloc(RD + 1, sizeof(double*)); x_full = (double**)calloc(RD + 1, sizeof(double*)); v_full = (double**)calloc(RD + 1, sizeof(double*));
		for (int i = 0; i < RD + 1; i++) {		
			t0 = T_min + i;
			x0 = state[t0].x;
			v0 = state[t0].v;
			te = ACCESS(te_opt, state[t0]);

			// // steps = (int)(ceil(te) - t0);
			// steps = 9;
			// a_oc[i] = (double*)calloc(steps + 1, sizeof(double)); x_oc[i] = (double*)calloc(steps + 1, sizeof(double)); v_oc[i] = (double*)calloc(steps + 1, sizeof(double));
			// dt = (te - t0) / (steps + 1);
			// v_oc[i][0] = v0;
			// for (int j = 0; j <= steps; j++) {
			// 	t = ((j * dt) + t0);
			// 	a_oc[i][j] = UP_control(t, t0, x0, v0, te);
			// 	// v_oc[i][j] = UP_speed(t, t0, x0, v0, te);
			// 	// x_oc[i][j] = UP_position(t, t0, x0, v0, te);
			// 	v_oc[i][j + 1] = v_oc[i][j] + a_oc[i][j] * (te - t0) / (steps);
			// 	fprintf(stderr,"%d \t (%.4f) \t %.4f \t %.4f \t %.4f \n", t0, t, a_oc[i][j], v_oc[i][j], x_oc[i][j]);
			// }
			// fprintf(stderr,"%d \t (%.4f) \t %.4f \t %.4f \n", t0, t, a_oc[i][steps], v_oc[i][steps]);

			// double fuel_det = arrb(v_oc[i], a_oc[i], steps, dt);
			// fprintf(stderr, "fuel deterministic: (%d) %.4f \n", t0, fuel_det);
			// fprintf(stderr,"\n");
			
			double xk = x0, vk = v0;
			// int opNumsteps = (int)(te - t0) + 1;
			int opNumsteps = P.numsteps;
			int count = 0;
			double a_up[100], v_up[100], x_up[100];

			double step = (te - t0) / opNumsteps;
			for (int i = 0; i < opNumsteps; i++) {
				double t = (i*(te - t0) / (opNumsteps - 1) + t0);
				a_up[count] = UP_control(t, t0, x0, v0, te);
				count++;
			}
			
			for (int i = 0; i < count; i++) {
				// x_up[i + 1] = x_up[i] + v_up[i] * step + 0.5 * pow(step, 2) * a_up[i];
				// v_up[i + 1] = v_up[i] + step * a_up[i];
				// xk += vk * step + 0.5 * pow(step, 2) * a_up[i];
				vk += step * a_up[i];

				v_up[i] = vk;
				// fprintf(stderr, "v: %.4f -- x: %.4f \n", v_up[i+1], x_up[i+1]);
			}
			double fuel_det = arrb(v_up, a_up, count, step);
			fprintf(stderr, "fuel deterministic: (%d) %.4f \n", t0, fuel_det);
		}
	#endif

	opt_J_it = ACCESS(J, state[0]);
	free(state);
	if (OUTPUT)
		fclose(f);
}

static void dp(void) {
	
	int k, ix, iv;
	double x, v;
	double te;
	if (DISPLAY) {
		if (STOCHASTIC)
			fprintf(stderr, "Stochastic DP \n");
		else
			fprintf(stderr, "Deterministic DP \n");
		fprintf(stderr, "NU: %d -- NV: %d -- NX: %d \n", D.NU, D.NV[2], D.NX[2]);
	}
	//fprintf(stderr, "NU: %d -- NV: %d -- NX: %d \n", D.NU, D.NV[2], D.NX[2]);
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
				here.ix = ((k == 0) ? 0 : ix);
				here.iv = ((k == 0) ? 0 : iv);
				
				if (here.ix > D.NX[k] || here.iv > D.NV[k]) {
					fprintf(stderr, "invalid state on the grid: k %d x %.1f/%.1f v %.1f/%.1f\n",
						here.k, here.x, D.X[D.NX[k] - 1][k], here.v, D.V[D.NV[k] - 1][k]);
					exit(1);
				}

				node_t next = { 0 };
				double u, phi;
				double *Jhere, *Jnext;
				ctrl_t *Chere;
				
				if (here.k == 0)
					assert(here.ix == 0 && here.iv == 0);

				double result = 0.0;
				if (T_min <= here.k && here.k <= T_max) {
					double te_init, temp_te, check_te, temp_cost = DBL_MAX;
					double checkSpeed = here.v;
					int negativeSpeed = 0;
					
					if (!FIXED_TIME) {
						/* calculate max te assuming the min speed is 1 */
						double max_te = (xe - here.x) / 1.0;
						te_init = here.k + 1;

						for (int i = 0; i < newtonIntervals; i++) {
							double te_cost = 0.0;

							check_te = newton_raphson(here.k, here.x, here.v, te_init);

							if (check_te > 2.0*max_te || isnan(check_te) || check_te < 0.0)
								continue;

							for (int i = 0; i < P.numsteps; i++) {
								double step = (check_te - here.k) / (P.numsteps);
								double t = (i*(check_te - here.k) / (P.numsteps - 1) + here.k);

								te_cost += 0.5*pow(UP_control(t, here.k, here.x, here.v, check_te), 2);
								/* skip optimal te that leads to negative speed trajectories */
								checkSpeed += UP_control(t, here.k, here.x, here.v, check_te)*P.T;
								if (checkSpeed < 0.0)
									negativeSpeed = 1;
							}

							if ((check_te > here.k) && (te_cost < temp_cost) && (negativeSpeed != 1)) {
								temp_cost = te_cost;
								temp_te = check_te;
							}
							te_init += newton_add_term;
						}
						te = temp_te;
					}
					else
						te = tf;

					te_opt[here.k][here.ix][here.iv] = te;

					// calculate the optimal control cost
					double step;
					for (int i = 0; i < P.numsteps; i++) {
						step = (te - here.k) / (P.numsteps);
						double t = (i*(te - here.k) / (P.numsteps - 1) + here.k);
						result += 0.5* (pow(UP_control(t, here.k, here.x, here.v, te), 2));
					}
					result *= step;

					if ((!STOCHASTIC && k == P.numsteps - 1) || STOCHASTIC) {
						J_opt[here.k][here.ix][here.iv] = wopt * result;
						if (!STOCHASTIC)
							escape[here.k][here.ix][here.iv] = wopt * result;
					}
					else {
						escape[here.k][here.ix][here.iv] = wopt*result;
						J_opt[here.k][here.ix][here.iv] = 0.0;
					}
		
					
					/* uniform distribution */
					p[here.k] = 1.0 / ((T_max - here.k) + 1.0);

					if (!STOCHASTIC) {
						if (here.k == (int)T_max) {
							p[here.k] = 1.0;
						}
						else {
							p[here.k] = 0.0;
						}	
					}
				}
				else
					p[here.k] = 0.0;

				/* final stage costs */
				if (here.k == P.numsteps - 1) {
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
						
					if (discretizexv(&next, next.k) != 0)
						continue;			

					if (crash(next) > 0 && OBSTACLES)
						continue;

					Jnext = &J[next.k][next.ix][next.iv];
					phi = ((1.0 - p[here.k]) * (0.5*pow(u, 2) + *Jnext)) + (p[here.k] * J_opt[here.k][here.ix][here.iv]);
					
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
  int k;

	P.numsteps = (int)((T_max - T_tl)) + 1;
	P.T = TIME_STEP;
  // du = init_du;

  double scenario_x[] = { 0.0,  10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0 };
  double scenario_v[] = { 1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0 };

	// double loop for different initial conditions
	for (int xi = 0; xi < 1; xi++) {
		for (int vi = 0; vi < 1; vi++) {

      du = init_du;
      // P.x0 = scenario_x[xi]; P.v0 = scenario_v[vi];
      P.x0 = 0.0; P.v0 = 5.0;
			fprintf(stderr, "x0: %.4f -- v0: %.4f \n", P.x0, P.v0);
      
      assert((opt_x = (double*)calloc(sizeof(double), P.numsteps)));
      assert((opt_v = (double*)calloc(sizeof(double), P.numsteps)));
      assert((opt_u = (double*)calloc(sizeof(double), P.numsteps)));

      assert((init_x = (double*)calloc(sizeof(double), P.numsteps)));
      assert((init_v = (double*)calloc(sizeof(double), P.numsteps)));
      assert((init_u = (double*)calloc(sizeof(double), P.numsteps)));

      assert((initial_x = (double*)calloc(sizeof(double), P.numsteps)));
      assert((initial_v = (double*)calloc(sizeof(double), P.numsteps)));
      assert((initial_u = (double*)calloc(sizeof(double), P.numsteps)));
      
			// calculate initial trajectory (deterministic solution of pessimistic case)
      negSpeedInitCon = 0;
			CP_problem(P.v0, P.x0, 0.0);
      if (negSpeedInitCon == 1) continue;
      
			// for (int k = 0; k < P.numsteps; k++)
			// 	fprintf(stderr, "(%d) \t init_x: %.4f \t init_v: %.4f \t init_a: %.4f \n", k, initial_x[k], initial_v[k], initial_u[k]);

      for (k = 0; k < P.numsteps; k++) {
        opt_x[k] = initial_x[k];
        opt_v[k] = initial_v[k];
        opt_u[k] = initial_u[k];
      }

      node_t initial = { 0 };
      initial.k = 0;
      initial.x = P.x0;
      initial.v = P.v0;
      initial.ix = initial.iv = 0;
      
      clock_t start_all, end_all;
      start_all = clock();
      double opt_J_prev_it = DBL_MAX, du_it_prev = du;
      int it = 0;
      while(true){
        
        for (k = 0; k < P.numsteps; k++) {
          init_x[k] = opt_x[k];
          init_v[k] = opt_v[k];
          init_u[k] = opt_u[k];
        }
        
        init_D(P.T, it);
        init_JC();

        /* run dp */
        clock_t start_it, end_it;
        start_it = clock();

        dp();

        end_it = clock();
        double cpu_time_used_it = ((double)(end_it - start_it)) / CLOCKS_PER_SEC;
        
        /*  prints */
        printsol(initial, it);
        if (DISPLAY) {
          fprintf(stderr, "\n--------------------------------\n");
          fprintf(stderr, "DU: %.2f -- DV: %.2f DX: %.2f \n", D.du, D.dv, D.dx);
          fprintf(stderr, "Iteration %d CPU time: %.4f \n", it + 1, cpu_time_used_it);
          fprintf(stderr, "--------------------------------\n");
        }
        else
          fprintf(stderr, "Iteration %d -- DU=DV: %.5f -- DX: %.5f -- cost: %.5f -- CPU time: %.4f \n", it + 1, D.du, D.dx, opt_J_it, cpu_time_used_it);
      
        /* terminal criterion: || u(k) - u(k-1) || < epsilon */
        double sum = 0.0;
        for (int i = 0; i < P.numsteps; i++)
          sum += pow(opt_u[i] - init_u[i], 2);
        double term_criterion = sqrt(sum);

        /* check if reduction of du is needed */
        du_it_prev = du;
        if (opt_J_it == opt_J_prev_it)
          du /= 2.0;

        free_JC();
        if (opt_J_it > opt_J_prev_it) {
          fprintf(stderr, "Termination error: current cost > previous cost !!! \n");
          break;
        }
        if (TERM_CRITERION_DU && (du < 0.125)) {
          fprintf(stderr, "Termination criterion: DU < 0.125 \n");
          break;
        }

        opt_J_prev_it = opt_J_it;
        it++;
        
      }
      end_all = clock();
      double cpu_time_used_all = ((double)(end_all - start_all)) / CLOCKS_PER_SEC;
      fprintf(stderr, "Total CPU time (in %d it.): %.4f \n", it + 1, cpu_time_used_all);

      #if (IDM)
        idm.aMax = 3.0; idm.b = 3.0; idm.d = 4; 
        idm.s0 = idm.aMax;

        FILE *f_idm;
        char buf[0x100];
        snprintf(buf, sizeof(buf), "outputs/IDM - x0=%.1f - v0=%.1f.txt", P.x0, P.v0);
        f_idm = fopen(buf,"w");
        // fprintf(f, "k  Position  Speed  Control  Cost  Real+escape   Escape \n");
        assert((idm_x = (double*)calloc(sizeof(double), 3*P.numsteps)));
        assert((idm_v = (double*)calloc(sizeof(double), 3*P.numsteps)));
        assert((idm_u = (double*)calloc(sizeof(double), 3*P.numsteps)));

        idm_x[0] = P.x0; idm_v[0] = P.v0;
        for (int k = 0; k < P.numsteps; k++) {
          // calculate IDM trajectory 
          idm_u[k] = idm_accel(idm_x[k], XMAX + idm.s0, idm_v[k], 0.0, VMAX, 1.0, 0);
          idm_v[k + 1] = idm_v[k] + P.T * idm_u[k];
          idm_x[k + 1] = idm_x[k] + P.T * idm_v[k] + 0.5 * pow(P.T, 2) * idm_u[k];
          double aCon;
          if (idm_v[k + 1] <= 0.0) {
            idm_u[k] = 0.0; aCon = idm.b;
            idm_x[k + 1] = idm_x[k] + pow(idm_v[k], 2) / (2.0 * aCon);
            idm_v[k + 1] = 0.0;
          }
          // fprintf(stderr, "(%d) idm_x: %.4f -- idm_v: %.4f -- idm_a: %.4f \n", k, idm_x[k], idm_v[k], idm_u[k]);
        }
        
        // idm_x[0] = idm_x[P.numsteps - 1]; idm_v[0] = idm_v[P.numsteps - 1];
        int count = 0;
        for (int k = P.numsteps - 1; k < 3*P.numsteps; k++) {
          idm_u[k] = idm_accel(idm_x[k], 220.0 + idm.s0, idm_v[k], 11.0, ve, 1.0, -2);
          idm_v[k + 1] = idm_v[k] + P.T * idm_u[k];
          idm_x[k + 1] = idm_x[k] + P.T * idm_v[k] + 0.5 * pow(P.T, 2) * idm_u[k];
          double aCon;
          if (idm_v[k + 1] <= 0.0) {
            idm_u[k] = 0.0; aCon = idm.b;
            idm_x[k + 1] = idm_x[k] + pow(idm_v[k], 2) / (2.0 * aCon);
            idm_v[k + 1] = 0.0;
          }
          // fprintf(stderr, "(%d) idm_x: %.4f -- idm_v: %.4f -- idm_a: %.4f \n", k, idm_x[k], idm_v[k], idm_u[k]);
          if (idm_x[k + 1] >= xe) {
            // fprintf(stderr, "(%d) idm_x: %.4f -- idm_v: %.4f \n", k + 1, idm_x[k + 1], idm_v[k + 1]);
            break;
          }
          count++;
        }

        for (int i = 0; i < P.numsteps + count + 1; i++)
          fprintf(stderr, "(%d) idm_x: %.4f -- idm_v: %.4f -- idm_a: %.4f  \n", i, idm_x[i], idm_v[i], idm_u[i]);

        fclose(f_idm);
        double idm_fuel = arrb(idm_v, idm_u, P.numsteps + count + 1, P.T);
        fprintf(stderr, "idm fuel: %.4f \n", idm_fuel);
      #endif

      free(opt_x); free(opt_v); free(opt_u);
      free(init_x); free(init_v); free(init_u);
      free(initial_x); free(initial_v); free(initial_u);


	    for (int i = 0; i < P.numsteps; i++) {
		    free(D.V[i]);
        free(D.X[i]);
      }
      free(D.V);
      free(D.X);

      free(D.vmax); free(D.vmin); free(D.xmax); free(D.xmin);
      free(D.NV); free(D.NX);
      free(D.NXlive_min); free(D.NXlive_max); free(D.NVlive_min); free(D.NVlive_max);
      
      fprintf(stderr, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n ");
    }
  }
	return 0;
}
