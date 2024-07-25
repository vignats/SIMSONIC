#include <stdio.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#ifdef MEX
#include "mex.h"
#endif


#define __min(a,b) ((a)>(b) ? (b) : (a))


#define _RECIPES_INFINITE 1.0E+100
#define SHFT(a,b,c,d)   (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)       ((b) > 0.0 ? fabs(a) : -fabs(a))
#define ITMAX           100
#define CGOLD           0.3819660
#define ZEPS            1.0e-10
#define TEST            fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)
#define GOLD            1.618034
#define GLIMIT          100.0
#define TINY            1.0e-20
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#define TOL             2.0e-4
#define ITMAX2          200
#define EPS             1.0e-10


inline
double sign (double x)

{
        return x >= 0. ? 1 : -1;
}

inline double
sqr (double x)

{
	return x * x;
}

class vecteur {

public:
        vecteur (double lx = 0.,
                 double lz = 0.);

        vecteur operator + (const vecteur & second);
        vecteur operator - (const vecteur & second);
        friend vecteur operator * (double l,
                                   const vecteur & second);
	vecteur operator / (double l);

        void disp (const char * header);

        double x;
        double z;
};

double
dot (const vecteur & first,
     const vecteur & second)

{
        return first.x * second.x + first.z * second.z;
}

double
module (const vecteur & v)

{
        return sqrt(sqr(v.x) + sqr(v.z));
}

void
vecteur::disp (const char * header)

{
        printf("%s: %.10f %.10f\n", header, x, z);
}

vecteur::vecteur (double lx,
                  double lz)

{
        x = lx;
        z = lz;
}

vecteur
vecteur::operator / (double l)

{
	return vecteur(x / l, z / l);
}

vecteur
vecteur::operator + (const vecteur & second)

{
        return vecteur(x + second.x, z + second.z);
}

vecteur
vecteur::operator - (const vecteur & second)

{
        return vecteur(x - second.x, z - second.z);
}

vecteur
operator * (double l,
            const vecteur & second)

{
        return vecteur(l * second.x, l * second.z);
}





double 
brent (double ax, 
       double bx, 
       double cx,
       double (*f)(double, void *), 
       void * localBuffer, 
       double tol, 
       double * xmin)

{
        int iter;
        double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
        double tmp1, tmp2, e = 0.0;

        a = ax < cx ? ax : cx;
        b = ax > cx ? ax : cx;
        x = bx;
	w = bx;
	v = bx;
	fx = (*f)(x, localBuffer);
        fw = fx;
	fv = fx;
        if (fw == _RECIPES_INFINITE) {
		return fw;
	}

	d = 0.0;
        for (iter = 1; iter <= ITMAX; iter++) {
                xm = 0.5 * (a + b);
		tol1 = tol * fabs(x) + ZEPS;
                tol2 = 2.0 * tol1;
                
                if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                        *xmin = x;
                        return fx;
                }
                if (fabs(e) > tol1) {
                        tmp1 = x - w; 
			tmp2 = x - v;
                        r = tmp1 * (fx - fv);
                        q = tmp2 * (fx - fw);
                        p = tmp2 * q - tmp1 * r;
                        q = 2.0 * (q - r);
                        if (q > 0.0) p = -p;
                        q = fabs(q);
                        etemp = e;
                        e = d;
                        if (TEST) {
				e = x >= xm ? a - x : b - x;
				d = CGOLD * e;
                        }
			else {
                                d = p / q;
                                u = x + d;
                                if (u - a < tol2 || b - u < tol2) d = SIGN(tol1, xm - x);
                        }
                }
                else {
			e = x >= xm ? a - x : b - x;
			d = CGOLD * e;
                }
		u = fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d);
                fu = (*f)(u, localBuffer);
                if (fu == _RECIPES_INFINITE) {
			return fu;
		}
                if (fu <= fx) {
                        if (u >= x) a = x; 
			else        b = x;
                        SHFT(v, w, x, u)
                        SHFT(fv, fw, fx, fu)
                }
                else {
                        if (u < x) a = u; 
			else       b = u;
                        if (fu <= fw || w == x) {
                                v = w;
                                w = u;
                                fv = fw;
                                fw = fu;
                        }
                        else if (fu <= fv || v == x || v == w) {
                                v = u;
                                fv = fu;
                        }
                }
        } 
        
        *xmin = x;
        return fx;
}






int
calcul_refraction (vecteur Ui,                  // vecteur incident norm??
                   vecteur N,                   // direction normale
                   double c12,                  // rapport vitesses avant et apr??s r??fraction
                   vecteur * Ur)                // vecteur r??fract?? norm??

{
        double cos1, sin1, sin2, cos2, sin22;
        int ret = 0;
        vecteur T(N.z, -N.x);

        cos1 = dot(Ui, N);                      // composante incidente normale
        sin1 = dot(Ui, T);                      // composante incidente tangentielle
        sin2 = sin1 / c12;
        sin22 = sqr(sin2);
        if (sin22 <= 1.0) {
                cos2 = sign(cos1) * sqrt(1.0 - sin22);
                *Ur = sin2 * T + cos2 * N;
                ret = 1;
        }
        return ret;
}



int
calcul_reflection (vecteur Ui,                  // vecteur incident norm??
                   vecteur N,                   // direction normale
                   vecteur * Ur)                // vecteur r??fract?? norm??

{
    double cos1, sin1;
    int ret = 0;
    vecteur T(N.z, -N.x);
    
    cos1 = dot(Ui, N);                      // composante incidente normale
    sin1 = dot(Ui, T);                      // composante incidente tangentielle
    if (fabs(sin1) < 1.0) {
        *Ur = sin1 * T - cos1 * N;
        ret = 1;
    }
    return ret;
}






typedef struct {

	double 	c_l;		// vitesse dans la lentille
	double 	h_l;		// ??paisseur lentille

	double 	c_t;		// vitesse dans le tissu mou


	double 	c_l_c_t;	// c_l / c_t
    
	double 	a, b, c;	// param??tres de la premi??re parabole
	double 	theta_min; 	// angles extr??mes ?? explorer
	double 	theta_max;

	double 	xmin, xmax;
	
} GPARAMS;

class LPARAMS {

public:
    LPARAMS () {
    }
	double 	xs, zs;		//??source
	double 	xr, zr;		//??recepteur
	int	check_snell;

};



class  LB {

public:
    LB () {
        verbose = 0;
    }
	GPARAMS * 	gp;
	LPARAMS *	lp;

    int         calc_distance;
    int         verbose;
    
} ;







double // find length of the ray between two interfaces
ray_length (double a,
            double b,
            double c,
            vecteur Pt,
            vecteur Ur)

{
	double alpha, beta, gamma, delta, lambda1, lambda2;

	alpha = a * sqr(Ur.x);
	beta = 2 * a * Pt.x * Ur.x + b * Ur.x - Ur.z;
	gamma = a * sqr(Pt.x) + b * Pt.x + c - Pt.z;
        if (alpha == 0) {                                       // parabola is a straight line OR vertical propagation
		if (beta == 0.) {				// ray is parallel to the interface
			return 0;
		}
        	gamma = a * sqr(Pt.x) + b * Pt.x + c - Pt.z;
		lambda1 = - gamma / beta;			// lambda1 is the length of the ray from Pt to interface
	}
	else {							// real parabola
		delta = sqr(beta) - 4. * alpha * gamma;
		if (delta < 0.) {				// no solution 
			return 0;
		}
		delta = sqrt(delta);
		lambda1 = 2. * gamma / (- beta - delta);
		lambda2 = 2. * gamma / (- beta + delta);
		if (lambda1 <= 0) {				// non acceptable solution
			lambda1 = lambda2;
		}
		else if (lambda2 < lambda1 && lambda2 > 0) {	// min value expected
			lambda1 = lambda2;
		}
	}

	return lambda1 >= 0 ? lambda1 : 0;
}

/*
double // find length of the ray between two interfaces
ray_length(double a,
           double b,
           double c,
           vecteur Pt, // exit point described by vector
           vecteur Ur) // exit ray direction decribed by unit vector

{
    
    double lambda1;
    if (fabs(a) < 1e-10) {			// cas particulier d'une parabole qui est une droite
        double gamma = b * Pt.x + c - Pt.z;
        double beta = b * Ur.x - Ur.z;
        if (beta == 0.) {
            return 0.;
        }
        lambda1 = -gamma / beta; // lambda1 is length of ray in soft tissue, the zero of linear function of lambda1 is searched
    }
    else {					// vraie parabole
        double gamma = a * sqr(Pt.x) + b * Pt.x + c - Pt.z;
        if (fabs(Ur.x) < 1e-10) {
            lambda1 = gamma / Ur.z;
        }
        else {
            double alpha = a * sqr(Ur.x);
            double beta = 2 * a * Pt.x * Ur.x + b * Ur.x - Ur.z;
            double delta = sqr(beta) - 4 * alpha * gamma;
            if (delta < 0.) {
                return 0.;
            }
            delta = sqrt(delta);
            lambda1 = (-beta - delta) / (2 * alpha);
            double lambda2 = (-beta + delta) / (2 * alpha);
            if (lambda1 < 0 && lambda2 >= 0) {
                lambda1 = lambda2; // lambda1 (or lambda2) is length of ray, the zeros of quadratic function of lambda1 are searched
            }
            else if (lambda2 < lambda1 && lambda2 > 0) {
                lambda1 = lambda2; // lambda1 (or lambda2) is length of ray, the zeros of quadratic function of lambda1 are searched
            }
//            else if (lambda1 < 0 && lambda2 < 0) {
//                return 0.;
//            }
        }
    }
    return lambda1;
}
*/








double 
tps (double theta,
     void * buffer)

{
	LB * lb = (LB *)buffer;
	double t = 0;
	GPARAMS * gp = lb->gp;
	LPARAMS * lp = lb->lp;
    int cd = lb->calc_distance;
   	int ret;
    double lambda1;
    
// premier trajet, de la source jusqu'?? la lentille
        vecteur Ur, Ui(sin(theta), cos(theta));
        vecteur N(0., 1.);
        double tan_theta = tan(theta);
        vecteur Pe(lp->xs, lp->zs);
        vecteur Pl1 = Pe + vecteur(gp->h_l * tan_theta, gp->h_l);
//    if (Pl1.x < 0) {
//        return 1. + fabs(theta);
//    }
        t = gp->h_l * sqrt(sqr(tan_theta) + 1.);
    
	if (!cd) {
		t /= gp->c_l;
	}


    // refraction at interface lens - soft tissue
    ret = calcul_refraction(Ui, N, gp->c_l_c_t, &Ur); // unit vector Ur (transmit ray) is updated
    if (ret == 0) {
        return 1. + fabs(theta);
    }
    lambda1 = ray_length(gp->a, gp->b, gp->c, Pl1, Ur);
    if (lambda1 == 0) {
        return 1. + fabs(theta);
    }
    vecteur Pb1 = Pl1 + lambda1 * Ur;
    t += cd ? lambda1 : lambda1 / gp->c_t;


    
    
    // reflection on periosteum
    Ui = Ur;
    N = vecteur(2. * gp->a * Pb1.x + gp->b, -1.);
    N = N / module(N);
    ret = calcul_reflection(Ui, N, &Ur); // unit vector Ur (transmit ray) is updated
    if (ret == 0) {
        return 1. + fabs(theta);
    }
    //printf(">> %g \n", Ur.z);
    if (Ur.z > 0) {
        return 1. + fabs(theta);
    }
    lambda1 = ray_length(0., 0., gp->h_l, Pb1, Ur);
    if (lambda1 == 0) {
        return 1. + fabs(theta);
    }
    vecteur Pl2 = Pb1 + lambda1 * Ur;
    t += cd ? lambda1 : lambda1 / gp->c_t;
    

    
    
    
    // de la deuxi??me parabole vers le point de l'image
    vecteur R(lp->xr, lp->zr);
    t += cd ? module(R - Pl2) : module(R - Pl2) / gp->c_l;

//    if (Pl2.x < 0) {
//        return 1. + fabs(theta);
//    }
    
	if (lp->check_snell) {
        
	FILE * toto = fopen("toto", "w");
	fprintf(toto, "%g %g\n", lp->xs, lp->zs);
	fprintf(toto, "%g %g\n", Pl1.x, Pl1.z);
	fprintf(toto, "%g %g\n", Pb1.x, Pb1.z);
	fprintf(toto, "%g %g\n", Pl2.x, Pl2.z);
	fprintf(toto, "%g %g\n", lp->xr, lp->zr);

	fclose(toto);
        
// 		N = vecteur(2. * gp->a * Pb1.x + gp->b, -1.);
// 		N = N / module(N);
// 		Ui = (Pb1 - Pl) / module(Pb1 - Pl);
// 		Ur = (R - Pb1) / module(R - Pb1);
// 		double cos_i = dot(Ui, N);
// 		double cos_r = dot(Ur, N);
// 		double sin_i = sqrt(1. - sqr(cos_i)) / gp->c_t;
// 		double sin_r = sqrt(1. - sqr(cos_r)) / gp->c_b;
// 		printf(">> %g %g / %g\n", sin_i, sin_r, sin_i - sin_r);
	}

	return t;
}

double
temps_de_vol (LB & lb,
	      double xs,
	      double zs,
	      double xr,
	      double zr,
	      double & thmin,
          double & dist)

{
//	double zp = lb.gp->a * sqr(xr) + lb.gp->b * xr + lb.gp->c;
//	if (zr < zp) {
//		thmin = 0;
//		dist = 0;
//		return 0;
//	}

	lb.lp->xs = xs;
	lb.lp->zs = zs;
	lb.lp->xr = xr;
	lb.lp->zr = zr;

    	lb.calc_distance = 0;
    	double tmin;
    

/*
	int N = 1001;
	FILE * out = fopen("out", "w");
	double th[N], t[N], dth = (lb.gp->theta_max - lb.gp->theta_min) / (N - 1);
	for (int i = 0; i < N; i++) {
		th[i] = lb.gp->theta_min + i * dth;
		t[i] = tps(th[i], &lb);
		//if (t[i] >= 1.) {
		//	t[i] = 0.;
		//}
		fprintf(out, "%g %g\n", th[i], t[i]);
	}
	fclose(out);
*/


/*
    lb.verbose = 1;
    FILE * fp = fopen("toto", "w");
    int N = 1000;
    double th, t, min = lb.gp->theta_min, max = lb.gp->theta_max, dth = (max - min) / N;
    for (int i = 0; i < N; i++) {
        th = min + i * dth;
        t = tps(th, &lb);
        if (t >= 1.) {
            t = 0.;
        }
        fprintf(fp, "%g %g\n", th, t);
    }
    fclose(fp);
*/
    
    int Nb_pts = 9;
    double theta_initial, theta_search, travel_time, travel_time_min, dtheta_search = (lb.gp->theta_max - lb.gp->theta_min) / (Nb_pts - 1);
    travel_time_min = 1.;
    for (int i = 0; i < Nb_pts; i++) {
        theta_search = lb.gp->theta_min + i * dtheta_search;
        travel_time = tps(theta_search, &lb);
        if (travel_time < travel_time_min) {
            theta_initial = theta_search;
            travel_time_min = travel_time;
        }
    }
    
    lb.verbose = 0;
	//tmin = brent(lb.gp->theta_min, 0., lb.gp->theta_max, tps, &lb, 1e-15, &thmin);
    tmin = brent(theta_initial-dtheta_search, theta_initial, theta_initial+dtheta_search, tps, &lb, 1e-15, &thmin);
    	lb.calc_distance = 1;
	lb.lp->check_snell = 0;
    	dist = tps(thmin, &lb);
    	lb.calc_distance = 0;
	lb.lp->check_snell = 0;

/*
	printf("theta = %e\n", thmin);
	printf("temps = %e\n", tmin);
	printf("dist  = %e\n", dist);
*/
 
    	return tmin;
}

typedef struct {

	GPARAMS *	p_gp;

	int		NS;
	double *	XS;
	int 		NX;
	double *	X;
	int 		NZ;
	double *	Z;

	int		sfrom;
	int		sto;

    	double *	T;
    	double *	theta;
    	double *	D;

} THDATA;

void *
th (void * buffer)

{
	THDATA * p_thdata = (THDATA *)buffer;
	LB lb;
	LPARAMS lp;
	int i, j, s, ofs;
	double thmin, dist;

	lb.lp = &lp;
	lb.gp = p_thdata->p_gp;
    

	lp.check_snell = 0;

#ifdef MEX
	for (j = 0; j < p_thdata->NZ; j++) {
		for (i = 0; i < p_thdata->NX; i++) {
			for (s = p_thdata->sfrom; s <= p_thdata->sto; s++) {
				ofs = s * p_thdata->NX * p_thdata->NZ + i * p_thdata->NZ + j;
				p_thdata->T[ofs] = temps_de_vol(lb, p_thdata->XS[s], 0., p_thdata->X[i], p_thdata->Z[j], thmin, dist);
                		if (p_thdata->theta != NULL) {
                    			p_thdata->theta[ofs] = thmin;
                		}
                		if (p_thdata->D != NULL) {
                    			p_thdata->D[ofs] = dist;
                		}
			}
		}
	}
#else
	ofs = p_thdata->sfrom * p_thdata->NX * p_thdata->NZ;
	for (s = p_thdata->sfrom; s <= p_thdata->sto; s++) {
		for (i = 0; i < p_thdata->NX; i++) {
			for (j = 0; j < p_thdata->NZ; j++) {
				p_thdata->T[ofs++] = temps_de_vol(lb, p_thdata->XS[s], 0., p_thdata->X[i], p_thdata->Z[j], thmin, dist);
			}
		}
	}
#endif

	delete p_thdata;
	return NULL;
}

#ifdef MEX

static double
extract_double (const mxArray * cell_ptr)

{
        return *mxGetPr(cell_ptr);
}

void
mexFunction (int nlhs, mxArray * plhs[],
             int nrhs, const mxArray * prhs[])

{
	GPARAMS gp;
	const mwSize * dims;
	double * p;
	int i, j, s;
	THDATA * p_thdata;

        if (nlhs < 1 || nlhs > 3) {
                printf(">> you must specify between 1 and 3 variable(s) to receive output(s)\n");
                return;
        }
	if (nrhs != 6) {
		printf(">> usage: rayons(...)\n");
        printf(">> parameter #1 must contain the 2 velocities cl, ct\n");
		printf(">> parameter #2 must contain the 4 geometrical parameters a, b, c and L\n");
		printf(">> parameter #3 must contain the horizontal positions of the emitters\n");
		printf(">> parameter #4 must contain the horizontal positions of the calculation points\n");
		printf(">> parameter #5 must contain the vertical positions of the calculation points\n");
		printf(">> parameter #6 must contain the number of computation cores\n");
		return;
	}

	dims = mxGetDimensions(prhs[0]);
	if (dims[0] != 1 && dims[1] != 2) {
        printf(">> parameter #1 must contain the 2 velocities cl and ct\n");
		return;
	}
	p = (double *)mxGetData(prhs[0]);
	gp.c_l = p[0];
	gp.c_t = p[1];

//     printf("<< %g %g\n", gp.c_l, gp.c_t);
	gp.c_l_c_t = gp.c_l /gp.c_t;

	dims = mxGetDimensions(prhs[1]);
	if (dims[0] != 1 && dims[1] != 4) {
		printf(">> parameter #2 must contain the 4 geometrical parameters a, b, c and L\n");
		return;
	}
	p = (double *)mxGetData(prhs[1]);
	gp.a = p[0];
	gp.b = p[1];
	gp.c = p[2];
    gp.h_l = p[3];
//     printf("<< %g %g %g %g\n", gp.a, gp.b, gp.c, gp.h_l);
	
	if ( (gp.c_l / gp.c_t) >= 1 ) {
		gp.theta_max = 1.57;
		}
		else {
		gp.theta_max = asin(gp.c_l / gp.c_t);
		}
	gp.theta_min = -gp.theta_max;

	dims = mxGetDimensions(prhs[2]);
	if (dims[0] != 1 && dims[1] <= 0) {
		printf(">> parameter #3 must contain the horizontal positions of the emitters\n");
		return;
	}
	int NS = dims[1];
	double XS[NS];
	p = (double *)mxGetData(prhs[2]);
	memcpy(XS, p, NS * sizeof(double));
// 	printf("<< ");
// 	for (i = 0; i < NS; i++) {
// 		printf("%g ", XS[i]);
// 	}
// 	printf("\n");

	dims = mxGetDimensions(prhs[3]);
	if (dims[0] != 1 && dims[1] <= 0) {
		printf(">> parameter #4 must contain the horizontal positions of the calculation points\n");
		return;
	}
	int NX = dims[1];
	double X[NX];
	p = (double *)mxGetData(prhs[3]);
	memcpy(X, p, NX * sizeof(double));
// 	printf("<< ");
// 	for (i = 0; i < NX; i++) {
// 		printf("%g ", X[i]);
// 	}
// 	printf("\n");

	dims = mxGetDimensions(prhs[4]);
	if (dims[0] != 1 && dims[1] <= 0) {
		printf(">> parameter #5 must contain the vertical positions of the calculation points\n");
		return;
	}
	int NZ = dims[1];
	double Z[NZ];
	p = (double *)mxGetData(prhs[4]);
	memcpy(Z, p, NZ * sizeof(double));
// 	printf("<< ");
// 	for (i = 0; i < NZ; i++) {
// 		printf("%g ", Z[i]);
// 	}
// 	printf("\n");

	dims = mxGetDimensions(prhs[5]);
	if (dims[0] != 1 && dims[1] != 1) {
		printf(">> parameter #6 must contain the number of computation cores\n");
		return;
	}
	p = (double *)mxGetData(prhs[5]);
	int NC = int(p[0] + 0.5);
// 	printf("<< %d\n", NC);

	const mwSize dims_out[] = { NZ, NX, NS };
	plhs[0] = mxCreateNumericArray(3, dims_out, mxDOUBLE_CLASS, mxREAL);
	double * T = (double *)mxGetData(plhs[0]);
    	double * p_theta = NULL;
    	double * p_dist = NULL;
    	if (nlhs >= 2) {
        	plhs[1] = mxCreateNumericArray(3, dims_out, mxDOUBLE_CLASS, mxREAL);
        	p_theta = (double *)mxGetData(plhs[1]);
		memset(p_theta, 0, NX * NS * NZ * sizeof(double));
		for (int toto = 0; toto < NX * NS * NZ; toto++) {
			p_theta[toto] = -100;
		}
    	}
    	if (nlhs == 3) {
        	plhs[2] = mxCreateNumericArray(3, dims_out, mxDOUBLE_CLASS, mxREAL);
        	p_dist = (double *)mxGetData(plhs[2]);
		memset(p_dist, 0, NX * NS * NZ * sizeof(double));
    	}
    
        int from, to, ex, core, coremax, ns;
        ns = NS / NC;
        pthread_t * p_tid = new pthread_t [NC];
    
        from = 0;
        ex = NS % NC;
        for (core = 0; core < NC; core++) {
            	to = core == NC-1? NS-1: __min(from + ns-1, NS-1);
            	if (ex > 0){
                	to++;
                	ex--;
            	}
                p_thdata = new THDATA;
                p_thdata->T = T;
                p_thdata->theta = p_theta;
                p_thdata->D = p_dist;
                p_thdata->p_gp = &gp;
                p_thdata->NS = NS;
                p_thdata->XS = XS;
                p_thdata->NX = NX;
                p_thdata->X = X;
                p_thdata->NZ = NZ;
                p_thdata->Z = Z;
                p_thdata->sfrom = from;
                p_thdata->sto = to;
            
//                 printf("core %d: %d --> %d\n", core, p_thdata->sfrom, p_thdata->sto);
                pthread_create(p_tid + core, NULL, th, p_thdata);
                coremax = core;
            	from = to + 1;
            	if (from >= NS) {
                	break;
            	}
        }
        for (core = 0; core <= coremax; core++) {
                pthread_join(p_tid[core], NULL);
        }
}

// #else
// 
// int
// main (void)
// 
// {
// 	GPARAMS gp;
// 
// 	gp.c_l = 1020;
// 	gp.c_t = 1540;
// 	//gp.c_b = 3200;
// 	gp.c_l_c_t = gp.c_l /gp.c_t;
// 	//gp.c_t_c_b = gp.c_t /gp.c_b;
// 
// 	gp.a = 12.6771;
// 	gp.b = -0.2274;
// 	gp.c = 0.0060;
// 	gp.h_l = 0.0016;
// 
// 	gp.xmin = 0;
// 	gp.xmax = 0.030;
// 	gp.theta_max = asin(gp.c_l / gp.c_t);
// 	gp.theta_min = -gp.theta_max;
// 
// 	double XS[1] = { 0.0177 };
// 	double X[1] = { 0.0160 };
// 	double Z[1] = { 0.0130 };
// 
// 	LPARAMS lp;
// 	GNUplot gnuplot, gnuplot2;
// 	LB lb;
// 
// 	lb.lp = &lp;
// 	lb.gp = &gp;
// 	
// 
// 	lp.xs = XS[0];
// 	lp.zs = 0.;
// 	lp.xr = X[0];
// 	lp.zr = Z[0];
// 	lp.check_snell = 1;
// 
// 
// 	printf("%g\n", tps(-9.527741e-02, &lb));
//           
// 	int dum;
// 	scanf("%d", &dum);
// 
// 	return 0;
// }

#endif
