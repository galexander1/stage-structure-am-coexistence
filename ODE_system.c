#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

//compile with:
// gcc -lgsl -lm -o O ODE_system.c -L /usr/local/lib -I /usr/local/include/

struct params_rec {
	double r;	/* maximum resource biomass */
	double K;
	double sigma1;	/* assimilation efficiency of species 1 */
	double sigma2;	/* assimilation efficiency of species 2 */
	double H1;	/* half-saturation species 1 functional response */
	double H2;	/* half-saturadiot species 2 functional response */
	double Imax1;	/* maximum ingestion rate species 1 */
	double Imax2;	/* maximum ingestion rate species 2 */
	double d1;
	double d2;
	double T;
	double sm;	/* mass of consumer */
};

/* this is a list of all the parameter combinations to test */
struct params_rec test_params[] = {
	{
	.r = 1.0,
	.K = 6.0,
	.sigma1 = 0.5,
	.sigma2 = 0.5,
	.H1 = 3.0,
	.H2 = 2.5,
	.Imax1 = 1.0,
	.Imax2 = 1.1,
	.d1 = 0.015,
	.d2 = 0.04,
	.T = 0.1,
	.sm = 1e-4,
	}
};

int func_ode(double t, const double y[], double f[], void *p_);
void simulate_ODE(struct params_rec *p, double *y, unsigned int Tf, int output);
double fitness(struct params_rec *p, double R, int species);


int main()
{
	FILE *file = fopen("CASE0_PDE.dat","w");

	for (int i = 0; i < sizeof test_params/sizeof test_params[0]; i++) {
		struct params_rec *p = &test_params[i];
		double y[3] = {1,1,1};
		simulate_ODE(p, y, 50000, /*output=*/0);
		simulate_ODE(p, y, 1000, /*output=*/1);


		double Rmin1 = 0.035;//sp2 alone.
		double Rmax1 = 3.98;
		double Rmin2 = 0.09;//sp1 alone.
		double Rmax2 = 3.27;
		double Rmin = 0.0;
		double Rmax = 4.0;

		int N = 10000;
		for (int i = 0; i < N; i++) {
			double R = Rmin + ((Rmax-Rmin)/(double)(N-1))*((double)i);
			if (R >= Rmin1 && R <= Rmax1) {
				if (R >= Rmin2 && R <= Rmax2) {
					fprintf(file,"%f %f %f\n",R, fitness(p, R, 1), fitness(p, R, 2));
				} else {
					fprintf(file,"%f %f %f\n",R, fitness(p, R, 1), nan(""));
				}
			} else {
				if (R >= Rmin2 && R <= Rmax2) {
					fprintf(file,"%f %f %f\n",R, nan(""), fitness(p, R, 2));
				} else {
					fprintf(file,"%f %f %f\n",R, nan(""), nan(""));
				}
			}
		}
		fclose(file);
	}

	return 0;
}


double fitness(struct params_rec *p, double R, int species)
{
	double fit = 0.0;
	if (species == 1) {
		fit = 0.05*(p->sigma1*R*p->Imax1/(p->H1+R)-p->T-p->d1);
	} else {
		fit = 0.05*(p->sigma2*R*p->Imax2/(p->H2+R)-p->T-p->d2);
	}
	return fit;
}


int func_ode(double t, const double y[], double f[], void *p_)
{
	struct params_rec *p = p_;
	double R = y[0];
	double N1 = y[1];
	double N2 = y[2];
	f[0] = R*(p->r*(1.0-R/p->K)-p->sm*N1*p->Imax1/(p->H1+R)-p->sm*N2*p->Imax2/(p->H2+R));
	f[1] = N1*(p->sigma1*R*p->Imax1/(p->H1+R)-p->T-p->d1);
	f[2] = N2*(p->sigma2*R*p->Imax2/(p->H2+R)-p->T-p->d2);
	return GSL_SUCCESS;
}

void simulate_ODE(struct params_rec *p, double *y, unsigned int Tf, int output)
{
	gsl_odeiv2_system sys = {func_ode, (void *)0, 3, (void *)p};
	gsl_odeiv2_driver *d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4,
					      1e-8, 1e-8, 0.0);

	double t = 0.0;
	double f[3];

	func_ode(t,y,f,(void *)p);
	if (output) {
		printf("%f %f %f %f\n", 0.0, y[0], y[1]*p->sm, y[2]*p->sm);
	}

	for (int i = 1; i <= Tf; i++) {
		int status = gsl_odeiv2_driver_apply (d, &t, i, y);

		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}
		if (output) {
			printf("%f %f %f %f\n", (double)i, y[0], y[1]*p->sm, y[2]*p->sm);
		}
	}
	gsl_odeiv2_driver_free (d);
}
