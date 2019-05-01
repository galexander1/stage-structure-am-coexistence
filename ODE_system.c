
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

//compile with:
// gcc -lgsl -lm -o O ODE_system.c -L /usr/local/lib -I /usr/local/include/

int func_ode(double t, const double y[], double f[], void *params);
void simulate_ODE(double *parms, double *y0, unsigned int Tf, int output);
double fitness(double *params, double R, int species);
void linspace(double min, double max, int points, double *c);

int main() {
	FILE *file = fopen("CASE0_PDE.dat","w");
	double *parms = (double*)calloc(12,sizeof(double));
	parms[0] = 1.0; //r
	parms[1] = 6.0; //K
	parms[2] = 0.5; //sigma1
	parms[3] = 0.5; //sigma2
	parms[4] = 3.0; //H1
	parms[5] = 2.5; //H2
	parms[6] = 1.0; //Imax1
	parms[7] = 1.1; //Imax2
	parms[8] = 0.015; //d1 
	parms[9] = 0.04; //d2
	parms[10] = 0.1; //T
	parms[11] = 1e-4; //sm
	
	double y0[3] = {1,1,1};
	simulate_ODE(parms, y0, 10000, /*output=*/0);
	simulate_ODE(parms, y0, 40000, /*output=*/0);
	simulate_ODE(parms, y0, 1000, /*output=*/1);
	
	
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
				fprintf(file,"%f %f %f\n",R, fitness(parms, R, 1), fitness(parms, R, 2));
			} else {
				fprintf(file,"%f %f %f\n",R, fitness(parms, R, 1), nan(""));
			}
		} else {
			if (R >= Rmin2 && R <= Rmax2) {
				fprintf(file,"%f %f %f\n",R, nan(""), fitness(parms, R, 2));
			} else {
				fprintf(file,"%f %f %f\n",R, nan(""), nan(""));
			}
		}
	}
	fclose(file);
	
	return 0;
}


double fitness(double *params, double R, int species) {
	double r = params[0];
	double K = params[1];
	double sigma1 = params[2];
	double sigma2 = params[3];
	double H1 = params[4];
	double H2 = params[5];
	double Imax1 = params[6];
	double Imax2 = params[7];
	double d1 = params[8];
	double d2 = params[9];
	double T = params[10];
	double sm = params[11];
	
	double fit = 0.0;
	if (species == 1) {
		fit = 0.05*(sigma1*R*Imax1/(H1+R)-T-d1);
	} else {
		fit = 0.05*(sigma2*R*Imax2/(H2+R)-T-d2);
	}
    return fit;
	
}


int func_ode(double t, const double y[], double f[], void *params_) {
	double *params = params_;
	double r = params[0];
	double K = params[1];
	double sigma1 = params[2];
	double sigma2 = params[3];
	double H1 = params[4];
	double H2 = params[5];
	double Imax1 = params[6];
	double Imax2 = params[7];
	double d1 = params[8];
	double d2 = params[9];
	double T = params[10];
	double sm = params[11];
	
	double R = y[0];
	double N1 = y[1];
	double N2 = y[2];
	f[0] = R*(r*(1.0-R/K)-sm*N1*Imax1/(H1+R)-sm*N2*Imax2/(H2+R));
	f[1] = N1*(sigma1*R*Imax1/(H1+R)-T-d1);
	f[2] = N2*(sigma2*R*Imax2/(H2+R)-T-d2);
    return GSL_SUCCESS;
}


void simulate_ODE(double *parms, double *y0, unsigned int Tf, int output) {
	
    gsl_odeiv2_system sys = {func_ode, (void *)0, 3, (void *)parms};
    gsl_odeiv2_driver *d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,1e-8, 1e-8, 0.0);
	
	double t = 0.0;
    double y[3];
	double f[3];
	y[0] = y0[0];
	y[1] = y0[1];
	y[2] = y0[2];
	
	func_ode(t,y,f,(void *)parms);
    if (output) {
       printf("%f %f %f %f\n", 0.0, y[0], y[1]*parms[11], y[2]*parms[11]);
    }

    for (int i = 1; i <= Tf; i++)
    {
        double ti = i;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        
        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }
	if (output) {
		printf("%f %f %f %f\n", ti, y[0], y[1]*parms[11], y[2]*parms[11]);
	}
    }
    gsl_odeiv2_driver_free (d);
}
