// Viscous burgers equation in Fourier space

#include <iostream>
#include <chrono>
#include <ctime>
#include <fftw3.h>
#include <stdio.h>
#include <complex>

#include "utilities.h" 

using namespace std;

int main () {
	
  	/* Parameters*/
	double dt, nu, Uc, p;
  	int N, Nf, Nc, Nt;
	//do 300
  	p = 5./3.; nu = 0.01; Uc = 1.; dt = 1.e-3;
  	Nf = 128; N = 2*(Nf-1); Nc = 1; Nt = 1;

	std::complex<double> alpha(dt,0.);
	
  	/* Structures */
  	std::complex<double> * UF  = new std::complex<double>[Nf];
	std::complex<double> * RHS = new std::complex<double>[Nf];
  	double *               U   = new double[N];
	double *               X   = new double[N];
		
	float t = 0.;
	
    // planner specific to 1D DFT of real data (real2complex)
	fftw_plan plan_fwd = build_plan_forward(N,U,UF);
	fftw_plan plan_inv = build_plan_inverse(N,UF,U);
	
    /* Set the intial conditions */
	set_domain(X,N);
	init_conditions(U, X, UF, plan_fwd, N, Nc);
	
	write_complex("Data/fourier.dat",UF,Nf);
	write_real("Data/physical.dat",U,X,N);
	
	/* Iteration loop */
	int i, k;
	for (i = 0; i < Nt; i++){
	
		show_modes(5, UF);
		/* Explicit Euler */
		
	    assemble_rhs(RHS, UF, Nf, nu);
		add(UF, alpha, RHS, Nf);
		
		fftw_execute(plan_inv);
		//show_modes(5, RHS);
		cout << i << " - UF Norm " << norm(UF,Nf) << " - RHS norm " << norm(RHS,Nf)  << endl;
		
		t += dt;
		
	}
	
	/* writing the physical solution */
	cout << " final time instant: tf  = " << t << " s ";
	write_real("Data/physical_final.dat",U,X,N);
	
	/* Memory deallocation */
	delete[] U; delete[] X; 
	delete[] UF; delete[] RHS;
	fftw_destroy_plan(plan_fwd); fftw_destroy_plan(plan_inv);

  	return 0;
}

/* ------ Thrash -------*/

	//std::complex<double>   UF1[Nf], UF2[Nf], UF3[Nf], UF_temp[Nf];
	//double                 U1[N], U2[N], U3[N], U_temp[N];

			/* Runge Kutta Jameson style */
		
		/*alph = dt/4.;
		
		copy(UF1,U1,UF,U,Nf,N);
        assemble_rhs(RHS, UF, U, Nf, N, nu, plan_fwd);
		add(UF1, &alph, RHS, Nf);
		fft_inverse(U1, UF1, plan_inv);
		
		alph = dt/3.;
		
		copy(UF2,U2,UF,U,Nf,N);
		assemble_rhs(RHS, UF1, U1, Nf, N, nu, plan_fwd);
		add(UF2, &alph, RHS, Nf);
		fft_inverse(U2, UF2, plan_inv);
		
		alph = dt/2.;
		
		copy(UF3,U3,UF,U,Nf,N);
		assemble_rhs(RHS, UF2, U2, Nf, N, nu, plan_fwd);
		add(UF3, &alph, RHS, Nf);
		fft_inverse(U3, UF3, plan_inv);
		
		alph = dt;
		
		assemble_rhs(RHS, UF3, U3, Nf, N, nu, plan_fwd);
		add(UF, &alph, RHS, Nf);
		fft_inverse(U, UF, plan_inv);*/
		
		/* RK4 Alternative - only one temporary thing*/
		
	    /*alph = dt/4.;
		
		copy(UF_temp,U_temp,UF,U,Nf,N);
        assemble_rhs(RHS, UF_temp, U_temp, Nf, N, nu, plan_fwd);
		add(UF_temp, &alph, RHS, Nf);
		fft_inverse(U_temp, UF_temp, plan_inv);
		
		alph = dt/3.;
		
		assemble_rhs(RHS, UF_temp, U_temp, Nf, N, nu, plan_fwd);
		copy(UF_temp,U_temp,UF,U,Nf,N);
		add(UF_temp, &alph, RHS, Nf);
		fft_inverse(U_temp, UF_temp, plan_inv);
		
		alph = dt/2.;
		
		assemble_rhs(RHS, UF_temp, U_temp, Nf, N, nu, plan_fwd);
		copy(UF_temp,U_temp,UF,U,Nf,N);
		add(UF_temp, &alph, RHS, Nf);
		fft_inverse(U_temp, UF_temp, plan_inv);
		
		alph = dt;
		
		assemble_rhs(RHS, UF_temp, U_temp, Nf, N, nu, plan_fwd);
		add(UF, &alph, RHS, Nf);
		fft_inverse(U, UF, plan_inv);*/
	
/*	end = std::chrono::system_clock::now();
	
	chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = chrono::system_clock::to_time_t(end);
 
    cout << "finished computation at " << ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
			  
			      /* Measuring CPU TIME*/
	/*chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();*/