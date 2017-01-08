// Viscous burgers equation in Fourier space

#include <iostream>
#include <chrono>
#include <ctime>
#include <fftw3.h>
#include <stdio.h>
#include <complex>

#include "utilities.h" 

using namespace std;

int main (int argc, char *argv[]) {
	
 	/* Parameters*/
	double dt, nu;
  	int N, Nf, Nc, Nt, RK;
	
	/* Reading them from the file path provided in argv */
	cout << argv[1] << endl;
	std::fstream infile(argv[1], std::ios_base::in);
	infile >> dt >> nu >> Nf >> Nc >> Nt >> RK;
	infile.close();
	
	float t = 0.;
	N = 2*(Nf-1); 	
	std::complex<double> alpha(dt,0.);
	
  	/* Base Structures */
  	std::complex<double> * UF  = new std::complex<double>[Nf];
	std::complex<double> * RHS = new std::complex<double>[Nf];
  	double *               U   = new double[N];
	double *               X   = new double[N];
	
	/* Temporary Fourier vector for RK4*/
	std::complex<double> * UF_temp  = new std::complex<double>[Nf];
	
    // planner specific to 1D DFT of real data (real2complex)
	fftw_plan plan_fwd = build_plan_forward(N,U,UF);
	fftw_plan plan_inv = build_plan_inverse(N,UF,U);
	
    /* Set the intial conditions */
	set_domain(X,N);
	init_conditions(U, X, UF, plan_fwd, N, Nc);
	
	/* Write the initial solution */
	write_complex("Data/fourier_init.dat",UF,Nf);
	write_real("Data/physical_init.dat",U,X,N);
	
	/* Iteration loop */
	int i, k;
	for (i = 0; i < Nt; i++){
			
		copy_complex(UF_temp,UF,Nf);
		/* Runge-Kutta Schemes - only one temporary thing*/
			
		for (k = RK; k > 1; k--){
			
			alpha = dt/double(k);
			assemble_rhs(RHS, UF_temp, Nf, nu);
			copy_complex(UF_temp,UF,Nf);
			add(UF_temp, alpha, RHS, Nf);
			
		}
		
		alpha = dt;
	    assemble_rhs(RHS, UF_temp, Nf, nu);
		add(UF, alpha, RHS, Nf);
		
		/* ----- */
		
		cout << i << " - UF Norm " << norm(UF,Nf) << " - RHS norm " << norm(RHS,Nf)  << endl;
		
		t += dt;
		
	}
	
	/* to recover the physical signal*/
	fftw_execute(plan_inv);
	
	/* writing the physical solution */
	cout << " final time instant: tf  = " << t << " s ";
	write_real("Data/physical_final.dat",U,X,N);
	
	/* Memory deallocation */
	delete[] U; delete[] X; delete[] UF_temp;
	delete[] UF; delete[] RHS;
	fftw_destroy_plan(plan_fwd); fftw_destroy_plan(plan_inv);

  	return 0;
}
	
/*	
    /* Measuring CPU TIME
	chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	end = std::chrono::system_clock::now();
	
	chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = chrono::system_clock::to_time_t(end);
 
    cout << "finished computation at " << ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
			  
			      */