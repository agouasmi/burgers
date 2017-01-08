#include "utilities.h"

using namespace std;

fftw_plan build_plan_forward(int N, double* U, std::complex<double>* UF){
	
    fftw_plan plan_fwd = fftw_plan_dft_r2c_1d(N, 
									U,
									reinterpret_cast<fftw_complex*> (UF),
									FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
	fftw_forget_wisdom();
	
	return plan_fwd;
	
}

fftw_plan build_plan_inverse(int N, std::complex<double>* UF, double* U){
	
    fftw_plan plan_inv = fftw_plan_dft_c2r_1d(N, 
									reinterpret_cast<fftw_complex*> (UF), 
									U, 
									FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);	
	fftw_forget_wisdom();
	
	return plan_inv;
	
}

void fft_forward(std::complex<double> * UF, double * U, int N, fftw_plan plan_fwd){
	/* Carry out the fftw plan on U to produce UF*/
	
	static std::complex<double> ZERO(0.,0.);
    fftw_execute_dft_r2c(plan_fwd, U, reinterpret_cast<fftw_complex*> (UF));
	
    static std::complex<double> fact(1./N,0.);
	
	/* Necessary scaling - not needed for the inverse*/
	scale(UF, int(N/2+1), &fact);
	UF[N/2] = ZERO;
}

void fft_inverse(double * U, std::complex<double> * UF, int N, fftw_plan plan_inv){
	
    fftw_execute_dft_c2r(plan_inv,reinterpret_cast<fftw_complex*> (UF),U);

    //static double fact = 1;	
	//scale(U, N, fact);
}

void square(double* U, int N){
	int i;
	for (i = 0; i < N; i++){
		U[i] *= U[i];	
	}
}

void square(double* U2, double* U, int N){
	int i;
	for (i = 0; i < N; i++){
		U2[i] = pow(U[i],2);	
	}
}

void convo_23(std::complex<double> * target, std::complex<double> * UF, int Nf){

	int N = 2*(Nf-1);
	
	std::complex<double> ZERO(0.,0.);
	std::complex<double> UF_23[Nf];
	double U_23[N];
	
	copy_complex(UF_23,UF,Nf);
	
	int N_23 = 2*Nf/3+1;
	
	/* remove the last 1/3 of the information to avoid aliasing*/
	int i;
	for (i = N_23; i < Nf; i++){
		UF_23[i] = ZERO;
	}
	
	fftw_plan plan_fwd    = build_plan_forward(N, U_23, UF_23);
	fftw_plan plan_inv    = build_plan_inverse(N, UF_23, U_23);
	
	/* compute the corresponding physical signal, don't need to scale */
	fft_inverse(U_23,UF_23,N,plan_inv);
	
	/* square it */
	square(U_23,N); 
	
	/* compute the Fourier transform of that and assign it to target*/
	fft_forward(target,U_23,N,plan_fwd);
	
	fftw_destroy_plan(plan_fwd); fftw_destroy_plan(plan_inv);
}


void convo(std::complex<double> * target, std::complex<double> * UF, int Nf){
    /* Computes the convolution and puts the result in target */
	/* If you always apply it to the same structure in your code, it's gonna be fine. otherwise...*/
	
	int i;
	
	/* Pad the damn fourier vector*/
	int Nf_pad = int (Nf);
	int N_pad  = 2*(Nf_pad-1);
	
	/* Temporary structures */
	std::complex<double> UF_pad[Nf_pad];
	double               U_pad[N_pad];
	
	fftw_plan plan_pad_fwd    = build_plan_forward(N_pad,   U_pad,   UF_pad);
	fftw_plan plan_pad_inv    = build_plan_inverse(N_pad,   UF_pad,   U_pad);
	
	/* Do the padding... basically copy UF into a larger one full of zeros*/
	copy_complex(UF_pad,UF,Nf);	

	/* compute the corresponding physical signal, don't need to scale */
	fft_inverse(U_pad,UF_pad,N_pad,plan_pad_inv);
	
	/* square it */
	square(U_pad,N_pad); 
	
	/* compute the Fourier transform of that squared signal*/
	fft_forward(UF_pad,U_pad,N_pad,plan_pad_fwd);
	
	/* put it into target*/
	copy_complex(target,UF_pad,Nf);	
	//cout << "convolution: \n";
	//show_modes(5,target);
	
	fftw_destroy_plan(plan_pad_fwd); fftw_destroy_plan(plan_pad_inv);
	
}


void assemble_rhs(std::complex<double> * RHS, std::complex<double> * UF, int Nf, double nu){
	/* Assembles the right hand side term of Burgers in Fourier Space */
  
	/* First put in the convection term*/
	convo_23(RHS, UF, Nf);
	
	//show_modes(5,RHS);
	
	int k;
	for(k = 0; k < Nf; k++){
		std::complex<double> ik2( 0. , -double(k)/2. );
		RHS[k] *= ik2;	
	}
	
	/* Add the Linear term */
	
	for (k = 0; k < Nf; k++){
		double a = -nu*pow((double) k, 2);
		RHS[k] +=  a * UF[k];
	}
	
}


void scale(std::complex<double>* signal, const int N, const void* fact){
    /* Uses BLAS to perform signal = signal * fact */
	
    cblas_zscal(N, fact, signal, 1);
 
} 

void scale(double* signal, const int N, double fact){
    /* Uses BLAS to perform signal = signal * fact */
	
    cblas_dscal(N, fact, signal, 1);
 
} 

/*void add(std::complex<double> * UF, std::complex<double> dt, std::complex<double> * RHS, int Nf){
	 Uses BLAS to perform UF = UF + dt.RHS 
	
	cblas_zaxpy(Nf, dt, RHS, 1, UF, 1);
	
}*/

void add(std::complex<double> * UF, std::complex<double> dt, std::complex<double> * RHS, int Nf){
	int i;
	for (i = 0; i < Nf; i++){
		UF[i] += dt*RHS[i];
	}
	
}


void copy_complex(std::complex<double>* target, std::complex<double>* source, int size){
    
	cblas_zcopy(size,source,1,target,1); 
	
}

void copy_real(double * target, double * source, int size){

	cblas_dcopy(size,source,1,target,1);
	
}

double norm(double * U, int N){
	return cblas_dnrm2(N, U, 1);	
}

double norm(std::complex<double> * UF, int Nf){
	return cblas_dznrm2(Nf, UF, 1);
}

void set_domain(double* X, int N){
	/* [0  2 M_PI] */
	
	double theta = 1. / ((double)N) * 2.0 * M_PI;
	int i;
	for (i = 0; i < N; ++i) {       	
        	X[i] = (double) i * theta;
    }
}

void init_conditions(double * U, double * X, std::complex<double> * UF, fftw_plan plan, int N, int Nc){
    /* Generate a real signal on [0 2Pi] */
	
    int i, k;
	for (k = 1; k < Nc+1; k++) {	
    	for (i = 0; i < N; ++i) {       	
        	U[i] += 2. * sin((double) k * X[i]);
    	}
	}
	
	fft_forward(UF, U, N, plan);
}

void show_modes(int N, std::complex<double>* UF){
	int i;
	cout << " ";
	for (i = 0; i < N-1; i++){
	   double r  = (abs(UF[i].real())<1e-8)? 0 : UF[i].real();
	   double im = (abs(UF[i].imag())<1e-8)? 0 : UF[i].imag();
	   cout << " Mode " << i << ": real " << r << ", imag " << im << " \n ";	
	}
	cout << " Mode " << i << ": real " << 0 << ", imag " << 0 << "      -  Oddball wavenumber \n ";
	cout << "\n\n";
	
}

void show_signal(int N, double* U){
	int i;
	cout << " ";
	for (i = 0; i < N; i++){
	   double r  = (U[i]<1e-8)? 0 : U[i];
	   cout << r << " \n ";	
	}
	cout << "\n";
	
}

void write_complex(string filename, std::complex<double>* data, int N){
	ofstream fout;
	fout.open(filename);
	for(int k = 0; k < N; k++){
		fout << k << "  " << data[k].real() << "  " << data[k].imag() << "\n";
	}
	fout.close();
	
}

void write_real(string filename, double* data, double* X, int N){
	
	ofstream fout;
	fout.open(filename);
	for(int k = 0; k < N; k++){
		fout << X[k] << "  " << data[k] << "\n";
	}
	fout.close();
}


/* -------------------- THRASH ------------------------ */
