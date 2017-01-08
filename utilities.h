#include <iostream>
#include <fftw3.h>
#include <stdio.h>
#include <string>
#include <cblas.h>
#include <complex>
#include <fstream>
#include <string>

/* FFT functions */
fftw_plan build_plan_forward(int N, double* U, std::complex<double>* UF);
fftw_plan build_plan_inverse(int N, std::complex<double>* UF, double* U);
void fft_forward(std::complex<double> * UF, double * U, int N, fftw_plan plan);
void fft_inverse(double * U, std::complex<double> * UF, int N, fftw_plan plan_inv);

/* Vectors operations (using BLAS, FFTW) */
void copy_complex(std::complex<double>* target, std::complex<double>* source, int size);
void copy_real(double * target, double * source, int size);
void square(double* U2, double* U, int N);
void scale(std::complex<double>* signal, const int length, const void* fact);
void scale(double* signal, const int N, double fact);
void add(std::complex<double> * UF, std::complex<double> dt, std::complex<double> * RHS, int Nf);
void convo(std::complex<double> * target, std::complex<double> * UF, int Nf);
void convo_23(std::complex<double> * target, std::complex<double> * UF, int Nf);
double norm(double* U, int N);
double norm(std::complex<double> * UF, int Nf);

/* Problem specific functions */
void set_domain(double* X, int N);
void assemble_rhs(std::complex<double> * RHS, std::complex<double> * UF, int Nf, double nu);
void init_conditions(double * U, double* X, std::complex<double> * UF, fftw_plan plan, int N, int Nc);

/* IO routines */
void show_modes(int N, std::complex<double>* UF);
void show_signal(int N, double* U);

/* Write data to file */
void write_complex(std::string filename, std::complex<double>* data, int N);
void write_real(std::string filename, double* data, double* x, int N);