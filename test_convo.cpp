/* Purpose: Doing the convolution through fft/ifft correctly */

/* Completed: */

#include <iostream>
#include <chrono>
#include <ctime>
#include <fftw3.h>
#include <stdio.h>
#include <complex>

#include "utilities.h" 

using namespace std;

int main () {
	
	static std::complex<double> ZERO(0.,0.);
	static std::complex<double> I(0.,1.);
	static std::complex<double> R(1.,0.);
	
    /* Build the test vectors */
  	int N, Nf;

  	Nf = 256; N = 2*(Nf-1); 
	
  	/* Two vectors - Theirs convolutions should be the same */
  	std::complex<double> UF[Nf];
	std::complex<double> UF_1[Nf];
	std::complex<double> conv[Nf]; 
	
	/* Initialize the data */
	// for (int i = 0; i<5;i++) {UF[i] = R + double(i)*I; UF_1[i] = double(i)*R - I ;} 
 	UF[1] = -I;
	
	/* Perform the 1st convolution */
	show_modes(5,UF);
	convo(conv, UF, Nf);
	
	cout << " /* the first convolution */ \n \n";
	
	/* show the results */
	//show_modes(5, conv); 
	
	/*----*/
	
	/* Perform the 2nd convolution */
	//show_modes(Nf,UF_1);
	//convo(conv, UF_1, Nf);
	
	//cout << " /* the second convolution */ \n \n";
	
	/* show the results */
	//show_modes(Nf, conv); 
	
  	return 0;
}

