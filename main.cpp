// Viscous burgers equation
#include <fftw3.h>
#include <rfftw.h>
#include <iostream>
#include <stdio.h>
#include <string>

#include "config.h" 

// 

using namespace std;

int main () {

  /* Parameters*/
  float p, nu, Uc, dt;
  int N, Nf, Nc;

  p = 5./3.; nu = 1.e-3; Uc = 1.; dt = 1.e-3;
  Nf = 128; N = 2*Nf; Nc = 32;

  Config conf = Config(p,nu,Uc,dt,Nf,N,Nc);

  cout << "Heyyyy!";

  return 0;
}

/*int main()
{
  
  // Problem parameters
  const float p, nu, Uc, dt;
  const int N, Nf, Nc;

  p = 5./3.; nu = 1.e-3; Uc = 1.; dt = 1.e-3;
  Nf = 128; N = 2*Nf; Nc = 32;

  Config conf = Config(p,nu,Uc,dt,Nf,N,Nc);
   
  cout << result;

  // process:
  a = 5;
  b = 2;
  a = a + 1;
  result = a - b;

  // Compound assignment
  a += b;
  a -= b;
  a /= b;
  a *= (b+1);

  // Conditional ternary operator
  c = (a>b) ? a : b; // evaluates to a if (a>b), b otherwise

  // print out the result:
  cout << result;

  cout << "Hello World! \n \n";
  cout << "I'm a C++ program \n \n";

  // string manipulation
  string mystring, x;
  mystring = "This is a string \n";
  cout << mystring;
  x = "string expressed in \
  two lines";
  cout << x;

  // terminate the program
  return 0;
}*/
