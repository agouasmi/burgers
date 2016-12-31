// Class for the problem configuration
#include <iostream>
using namespace std;

class Config {
  public:
    float p, nu, Uc, dt;
    int N, Nf, Nc;
 
    Config (float, float, float,
               float, int, int, int);
};

Config::Config (float p_, float nu_, float Uc_,
                float dt_, int N_, int Nf_, int Nc_) {
  p = p_;
  nu = nu_;
  Uc = Uc_;
  dt = dt_;
  N  = N_;
  Nf = Nf_;
  Nc = Nc_;
}

