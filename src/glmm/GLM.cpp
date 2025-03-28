#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <string>
#include <random>

using namespace std;
using namespace arma;


class Family {
public:
  virtual vec link(const vec &mu) const = 0;
  virtual vec linkinv(const vec &eta) const = 0;
  virtual vec var(const vec &mu) const = 0;
  virtual vec mu_start(const vec &y) const = 0;
  virtual string getFamilyName() const = 0;
  virtual ~Family() {}
};


class Gaussian : public Family {
public:
  vec link(const vec &mu) const override { return mu; }
  vec linkinv(const vec &eta) const override { return eta; }
  vec var(const vec &mu) const override { return ones<vec>(mu.n_elem); }
  vec mu_start(const vec &y) const override { return y; }
  string getFamilyName() const override { return "Gaussian"; }
};


class Binomial : public Family {
public:
  vec link(const vec &mu) const override { return log(mu) - log(1 - mu); }
  vec linkinv(const vec &eta) const override { return exp(eta) / (1 + exp(eta)); }
  vec var(const vec &mu) const override { return mu % (1 - mu); }
  vec mu_start(const vec &y) const override {
    if(any(y < 0) || any(y > 1))
      throw std::invalid_argument("y values must be between 0 and 1 for the Binomial family.");
    return (y + 0.5) / 2.0;
  }
  string getFamilyName() const override { return "Binomial"; }
};


class Poisson : public Family {
public:
  vec link(const vec &mu) const override { return log(mu); }
  vec linkinv(const vec &eta) const override { return exp(eta); }
  vec var(const vec &mu) const override { return mu; }
  vec mu_start(const vec &y) const override {
    if(any(y < 0))
      throw std::invalid_argument("Negative y values not allowed for the Poisson family.");
    return y + 0.1;
  }
  string getFamilyName() const override { return "Poisson"; }
};


vec FitIRWLS_GetBeta(const vec &y, const mat &X, const Family &family, double tol = 1e-5, int j_max = 25) {
  int m = X.n_cols;                     
  vec beta = zeros<vec>(m);             
  vec mu = family.mu_start(y);          
  vec eta = family.link(mu);            

  for (int j = 0; j < j_max; ++j) {
    vec v = family.var(mu);            
    vec Z = eta + (y - mu) / v;        
    mat W = diagmat(v);                
    vec beta_new = solve(X.t() * W * X, X.t() * W * Z); 
    double epsilon = norm(beta_new - beta, 2); 
    beta = beta_new;                    

    if(epsilon < tol)                  
      break;
    if(family.getFamilyName() == "Gaussian") 
      break;
    
    eta = X * beta;                    
    mu = family.linkinv(eta);          
  }
  return beta;                        
}




