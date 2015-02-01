#include <iostream>
#include <cmath>
#include "eigen3/Eigen/Dense"
#include "util.h"
#include "emsf.h"


using namespace std;
using namespace Eigen;
using namespace util;
using namespace emsf;


inline mat cwiseInverse(mat A)
{
  mat B(A.rows(), A.cols());
  for (Natural i = 0; i < A.rows(); ++i)
    for (Natural j = 0; j < A.cols(); ++j)
      B(i, j) = 1.0 / A(i, j);

  return B;
}


inline mat cwiseLog(mat A)
{
  mat B(A.rows(), A.cols());
  for (Natural i = 0; i < A.rows(); ++i)
    for (Natural j = 0; j < A.cols(); ++j)
      B(i, j) = log(A(i, j));

  return B;
}


inline void normalize(mat &A)
{
  for (Natural i = 0; i < A.rows(); ++i) {
    Real sum = A.row(i).sum();
    if (sum > 0) for (Natural j = 0; j < A.cols(); ++j) A(i,j) /= sum;
  }
}


stoch_mat generate_stochastic_matrix(const Natural nrows, const Natural ncols)
{
  stoch_mat M = mat::Random(nrows,ncols).cwiseAbs();
  normalize(M);

  return M;
}


v_stoch_mat generate_stochastic_matrices(const Natural nrows,const Natural ncols, const Natural na)
{
  v_stoch_mat MS;
  MS.resize(na);
  for(Natural a = 0; a < na; ++a)
    MS[a] = generate_stochastic_matrix(nrows,ncols);

  return MS;
}


v_mat generate_zero_matrices(const Natural nrows, const Natural ncols, const Natural na)
{
  v_mat vm;
  vm.resize(na);
  for (Natural a = 0; a < na; ++a)
    vm[a] = mat::Zero(nrows, ncols);

  return vm;
}


Natural sample_from_dist(vec dist)
{
  Real v = random_Real();
  Real vv = 0.0;
  Natural ind = 0;

  while (1) {
    vv += dist[ind];
    if (vv > v)
      break;
    else
      ++ind;
  }

  return ind;
}


model generate_model(const Natural n, const Natural m, const Natural na)
{
  v_stoch_mat D = generate_stochastic_matrices(n, m, na);
  v_stoch_mat K = generate_stochastic_matrices(m, n, na);

  model md;
  md.P.resize(na);
  for (Natural a = 0; a < na; ++a)
    md.P[a] = D[a] * K[a];

  md.mu = generate_stochastic_matrix(1, na).transpose();
  md.pi = generate_stochastic_matrix(n, na);

  return md;
}


data generate_data(model md, const Natural T)
{
  data dt;
  dt.y.resize(T);
  dt.a.resize(T-1);

  dt.y[0] = sample_from_dist(md.mu);
  for (Natural i = 0; i < T-1; ++i) {
    dt.a[i] = sample_from_dist(md.pi.row(dt.y[i]).transpose());
    dt.y[i+1] = sample_from_dist(md.P[dt.a[i]].row(dt.y[i]).transpose());
  }

  return dt;
}


v_data generate_batch_data(model md, const Natural T, const Natural num_batches)
{
  v_data v_dt;
  v_dt.resize(num_batches);

  for (Natural i = 0; i < num_batches; ++i)
    v_dt[i] = generate_data(md, T);

  return v_dt;
}


void em_sf(model md, v_data dt, const Natural n, const Natural m, const Natural na, const Real eps = 1e-20, const Natural max_it = 10)
{
  v_stoch_mat P = md.P;
  vec mu = md.mu;
  stoch_mat pi = md.pi;

  v_stoch_mat D = generate_stochastic_matrices(n, m, na);
  v_stoch_mat K = generate_stochastic_matrices(m, n, na);

  Real old_score = minf;
  Real score = 0.0;
  Natural it = 0;

  std::vector<Real> error;

  while (abs(score - old_score) > eps || it < max_it) {
    old_score = score;
    score = 0.0;

    v_mat D2 = generate_zero_matrices(n, m, na);
    v_mat K2 = generate_zero_matrices(m, n, na);

    for (Natural batch = 0; batch < dt.size(); ++batch) {
      vecn y = dt[batch].y;
      vecn a = dt[batch].a;
      
      Natural T = y.cols();

      mat A = mat::Zero(T-1, m);
      mat B = mat::Zero(T-1, m);
      vec NF = vec::Zero(T-1);

      A.row(0) = mu(y(0)) * pi(y(0), a(0)) * D[a(0)].row(y(0));
      NF(0) = A.row(0).sum();
      A.row(0) = A.row(0) / NF(0);

      for (Natural t = 1; t < T-1; ++t) {
        row_vec k_row = K[a(t-1)].col(y(t)).transpose();
        Real aa = A.row(t-1).cwiseProduct(k_row).sum();

        A.row(t) = aa * pi(y(t), a(t)) * D[a(t)].row(y(t));
        NF(t) = A.row(t).sum();
        A.row(t) = A.row(t) / NF(t);
      }

      Natural B_T = T-2;                                        // Different index strategie on the R language
      B.row(B_T) = K[a(B_T)].col(y(B_T)).transpose();
      B.row(B_T) = B.row(B_T) / NF(B_T);
      for (Natural t = B_T-1; t >= 0; --t) {
        row_vec d_row = D[a(t+1)].row(y(t+1));
        Real bb = B.row(t+1).cwiseProduct(d_row).sum();

        row_vec k_row = K[a(t)].col(y(t+1)).transpose();
        B.row(t) = bb * pi(y(t+1), a(t+1)) * k_row;
        B.row(t) = B.row(t) / NF(t);
      }

      mat C = A.cwiseProduct(B);
      normalize(C);

      for (Natural t = 0; t < T-1; ++t) {
        D2[a(t)].row(y(t)) = D2[a(t)].row(y(t)).array() + C.row(t).array();
        K2[a(t)].col(y(t+1)) = K2[a(t)].col(y(t+1)).array() + C.row(t).transpose().array();
      }

      score = score - cwiseLog(cwiseInverse(NF)).sum();
    }

    for (Natural i = 0; i < na; ++i) {
      normalize(D2[i]);
      normalize(K2[i]);
    }

    D = D2;
    K = K2;

    Real e = 0.0;
    for (Natural i = 0; i < na; ++i) {
      mat P_bar = D[i] * K[i];
      mat diff = P_bar.array() - P[i].array();
      mat sqr_diff = diff.array() * diff.array();
      e += sqr_diff.sum();
    }

    error.push_back(e);

    cout << "e: " << e << endl;
    cout << "score: " << score << endl;
  }
}


int main()
{
  const Natural n = 50;
  const Natural m = 5;
  const Natural na = 1;
  const Natural T = 5000;
  const Natural num_batches = 10;

  model md = generate_model(n, m, na);
  v_data dt = generate_batch_data(md, T, num_batches);
  em_sf(md, dt, n, m, na);
  
  return 0;
}
