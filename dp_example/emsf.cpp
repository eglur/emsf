#include <iostream>
#include "eigen3/Eigen/Dense"
#include "util.h"
#include "emsf.h"


using namespace std;
using namespace Eigen;
using namespace util;
using namespace emsf;


stoch_mat generate_stochastic_matrix(const Natural nrows, const Natural ncols)
{
  stoch_mat M = mat::Random(nrows,ncols).cwiseAbs();
  vec result = M.rowwise().sum();
  for(Natural i = 0; i < ncols; ++i) {
    M.col(i) = M.col(i).array() / result.array();
  }

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


int main()
{
  const Natural n = 3;
  const Natural m = 3;
  const Natural na = 3;

  model md = generate_model(n, m, na);
  for (Natural i = 0; i < na; ++i)
    cout << "md.P[" << i << "] =" << endl << md.P[i] << endl;
  cout << "md.mu =" << endl << md.mu << endl;
  cout << "md.pi =" << endl << md.pi << endl;

  return 0;
}
