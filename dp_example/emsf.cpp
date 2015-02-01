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


stoch_mats generate_stochastic_matrices(const Natural nrows,const Natural ncols, const Natural na)
{
  stoch_mats MS;
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

  while (vv < v) {
    ++ind;
    vv += dist[ind];
  }

  return ind;
}


int main()
{
  const Natural na = 3;
  const Natural nrows = 3;
  const Natural ncols = 3;

  stoch_mat m = generate_stochastic_matrix(nrows, ncols);
  stoch_mats ms = generate_stochastic_matrices(na, nrows, ncols);

  for (Natural i = 0; i < na; ++i)
    cout << "m[" << i << "] =" << endl << ms[i] << endl;

  vec test = m.col(1);
  cout << "test =" << endl << test << endl;

  model my_model;
  my_model.P = ms;

  for (Natural i = 0; i < na; ++i)
    cout << "my_model.P[" << i << "] =" << endl << my_model.P[i] << endl;

  return 0;
}
