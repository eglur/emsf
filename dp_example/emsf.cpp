#include <iostream>
#include "eigen3/Eigen/Dense"
#include "util.h"

using namespace std;
using namespace Eigen;
using namespace util;

typedef mat stoch_mat;
typedef std::vector<stoch_mat> stoch_mats;


stoch_mat generate_stochastic_matrix(const Natural nrows, const Natural ncols)
{
  stoch_mat M = mat::Random(nrows,ncols).cwiseAbs();
  vec result = M.rowwise().sum();
  for(Natural i = 0; i < ncols; ++i) {
    M.col(i) = M.col(i).array() / result.array();
  }

  return M;
}


stoch_mats generate_stochastic_matrices(const Natural na, const Natural nrows,
                                        const Natural ncols)
{
  stoch_mats MS;
  MS.resize(na);
  for(Natural a = 0; a < na; ++a)
    MS[a] = generate_stochastic_matrix(nrows,ncols);

  return MS;
}


int sample_from_dist(MatrixXd dist)
{
  double v = 0.0;                                              // TODO
  double vv = 0.0;
  int ind = 0;
  while (vv < v) {
    ind++;
    vv = vv + dist[ind];
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

  return 0;
}
