#ifndef EMSF_H
#define EMSF_H


#include "util.h"

using namespace util;

namespace emsf
{
  typedef mat stoch_mat;
  typedef std::vector<stoch_mat> v_stoch_mat;

  class model
  {
    public:
      v_stoch_mat P;
      vec mu;
      stoch_mat pi;
  };

  class data
  {
    public:
      vecn y;
      vecn a;
  };

  typedef std::vector<data> v_data;
  typedef std::vector<mat> v_mat;


  Natural sample_from_dist(vec dist)
  {
    Real v = random_Real();
    Real vv = 0.0;
    Natural ind = 0;

    while (1) {
      vv += dist[ind];
      if (vv > v || abs(vv - v) < 1e-5)
        break;
      else
        ++ind;
    }

    return ind;
  }


  inline void normalize(mat &A)
  {
    for (Natural i = 0; i < A.rows(); ++i) {
      Real sum = A.row(i).sum();
      if (sum > 0) for (Natural j = 0; j < A.cols(); ++j) A(i,j) /= sum;
    }
  }


  stoch_mat generate_stochastic_matrix(const Natural nrows, const Natural ncols, const bool constant = false, const bool one = false)
  {
    stoch_mat M;
    if (one) {
      M = mat::Zero(nrows, ncols);
      for (Natural i = 0; i < nrows; ++i) {
        M(i, util::random_Natural(0, ncols-1)) = 1.0;
      }
    }
    else {
      if (constant)
        M = mat::Constant(nrows,ncols,1);
      else
        M = mat::Random(nrows,ncols).cwiseAbs();
      normalize(M);
    }

    return M;
  }
}


#endif
