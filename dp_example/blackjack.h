#ifndef BLACKJACK_H
#define BLACKJACK_H


#include "util.h"

using namespace util;

namespace blackjack
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
    std::vector<Natural> y;
    std::vector<Natural> a;
    std::vector<Natural> r;
  };

  typedef std::vector<data> v_data;
  typedef std::vector<mat> v_mat;

  Natural sample_from_dist(vec dist);
  void normalize(mat &A);
  stoch_mat generate_stochastic_matrix(const Natural nrows, const Natural ncols, const bool constant = false, const bool one = false);
  v_stoch_mat generate_stochastic_matrices(const Natural nrows,const Natural ncols, const Natural na, const bool constant = false, const bool one = false);
}


#endif
