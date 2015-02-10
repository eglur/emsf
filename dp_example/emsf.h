#ifndef EMSF_H
#define EMSF_H


#include "util.h"
#include "blackjack.h"

using namespace util;
using namespace blackjack;

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
    vecn r;
  };

  typedef std::vector<data> v_data;
  typedef std::vector<mat> v_mat;

  Natural sample_from_dist(vec dist);
  void normalize(mat &A);
  stoch_mat generate_stochastic_matrix(const Natural nrows, const Natural ncols, const bool constant = false, const bool one = false);
  v_stoch_mat generate_stochastic_matrices(const Natural nrows,const Natural ncols, const Natural na, const bool constant = false, const bool one = false);
  v_mat generate_zero_matrices(const Natural nrows, const Natural ncols, const Natural na);
  model generate_model(const Natural n, const Natural sr, const Natural na);
  void em_sf_sk(model &md, v_data_bj &dt, const Natural n, const Natural m, const Natural na, const Natural num_batches, v_stoch_mat D, v_stoch_mat K, const Natural max_it);
}


#endif
