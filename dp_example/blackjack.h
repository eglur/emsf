#ifndef BLACKJACK_H
#define BLACKJACK_H


#include "util.h"

using namespace util;

namespace blackjack
{
  typedef mat stoch_mat;
  typedef std::vector<stoch_mat> v_stoch_mat;

  class model_bj
  {
  public:
    v_stoch_mat P;
    vec mu;
    vecn pi;
  };

  class data_bj
  {
  public:
    std::vector<Natural> y;
    std::vector<Natural> a;
    std::vector<Natural> r;
  };

  typedef std::vector<data_bj> v_data_bj;
}


#endif
