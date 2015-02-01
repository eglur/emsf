#ifndef EMSF_H
#define EMSF_H


#include "util.h"

using namespace util;

namespace emsf
{
  typedef mat stoch_mat;
  typedef std::vector<stoch_mat> v_stoch_mat;
  typedef std::vector<vecn> vecns;

  class model
  {
    public:
      v_stoch_mat P;
      stoch_mat mu;
      stoch_mat pi;
  };
}


#endif
