#ifndef EMSF_H
#define EMSF_H


#include "util.h"

using namespace util;

namespace emsf
{
  typedef mat stoch_mat;
  typedef std::vector<stoch_mat> stoch_mats;

  class model
  {
    public:
      stoch_mats P;
      stoch_mat mu;
      stoch_mat pi;
  };
}


#endif
