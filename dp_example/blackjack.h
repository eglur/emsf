#ifndef BLACKJACK_H
#define BLACKJACK_H


#include "util.h"

using namespace util;

namespace blackjack
{
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
