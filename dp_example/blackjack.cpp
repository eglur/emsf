#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include "eigen3/Eigen/Dense"
#include "util.h"
#include "emsf.h"


using namespace std;
using namespace Eigen;
using namespace util;
using namespace emsf;


Natural card(vec &card_dist)
{
  Natural c = sample_from_dist(card_dist) + 1;
  if (c > 10) c = 10;

  return c;
}


void draw_card(Natural &xc, Natural &x_ace, vec &card_dist)
{
  Natural c = card(card_dist);

  xc += c;

  if (!x_ace && c == 1) {
    xc += 10;
    x_ace = 1;
  }

  if (x_ace && xc > 21) {
    xc -= 10;
    x_ace = 0;
  }
}


void transition(Natural &pc, Natural &p_ace, Natural &dc, Natural &d_ace, Natural &r, Natural a, bool &end, vec &card_dist)
{
  Natural stick = 0;
  Natural hit = 1;
  
  if (a == hit) {
    draw_card(pc, p_ace, card_dist);

    if (pc > 21) {
      r = -1;
      end = true;
    }
    else {
      end = false;
    }
  }
  else if (a == stick) {
    while (dc < 17)                                             /* TODO: checar se é < ou <= */
      draw_card(dc, d_ace, card_dist);

    if (dc > 21)
      r = 1;
    else {
      Natural p_diff = 21 - pc;
      Natural d_diff = 21 - dc;

      if (p_diff < d_diff)
        r = 1;
      else if (p_diff > d_diff)
        r = -1;
      else
        r = 0;
    }

    end = true;
  }
}


Natural get_s(Natural pc, Natural p_ace, Natural dc)
{
  return p_ace * 100 + pc * 10 + dc;
}


inline Natural get_a(Natural s, mat &pi)
{
  //TODO

  return rand() < 0.5;
}


Natural episode(mat &pi, vec &card_dist)
{
  Natural pc = 0;
  Natural p_ace = 0;

  Natural dc = 0;
  Natural d_ace = 0;

  Natural r;
  
  draw_card(pc, p_ace, card_dist);
  draw_card(pc, p_ace, card_dist);

  draw_card(dc, d_ace, card_dist);

  bool end = false;
  while (!end) {
    Natural s = get_s(pc, p_ace, dc);
    Natural a = get_a(s, pi);
    transition(pc, p_ace, dc, d_ace, r, a, end, card_dist);

    // cout <<
    //   pc <<
    //   " " <<
    //   p_ace <<
    //   " " <<
    //   dc <<
    //   " " <<
    //   d_ace <<
    //   " " <<
    //   a <<
    //   " " <<
    //   r <<
    //   endl;
  }

  return r;
}


Real evaluation(Natural n_eval, mat &pi, vec &card_dist)
{
  Real E;

  Natural R = 0;
  for (Natural i = 0; i < n_eval; ++i)
    R += episode(pi, card_dist);

  E = (double) R / (double) n_eval;
  
  return E;
}


int main()
{
  Natural n = 200;
  Natural na = 2;

  vec card_dist = generate_stochastic_matrix(1, 13, true).transpose();
  stoch_mat pi = generate_stochastic_matrix(n, na, true);

  Real E = evaluation(1000, pi, card_dist);

  std::cout << E << std::endl;
  
  return 0;
}
