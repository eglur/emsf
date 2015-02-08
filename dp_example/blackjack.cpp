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

  if (ind > dist.size() - 1)
    cout << "ERRO MAIOR!!!";
  
  if (ind < 0)
    cout << "ERRO MENOR!!!";
  
  return ind;
}


Natural card(mat &D)
{
  Natural c = sample_from_dist(D) + 1;
  if (c > 10) c = 10;

  return c;
}


Natural draw_card(Natural &xc, Natural &x_ace, mat &D)
{
  Natural c = card(D);

  xc += c;

  if (c == 1 && !x_ace) {
    xc += 10;
    x_ace = 1;
  }

  if (xc > 21 && x_ace) {
    xc -= 10;
    x_ace = 0;
  }
}


bool transition(Natural &pc, Natural &p_ace, Natural &dc, Natural &d_ace, Natural &r, Natural a, bool &end, mat &D)
{
  Natural stick = 0;
  Natural hit = 1;
  
  if (a == hit) {
    draw_card(pc, p_ace, D);

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
      draw_card(dc, d_ace, D);

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

  return 0;
}


Natural episode(mat &pi, mat &D)
{
  Natural pc = 0;
  Natural p_ace = 0;
  Natural dc = 0;
  Natural d_ace = 0;
  Natural r;
  
  draw_card(pc, p_ace, D);
  draw_card(pc, p_ace, D);

  draw_card(dc, d_ace, D);

  bool end = false;
  while (!end) {
    Natural s = get_s(pc, p_ace, dc);
    Natural a = get_a(s, pi);
    transition(pc, p_ace, dc, d_ace, a, r, end, D);
  }

  return r;
}


Real evaluation(Natural n_eval, mat &pi, mat &D)
{
  Real E;

  Natural R = 0;
  for (Natural i = 0; i < n_eval; ++i)
    R += episode(pi, D);

  E = (double) R / (double) n_eval;
  
  return E;
}
