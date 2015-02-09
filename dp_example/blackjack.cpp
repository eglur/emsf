#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
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


Natural transition(Natural &pc, Natural &p_ace, Natural &dc, Natural &d_ace, Natural &r, Natural a, vec &card_dist, Natural &sf)
{
  Natural stick = 0;
  Natural hit = 1;
  
  if (a == hit) {
    draw_card(pc, p_ace, card_dist);

    if (pc > 21) {
      r = -1;
      sf = 200;
    }
    else
      r = 0;
  }
  else if (a == stick) {
    while (dc < 17)                                             /* TODO: checar se é < ou <= */
      draw_card(dc, d_ace, card_dist);

    if (dc > 21) {
      r = 1;
      sf = 202;
    }
    else {
      Natural p_diff = 21 - pc;
      Natural d_diff = 21 - dc;

      if (p_diff < d_diff) {
        r = 1;
        sf = 202;
      }
      else if (p_diff > d_diff) {
        r = -1;
        sf = 200;
      }
      else {
        r = 0;
        sf = 201;
      }
    }
  }
}


Natural get_s(Natural pc, Natural p_ace, Natural dc)
{
  return (p_ace * 100) + ((pc - 12) * 10) + (dc - 2);
}


Natural get_f(Natural s, Natural &pc, Natural &p_ace, Natural &dc)
{
  p_ace = s / 100;

  s -= p_ace * 100;
  pc = s / 10;

  s -= pc * 10;
  dc = s;
}


inline Natural get_a(Natural s, mat &pi)
{
  return sample_from_dist(pi.row(s).transpose());
}


Natural episode(mat &pi, vec &card_dist, const bool save, std::vector<Natural> &yv, std::vector<Natural> &av, std::vector<Natural> &rv)
{
  Natural pc = 0, p_ace = 0;
  Natural dc = 0, d_ace = 0;
  Natural s, a, r, sf;
  
  draw_card(pc, p_ace, card_dist);
  while (pc < 12)
    draw_card(pc, p_ace, card_dist);

  draw_card(dc, d_ace, card_dist);

  sf = 0;
  while (sf < 200) {
    s = get_s(pc, p_ace, dc);
    if (save) yv.push_back(s);

    a = get_a(s, pi);
    if (save) av.push_back(a);

    transition(pc, p_ace, dc, d_ace, r, a, card_dist, sf);

    if (save) rv.push_back(r);
  }

  if (save) yv.push_back(sf);

  return r;
}


Real evaluation(Natural n_eval, mat &pi, vec &card_dist)
{
  Real E;
  std::vector<Natural> yv, av, rv;

  Natural R = 0;
  for (Natural i = 0; i < n_eval; ++i)
    R += episode(pi, card_dist, false, yv, av, rv);

  E = (double) R / (double) n_eval;
  
  return E;
}


data generate_data_bj(model &md, vec &card_dist, const bool save)
{
  data dt;
  episode(md.pi, card_dist, true, dt.y, dt.a, dt.r);

  return dt;
}


v_data generate_batch_data_bj(model &md, const Natural T, const Natural num_batches)
{
  v_data v_dt;
  v_dt.resize(num_batches);

  for (Natural i = 0; i < num_batches; ++i)
    v_dt[i] = generate_data(md, T);

  return v_dt;
}


int main(int argc, char* argv[])
{
  Natural nargs = 3;
  if (argc != nargs) {
    cout << "Usage: blackjack eval_qty eval_size" << endl;
    exit(EXIT_FAILURE);
  }

  const Natural n = 200;
  const Natural na = 2;
  const Natural num_batches = atoi(argv[1]);
  const Natural num_episodes = atoi(argv[2]);

  srand(time(NULL));

  vec card_dist = generate_stochastic_matrix(1, 13, true).transpose();
  stoch_mat pi = generate_stochastic_matrix(n, na, true);

  for (Natural i = 0; i < num_batches; ++i) {
    Natural R = 0;
    for (Natural j = 0; j < num_episodes; ++j) {
      std::vector<Natural> yv, av, rv;
      R += episode(pi, card_dist, true, yv, av, rv);

      Natural k = 0, pc, p_ace, dc;
      while (k < yv.size() - 1) {
        get_f(yv[k], pc, p_ace, dc);
        cout <<"s= ["
             << std::setw(1) << p_ace << ", "
             << std::setw(2) << pc << ", "
             << std::setw(2) << dc << "] ("
             << std::setw(3) << yv[k] << "), a="
             << std::setw(1) << av[k] << ", r="
             << std::setw(2) << rv[k] << endl;;
        ++k;
      }
      cout << std::setw(3) << yv[k] << "." << endl << endl;
    }

    Real E = (double) R / (double) num_episodes;
    cout << "E: " << E << endl << endl;
  }
  cout << "================================================================================" << endl << endl;

  return 0;
}
