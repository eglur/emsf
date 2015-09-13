#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include "blackjack.h"
#include "emsf.h"
#include "policy_iteration.h"


using namespace std;
using namespace dp;
using namespace std;
using namespace Eigen;
using namespace blackjack;
using namespace util;
using namespace emsf;
using namespace dp;


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


Natural episode(mat &pi, Real epsilon, vec &card_dist, const bool save, std::vector<Natural> &yv, std::vector<Natural> &av, std::vector<Natural> &rv)
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

    if (((double) rand() / (RAND_MAX)) <= epsilon)
      a = rand() % 2; // stick == 0; hit == 1
    else
      a = get_a(s, pi);
    
    if (save) av.push_back(a);

    transition(pc, p_ace, dc, d_ace, r, a, card_dist, sf);

    if (save) rv.push_back(r);
  }

  if (save) yv.push_back(sf);

  return r;
}


v_mat get_P_by_counting_bj(v_data_bj &dt, const Natural num_batches, const Natural n, const Natural na)
{
  v_mat P = generate_zero_matrices(n, n, na);

  for (Natural batch = 0; batch < num_batches; ++batch) {
    std::vector<Natural> y = dt[batch].y;
    std::vector<Natural> a = dt[batch].a;

    Natural T = y.size();
    for (Natural t = 0; t < T-1; ++t)
      P[a[t]](y[t], y[t+1]) = P[a[t]](y[t], y[t+1]) + 1;
  }

  for (Natural i = 0; i < na; ++i) {
    normalize(P[i]);
    P[i](200, 200) = 1.0;
    P[i](201, 201) = 1.0;
    P[i](202, 202) = 1.0;
  }

  return P;
}


v_mat get_R_by_counting_bj(v_data_bj &dt, const Natural num_batches, const Natural n, const Natural na)
{
  v_mat r_count = generate_zero_matrices(n, 1, na);
  v_mat r_sum = generate_zero_matrices(n, 1, na);
  v_mat r_mean = generate_zero_matrices(n, 1, na);

  for (Natural batch = 0; batch < num_batches; ++batch) {
    std::vector<Natural> y = dt[batch].y;
    std::vector<Natural> a = dt[batch].a;
    std::vector<Natural> r = dt[batch].r;

    Natural T = y.size();
    for (Natural t = 0; t < T-1; ++t) {
      r_count[a[t]](y[t]) += 1;
      r_sum[a[t]](y[t]) += r[t];
    }
  }

  for (Natural a = 0; a < na; ++a) {
    for (Natural s = 0; s < n; ++s) {
      if (r_count[a](s) > 0.0) {
        r_mean[a](s) = r_sum[a](s) / r_count[a](s);
      }
    }
  }

  return r_mean;
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


data_bj generate_data_bj(mat &pi, vec &card_dist)
{
  data_bj dt;
  episode(pi, card_dist, true, dt.y, dt.a, dt.r);

  return dt;
}


data_bj generate_data_bj(mat &pi, Real epsilon, vec &card_dist)
{
  data_bj dt;
  episode(pi, epsilon, card_dist, true, dt.y, dt.a, dt.r);

  return dt;
}


v_data_bj generate_batch_data_bj(model &md, vec &card_dist, const Natural num_batches)
{
  v_data_bj v_dt;
  v_dt.resize(num_batches);

  for (Natural i = 0; i < num_batches; ++i)
    v_dt[i] = generate_data_bj(md.pi, card_dist);

  return v_dt;
}


void print_data_bj(data_bj dt)
{
  Natural k = 0, pc, p_ace, dc;

  while (k < dt.y.size() - 1) {
    get_f(dt.y[k], pc, p_ace, dc);
    cout <<"s= ["
         << std::setw(1) << p_ace << ", "
         << std::setw(2) << pc << ", "
         << std::setw(2) << dc << "] ("
         << std::setw(3) << dt.y[k] << "), a="
         << std::setw(1) << dt.a[k] << ", r="
         << std::setw(2) << dt.r[k] << endl;;
    ++k;
  }
  cout << std::setw(3) << dt.y[k] << "." << endl << endl;
}


void print_batch_data_bj(v_data_bj dt, const Natural num_batches)
{
  for (Natural i = 0; i < num_batches; ++i)
    print_data_bj(dt[i]);
}


void count_transitions(v_mat &C, data_bj &dt)
{
  for (Natural t = 0; t < dt.y.size() - 1; ++t)
    C[dt.a[t]](dt.y[t], dt.y[t+1]) += 1;
}


void get_P_from_C(v_mat &C, Natural n, Natural na) {
  v_mat P = generate_zero_matrices(n, n, na);
  
  for (Natural a = 0; a < na; ++a) {
    P[a] = C[a];
    normalize(P[a]);
  }
}


int main(int argc, char* argv[])
{
  ofstream file;
  stringstream filename;
  stringstream id;

  Natural nargs = 6;
  if (argc != nargs) {
    cout << "Usage: blackjack run num_batches num_episodes min_batches num_points" << endl;
    exit(EXIT_FAILURE);
  }

  const Natural n = 203;
  const Natural sr = 20;
  const Natural na = 2;
  const Natural run = atoi(argv[1]);
  const Natural num_batches = atoi(argv[2]);
  const Natural num_episodes = atoi(argv[3]);
  const Natural min_batches = atoi(argv[4]);
  const Natural num_points = atoi(argv[5]);
  const Natural inc_batches = (double) (num_batches - min_batches) / (double) num_points + 1;

  srand(run);

  vec card_dist = generate_stochastic_matrix(1, 13, true).transpose();
  stoch_mat pi = generate_stochastic_matrix(n, na, true);

  model md = generate_model(n, sr, na);
  data_bj dt = generate_data_bj(pi, card_dist);
  v_mat C = generate_zero_matrices(n, n, na);
  
  for (Natural nb = min_batches; nb <= num_batches; nb += inc_batches) {
    clock_t begin, end;
    Real v_cnt;
    double t_cnt;

    begin = clock();
    v_mat P = get_P_from_C(C, n, na);
    v_mat r = get_R_by_counting_bj(dt, nb, n, na);

    mdp M(n, na);
    for (Natural a = 0; a < na; ++a) {
      M.P(a) = P[a];
      M.r(a) = r[a];
    }

    // Solve the MDP using policy iteration
    Real gamma = 1.0;
    pt_agent agt = policy_iteration(M, gamma);
    end = clock();
    t_cnt = double(end - begin) / CLOCKS_PER_SEC;

    vecn pi_det = *agt->pi();
    stoch_mat pi_stc = mat::Zero(n, na);
    for (Natural s = 0; s < n; ++s)
      pi_stc(s, pi_det[s]) = 1.0;

    v_cnt = evaluation(num_episodes, pi_stc, card_dist);

    id.str(std::string());
    id << sr << "_"
       << na << "_"
       << num_batches << "_"
       << num_episodes << "_"
       << min_batches << "_"
       << num_points << "_"
       << std::setw(2) << std::setfill('0') << run;

    // Log time
    filename.str(std::string());
    filename << "t_cnt_bj_" << id.str() << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << t_cnt << " ";
    file.close();

    // Log value
    filename.str(std::string());
    filename << "v_cnt_bj_" << id.str() << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << v_cnt << " ";
    file.close();
  }

  return 0;
}
