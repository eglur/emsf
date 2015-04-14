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
#include "pisf.h"


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
    while (dc < 17)
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


Natural episode(mat &pi, vec &card_dist, vec &mu, const bool save, std::vector<Natural> &yv, std::vector<Natural> &av, std::vector<Natural> &rv)
{
  Natural pc = 0, p_ace = 0;
  Natural dc = 0, d_ace = 0;
  Natural s, a, r, sf;
  
  draw_card(pc, p_ace, card_dist);
  while (pc < 12)
    draw_card(pc, p_ace, card_dist);

  draw_card(dc, d_ace, card_dist);

  mu(get_s(pc, p_ace, dc)) += 1;

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
  vec mu(203);
  
  Natural R = 0;
  for (Natural i = 0; i < n_eval; ++i)
    R += episode(pi, card_dist, mu, false, yv, av, rv);

  E = (double) R / (double) n_eval;
  
  return E;
}


data_bj generate_data_bj(model &md, vec &card_dist)
{
  data_bj dt;
  episode(md.pi, card_dist, md.mu, true, dt.y, dt.a, dt.r);

  return dt;
}


v_data_bj generate_batch_data_bj(model &md, vec &card_dist, const Natural num_batches)
{
  v_data_bj v_dt;
  v_dt.resize(num_batches);

  for (Natural i = 0; i < num_batches; ++i)
    v_dt[i] = generate_data_bj(md, card_dist);

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


int main(int argc, char* argv[])
{
  ofstream file;
  stringstream filename;                                           // Nome do arquivo com os resultados
  stringstream id;                                                 // ID discriminante de cada experimento (contera os parametros utilizados)

  // Testa se os parâmetros foram informados corretamente
  Natural nargs = 10;
  if (argc != nargs) {
    cout << "Usage: blackjack RUN NUM_BATCHES NUM_EPISODES MIN_BATCHES NUM_POINTS MAX_IT GAMMA MAX_IT_PISF M" << endl;
    exit(EXIT_FAILURE);
  }

  // Parâmetros do experimento
  const Natural n = 203;
  const Natural sr = 20;
  const Natural na = 2;
  const Natural run = atoi(argv[1]);
  const Natural num_batches = atoi(argv[2]);
  const Natural num_episodes = atoi(argv[3]);
  const Natural min_batches = atoi(argv[4]);
  const Natural num_points = atoi(argv[5]);
  const Natural inc_batches = (double) (num_batches - min_batches) / (double) num_points + 1;
  const Natural max_it = atoi(argv[6]);
  const Real gamma_pisf = atof(argv[7]);
  const Real max_it_pisf = atof(argv[8]);
  const Natural m = atoi(argv[9]);

  // "Aleatoriza" com base no número do experimento (o número da
  // repetição do experimento que está sendo executada)
  srand(run);

  // Distribuição de probabilidade para cada carta do baralho
  // (Blackjack)
  vec card_dist = generate_stochastic_matrix(1, 13, true).transpose();

  // D e K utilizadas em todos os experimentos (D e K são modificadas
  // em em_sf, portanto é necessário guardá-las para realizar um novo,
  // e.g., com diferentes valores de batches)
  v_stoch_mat D_original = generate_stochastic_matrices(n, m, na);
  v_stoch_mat K_original = generate_stochastic_matrices(m, n, 1);

  model md;
  md.mu = vec::Zero(n, 1);
  md.pi = generate_stochastic_matrix(n, na, true);

  v_data_bj dt = generate_batch_data_bj(md, card_dist, num_batches);

  // Normalize md.mu
  mat mu2 = md.mu.transpose();
  normalize(mu2);
  md.mu = mu2.transpose();

  for (Natural nb = min_batches; nb <= num_batches; nb += inc_batches) {
    clock_t begin, end;                                            // Para cronometrar experimento
    v_stoch_mat D = D_original;
    v_stoch_mat K = K_original;

    begin = clock();                                               // Inicia cronômetro
    em_sf_sk(md, dt, n, m, na, nb, D, K, max_it);                  // Executa em_sf_sk

    vec r_hat = vec::Zero(n, 1);                                   // Obtém dados para o PISF
    r_hat(200) = -1.0;
    r_hat(201) = 0.0;
    r_hat(202) = 1.0;
  
    pt_agent agt = pisf(D, K[0], K[0] * r_hat, gamma_pisf, max_it_pisf); // Executa PISF
    end = clock();                                                 // Pára o cronômetro
    double t_emsf_bj = double(end - begin) / CLOCKS_PER_SEC;       // Calcula tempo gasto em segundos

    // Obtém política estocástica (PISF retorna uma política
    // determinística)
    stoch_mat pi_stc = stoch_mat::Zero(n, na);
    for (Natural i = 0; i < n; ++i)
      pi_stc(i, agt->pi(i)) = 1;

    // Avalia a política no domínio
    Real v_emsf_bj = evaluation(num_episodes, pi_stc, card_dist);

    // Obtém o ID para o arquivo com os resultados do experimento
    id.str(std::string());
    id << num_batches << "_"
       << num_episodes << "_"
       << min_batches << "_"
       << num_points << "_"
       << num_points << "_"
       << max_it << "_"
       << gamma_pisf << "_"
       << max_it_pisf << "_"
       << m << "_"
       << std::setw(2) << std::setfill('0') << run;

    // Log time (persiste tempo gasto)
    filename.str(std::string());
    filename << "t_emsf_bj_" << id.str() << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << t_emsf_bj << " ";
    file.close();

    // Log value (persiste valor da política obtida com o EM-SF-SK)
    filename.str(std::string());
    filename << "v_emsf_bj_" << id.str() << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << v_emsf_bj << " ";
    file.close();
  }

  return 0;
}
