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


inline Real frobenius_norm(mat &A, mat &B)
{
  return (A.array() - B.array()).array().square().sum();
}


inline Real frobenius_norm_v(v_mat &A, v_mat &B, const Natural size)
{
  Real e = 0.0;
  for (Natural i = 0; i < size; ++i)
    e += frobenius_norm(A[i], B[i]);

  return e;
}


inline std::string date_time_str()
{
  time_t rawtime;
  struct tm *timeinfo;
  char buffer[80];

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(buffer, 80, "%F_%T", timeinfo);
  std::string str(buffer);

  return str;
}


inline void normalize(mat &A)
{
  for (Natural i = 0; i < A.rows(); ++i) {
    Real sum = A.row(i).sum();
    if (sum > 0) for (Natural j = 0; j < A.cols(); ++j) A(i,j) /= sum;
  }
}


stoch_mat generate_stochastic_matrix(const Natural nrows, const Natural ncols, const bool constant = false, const bool one = false)
{
  stoch_mat M;
  if (one) {
    M = mat::Zero(nrows, ncols);
    for (Natural i = 0; i < nrows; ++i) {
      M(i, util::random_Natural(0, ncols-1)) = 1.0;
    }
  }
  else {
    if (constant)
      M = mat::Constant(nrows,ncols,1);
    else
      M = mat::Random(nrows,ncols).cwiseAbs();
    normalize(M);
  }

  return M;
}


v_stoch_mat generate_stochastic_matrices(const Natural nrows,const Natural ncols, const Natural na, const bool constant = false, const bool one = false)
{
  v_stoch_mat M;
  M.resize(na);
  for(Natural a = 0; a < na; ++a)
    if (one) {
      M[a] = generate_stochastic_matrix(nrows, ncols, false, true);
    }
    else {
      if (constant)
        M[a] = generate_stochastic_matrix(nrows, ncols, true);
      else
        M[a] = generate_stochastic_matrix(nrows, ncols, false);
    }

  return M;
}


v_mat generate_zero_matrices(const Natural nrows, const Natural ncols, const Natural na)
{
  v_mat vm;
  vm.resize(na);
  for (Natural a = 0; a < na; ++a)
    vm[a] = mat::Zero(nrows, ncols);

  return vm;
}


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


model generate_model(const Natural n, const Natural sr, const Natural na)
{
  v_stoch_mat D = generate_stochastic_matrices(n, sr, na);
  v_stoch_mat K = generate_stochastic_matrices(sr, n, na);

  model md;
  md.P.resize(na);
  for (Natural a = 0; a < na; ++a)
    md.P[a] = D[a] * K[a];

  md.mu = generate_stochastic_matrix(1, n).transpose();
  md.pi = generate_stochastic_matrix(n, na);

  return md;
}


data generate_data(model &md, const Natural T)
{
  data dt;
  dt.y.resize(T);
  dt.a.resize(T-1);

  dt.y[0] = sample_from_dist(md.mu);
  for (Natural i = 0; i < T-1; ++i) {
    dt.a[i] = sample_from_dist(md.pi.row(dt.y[i]));
    dt.y[i+1] = sample_from_dist(md.P[dt.a[i]].row(dt.y[i]));
  }

  return dt;
}


v_data generate_batch_data(model &md, const Natural T, const Natural num_batches)
{
  v_data v_dt;
  v_dt.resize(num_batches);

  for (Natural i = 0; i < num_batches; ++i)
    v_dt[i] = generate_data(md, T);

  return v_dt;
}


Real em_sf(model &md, v_data &dt, const Natural n, const Natural m, const Natural na, const Natural T, v_stoch_mat D, v_stoch_mat K, const Real eps = 1e-20, const Natural max_it = 10)
{
  v_stoch_mat P = md.P;
  vec mu = md.mu;
  stoch_mat pi = md.pi;

  Real old_score = minf;
  Real score = 0.0;
  Natural it = 0;

  while (abs(score - old_score) > eps || it < max_it) {
    old_score = score;
    score = 0.0;

    v_mat D2 = generate_zero_matrices(n, m, na);
    v_mat K2 = generate_zero_matrices(m, n, na);

    for (Natural batch = 0; batch < dt.size(); ++batch) {
      vecn y = dt[batch].y;
      vecn a = dt[batch].a;
      
      mat A = mat::Zero(T-1, m);
      mat B = mat::Zero(T-1, m);
      vec NF = vec::Zero(T-1);

      A.row(0) = mu[y[0]] * pi(y[0], a[0]) * D[a[0]].row(y[0]);
      NF[0] = A.row(0).sum();
      A.row(0) = A.row(0) / NF[0];

      for (Natural t = 1; t < T-1; ++t) {
        row_vec k_row = K[a[t-1]].col(y[t]).transpose();
        Real aa = A.row(t-1).cwiseProduct(k_row).sum();

        A.row(t) = aa * pi(y[t], a[t]) * D[a[t]].row(y[t]);
        NF[t] = A.row(t).sum();
        A.row(t) = A.row(t) / NF[t];
      }

      Natural B_T = T-2;                                        // Different index strategy on the R language
      B.row(B_T) = K[a[B_T]].col(y[B_T]).transpose();
      B.row(B_T) = B.row(B_T) / NF[B_T];
      for (Natural t = B_T-1; t >= 0; --t) {
        row_vec d_row = D[a[t+1]].row(y[t+1]);
        Real bb = B.row(t+1).cwiseProduct(d_row).sum();

        row_vec k_row = K[a[t]].col(y[t+1]).transpose();
        B.row(t) = bb * pi(y[t+1], a[t+1]) * k_row;
        B.row(t) = B.row(t) / NF[t];
      }

      mat C = A.cwiseProduct(B);
      normalize(C);

      for (Natural t = 0; t < T-1; ++t) {
        D2[a[t]].row(y[t]) = D2[a[t]].row(y[t]).array() + C.row(t).array();
        K2[a[t]].col(y[t+1]) = K2[a[t]].col(y[t+1]).array() + C.row(t).transpose().array();
      }

      score = score - NF.array().inverse().log().sum();
    }

    for (Natural i = 0; i < na; ++i) {
      normalize(D2[i]);
      normalize(K2[i]);
    }

    D = D2;
    K = K2;

    ++it;
  }

  v_stoch_mat P_bar(na);
  for (Natural i = 0; i < na; ++i)
    P_bar[i] = D[i] * K[i];

  return frobenius_norm_v(P_bar, P, na);
}


v_mat get_P_by_counting(v_data &dt, const Natural num_batches, const Natural T, const Natural n, const Natural na)
{
  v_mat P = generate_zero_matrices(n, n, na);

  for (Natural batch = 0; batch < num_batches; ++batch) {
    vecn y = dt[batch].y, a = dt[batch].a;
    for (Natural t = 0; t < T-1; ++t)
      P[a[t]](y[t], y[t+1]) = P[a[t]](y[t], y[t+1]) + 1;
  }

  for (Natural i = 0; i < na; ++i)
    normalize(P[i]);

  return P;
}


Real counting(v_data &dt, const Natural num_batches, const Natural T, const Natural n, const Natural na, v_mat &P)
{
  v_mat P_cnt = get_P_by_counting(dt, num_batches, T, n, na);
  return frobenius_norm_v(P_cnt, P, na);
}


int main(int argc, char* argv[])
{
  ofstream file;
  stringstream filename;
  stringstream id;

  Natural nargs = 3;
  if (argc != nargs) {
    cout << "Usage: emsf run n" << endl;
    exit(EXIT_FAILURE);
  }

  const Natural run = atoi(argv[1]);
  const Natural n = atoi(argv[2]);
  const Natural na = 1;
  const Natural T = (Natural) 10 * (n * n);
  const Natural num_batches = 10;
  const Real eps = 1e20;

  Natural max_it_min = 1;
  Natural max_it_max = 10;
  Natural max_it_inc = 1;

  srand(run);

  // sr values: factors, increment, quantity and the own vector
  Real srf_min = 0.2;
  Real srf_max = 0.2;
  Real srf_inc = 0.1;

  // Just can't see why, but only works this way
  Real srf_qty_real = ((srf_max - srf_min) / srf_inc) + 1.0;
  Natural srf_qty = (Natural) srf_qty_real;
  
  std::vector<Natural> sr;
  for (Real srf = srf_min; srf <= srf_max; srf += srf_inc)
    sr.push_back((Natural) (srf * (Real) n));

  std::vector<model> md(srf_qty);
  for (Natural i = 0; i < srf_qty; ++i)
    md[i] = generate_model(n, sr[i], na);

  std::vector<v_data> dt(srf_qty);
  for (Natural i = 0; i < srf_qty; ++i)
    dt[i] = generate_batch_data(md[i], T, num_batches);

  // m values: factors, increment, quantity and the own vector
  Real mf_min = 0.25;
  Real mf_max = 1.25;
  Real mf_inc = 0.25;

  Natural q_min = n;
  Natural q_max = T;
  Natural n_points = 100;
  Natural q_inc = (q_max - q_min) / (n_points - 1);
  for (Natural i = 0; i < srf_qty; ++i) {
    clock_t begin, end;

    // EM-SF
    for (Real mf = mf_min; mf <= mf_max; mf += mf_inc) {
      Real m_real = mf * (Real) sr[i];
      Natural m = (Natural) m_real;

      v_stoch_mat D = generate_stochastic_matrices(n, m, na);
      v_stoch_mat K = generate_stochastic_matrices(m, n, na);
      for (Natural max_it = 1; max_it <= max_it_max; max_it += max_it_inc) {
        for (Natural q = q_min; q <= q_max; q += q_inc) {
          Real e_emsf;
          double t_emsf;

          // Calculate
          begin = clock();
          e_emsf = em_sf(md[i], dt[i], n, m, na, q, D, K, eps, max_it);
          end = clock();
          t_emsf = double(end - begin) / CLOCKS_PER_SEC;

          // Log error
          id.str(std::string());
          id << n << "_"
             << na << "_"
             << T << "_"
             << num_batches << "_"
             << eps << "_"
             << max_it << "_"
             << std::setw(4) << std::setfill('0') << sr[i] << "_"
             << std::setw(4) << std::setfill('0') << m << "_"
             << std::setw(2) << std::setfill('0') << run;

          filename.str(std::string());
          filename << "e_emsf_" << id.str() << ".log";
          file.open(filename.str().c_str(), ios::app);
          file << e_emsf << " ";
          file.close();

          // Log time
          filename.str(std::string());
          filename << "t_emsf_" << id.str() << ".log";
          file.open(filename.str().c_str(), ios::app);
          file << t_emsf << " ";
          file.close();
        }
      }
    }

    // Counting
    for (Natural q = q_min; q <= q_max; q += q_inc) {
      Real e_cnt;
      double t_cnt;

      begin = clock();
      e_cnt = counting(dt[i], num_batches, q, n, na, md[i].P);
      end = clock();
      t_cnt = double(end - begin) / CLOCKS_PER_SEC;

      // Log error
      id.str(std::string());
      id << n << "_"
         << na << "_"
         << T << "_"
         << num_batches << "_"
         << eps << "_"
         << std::setw(4) << std::setfill('0') << sr[i] << "_"
         << std::setw(2) << std::setfill('0') << run;

      filename.str(std::string());
      filename << "e_cnt_" << id.str() << ".log";
      file.open(filename.str().c_str(), ios::app);
      file << e_cnt << " ";
      file.close();

      // Log time
      filename.str(std::string());
      filename << "t_cnt_" << id.str() << ".log";
      file.open(filename.str().c_str(), ios::app);
      file << t_cnt << " ";
      file.close();
    }
  }

  return 0;
}
