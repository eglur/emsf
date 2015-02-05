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


stoch_mat generate_stochastic_matrix(const Natural nrows, const Natural ncols, const bool constant = false)
{
  stoch_mat M;
  if (constant)
    M = mat::Constant(nrows,ncols,1);
  else
    M = mat::Random(nrows,ncols).cwiseAbs();
  normalize(M);

  return M;
}


v_stoch_mat generate_stochastic_matrices(const Natural nrows,const Natural ncols, const Natural na, const bool constant = false)
{
  v_stoch_mat M;
  M.resize(na);
  for(Natural a = 0; a < na; ++a)
    if (constant)
      M[a] = generate_stochastic_matrix(nrows,ncols,true);
    else
      M[a] = generate_stochastic_matrix(nrows,ncols);

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
  pt_vecn pt = sample(0, dist.size() - 1, 1);
  return(*pt)[0];
}


model generate_model(const Natural n, const Natural sr, const Natural na)
{
  v_stoch_mat D = generate_stochastic_matrices(n, sr, na);
  v_stoch_mat K = generate_stochastic_matrices(sr, n, na);

  model md;
  md.P.resize(na);
  for (Natural a = 0; a < na; ++a)
    md.P[a] = D[a] * K[a];

  md.mu = generate_stochastic_matrix(1, na).transpose();
  md.pi = generate_stochastic_matrix(n, na);

  return md;
}


data generate_data(model &md, const Natural T)
{
  try {
    data dt;
    dt.y.resize(T);
    dt.a.resize(T-1);

    // cout << "a" << endl; 
    dt.y[0] = sample_from_dist(md.mu);
    for (Natural i = 0; i < T-1; ++i) {
      // cout << "b" << endl; 
      dt.a[i] = sample_from_dist(md.pi.row(dt.y[i]));
      // cout << "c" << endl; 
      dt.y[i+1] = sample_from_dist(md.P[dt.a[i]].row(dt.y[i]));
    }

    return dt;
  }
  catch (int e) {
    cout << "yahoo" << endl;
  }
}


v_data generate_batch_data(model &md, const Natural T, const Natural num_batches)
{
  v_data v_dt;
  v_dt.resize(num_batches);

  for (Natural i = 0; i < num_batches; ++i)
    v_dt[i] = generate_data(md, T);

  return v_dt;
}


Real em_sf(model &md, v_data &dt, const Natural n, const Natural m, const Natural na, const Natural T, const Real eps = 1e-20, const Natural max_it = 10)
{
  v_stoch_mat P = md.P;
  vec mu = md.mu;
  stoch_mat pi = md.pi;

  v_stoch_mat D = generate_stochastic_matrices(n, m, na, true);
  v_stoch_mat K = generate_stochastic_matrices(m, n, na, true);

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

        row_vec k_row = K[a(t)].col(y(t+1)).transpose();
        B.row(t) = bb * pi(y(t+1), a(t+1)) * k_row;
        B.row(t) = B.row(t) / NF(t);
      }

      mat C = A.cwiseProduct(B);
      normalize(C);

      for (Natural t = 0; t < T-1; ++t) {
        D2[a(t)].row(y(t)) = D2[a(t)].row(y(t)).array() + C.row(t).array();
        K2[a(t)].col(y(t+1)) = K2[a(t)].col(y(t+1)).array() + C.row(t).transpose().array();
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
      ++P[a(t)](y(t), y(t+1));
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

  Natural nargs = 3;
  if (argc != nargs) {
    cout << "Usage: emsf run n" << endl;
    exit(EXIT_FAILURE);
  }

  const Natural run = atoi(argv[1]);
  const Natural n = atoi(argv[2]);
  const Natural sr = (Natural) (0.1 * n);
  const Natural na = 1;
  const Natural T = 10 * n * n;
  const Natural num_batches = 10;
  const Real eps = 1e20;
  const Natural max_it = 30;

  srand(run);

  model md = generate_model(n, sr, na);
  v_data dt = generate_batch_data(md, T, num_batches);

  Natural inc = ((10 * n * n) - n) / 9;
  for (Natural q = n; q <= T; q += inc) {
    Real e_cnt, e_emsf_a, e_emsf_b, e_emsf_c;
    clock_t begin, end;
    double t_cnt, t_emsf_a, t_emsf_b, t_emsf_c;

    begin = clock();
    e_cnt = counting(dt, num_batches, q, n, na, md.P);
    end = clock();
    t_cnt = double(end - begin) / CLOCKS_PER_SEC;
    
    begin = clock();
    e_emsf_a = em_sf(md, dt, n, (Natural) 0.8 * sr, na, q, eps, max_it);
    end = clock();
    t_emsf_a = double(end - begin) / CLOCKS_PER_SEC;
    
    begin = clock();
    e_emsf_b = em_sf(md, dt, n, sr, na, q, eps, max_it);
    end = clock();
    t_emsf_b = double(end - begin) / CLOCKS_PER_SEC;
    
    begin = clock();
    e_emsf_c = em_sf(md, dt, n, 1.25 * sr, na, q, eps, max_it);
    end = clock();
    t_emsf_c = double(end - begin) / CLOCKS_PER_SEC;
    

    // Log error
    filename.str(std::string());
    filename << "e_cnt_" << std::setw(2) << std::setfill('0') << run << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << e_cnt << " ";
    file.close();

    filename.str(std::string());
    filename << "e_emsf_a_" << std::setw(2) << std::setfill('0') << run << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << e_emsf_a << " ";
    file.close();

    filename.str(std::string());
    filename << "e_emsf_b_" << std::setw(2) << std::setfill('0') << run << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << e_emsf_b << " ";
    file.close();

    filename.str(std::string());
    filename << "e_emsf_c_" << std::setw(2) << std::setfill('0') << run << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << e_emsf_c << " ";
    file.close();


    // Log time
    filename.str(std::string());
    filename << "t_cnt_" << std::setw(2) << std::setfill('0') << run << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << t_cnt << " ";
    file.close();

    filename.str(std::string());
    filename << "t_emsf_a_" << std::setw(2) << std::setfill('0') << run << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << t_emsf_a << " ";
    file.close();

    filename.str(std::string());
    filename << "t_emsf_b_" << std::setw(2) << std::setfill('0') << run << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << t_emsf_b << " ";
    file.close();

    filename.str(std::string());
    filename << "t_emsf_c_" << std::setw(2) << std::setfill('0') << run << ".log";
    file.open(filename.str().c_str(), ios::app);
    file << t_emsf_c << " ";
    file.close();
  }

  return 0;
}
