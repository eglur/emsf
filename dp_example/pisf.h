#ifndef PISF_H
#define PISF_H

#include "policy_iteration.h"


namespace dp{


inline void compute_policy(const vector<mat> &Da, const vec &vb, policy &pi, const Real gamma, mat &Q)
{
    for (size_t a = 0; a < Da.size(); ++a) Q.col(a) = Da[a] * vb;
    for (Natural i = 0; i < pi.size(); ++i) Q.row(i).maxCoeff(&pi[i]);
}


pt_agent pisf(const vector<mat> &Da, const mat &K, const vec &rb, const Real gamma, const int max_it, const Real epsilon = const_epsilon_pi);
void     pisf(vec &vb, const vector<mat> &Da, const mat &K, const vec &rb, const Real gamma, const int max_it, const Real epsilon = const_epsilon_pi);


// Sparse version

inline void compute_policy(const vector<smat> &Da, const vec &vb, policy &pi, const Real gamma, mat &Q)
{
    for (size_t a = 0; a < Da.size(); ++a) Q.col(a) = Da[a] * vb;
    for (Natural i = 0; i < pi.size(); ++i) Q.row(i).maxCoeff(&pi[i]);
}

pt_agent pisf(const vector<smat> &Da, const smat &K, const vec &rb, const Real gamma, const int max_it, const Real epsilon = const_epsilon_pi);
void     pisf(vec &vb, const vector<smat> &Da, const smat &K, const vec &rb, const Real gamma, const int max_it, const Real epsilon = const_epsilon_pi);


}

#endif // PISF_H
