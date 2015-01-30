#ifndef POLICY_ITERATION_H
#define POLICY_ITERATION_H

#include "mdp.h"
#include "agent.h"

namespace dp
{

// tolerance for computing the value function
const Real const_epsilon_pi = 1e-6;


/* Auxiliary functions */
template <class MDP_TYPE>
inline void compute_policy(const MDP_TYPE &M, const vec &v, policy &pi, const Real gamma, mat &Q)
{
    for (Natural a = 0; a < M.num_actions(); ++a ) Q.col(a) = M.r(a) + gamma * M.P(a) * v;
    for (Natural i = 0; i < M.num_states(); ++i) Q.row(i).maxCoeff(&pi[i]);
}


template <class MAT_TYPE>
void compute_v_iteratively(const MAT_TYPE &Ppi, const vec &rpi, vec &v, const Real gamma, const Real epsilon)
{
    /* This stop criterion is based on Puterman, Proposition 6.6.5, Eq. (6.6.11). It guarantees that
     ||v_true - v_approximate|| < epsilon */

    Real tolerance = epsilon * (1 - gamma) / gamma;
//    Natural it = 0;
    v = rpi;
    vec v_tmp(v.size());
    do
    {
        v_tmp = v;
        v = rpi + gamma * Ppi * v;
        v_tmp = v - v_tmp;
    }
    while (v_tmp.maxCoeff() - v_tmp.minCoeff() > tolerance); // && ++it < rpi.size()); /// review the upper bound for "it"
}


// Dense MDPs
inline mat *generate_Markov_process(const mdp &M, const policy &pi, vec &rpi)
{
    mat *Ppi = new mat(M.num_states(), M.num_states());
    for (Natural i = 0; i < M.num_states(); ++i) // Não precisava da matriz Ppi...
    {
       Ppi->row(i) = M.P(pi[i]).row(i);
       rpi[i] = M.r(pi[i],i);
    }
    return Ppi;
}



// Sparse MDPs
inline smat *generate_Markov_process(const smdp &M, const policy &pi, vec &rpi)
{

    smat *Ppi = new smat(M.num_states(), M.num_states());
    for (Natural i = 0; i < M.num_states(); ++i)
    {
       for (smat::InnerIterator it(M.P(pi[i]),i); it; ++it) Ppi->insert(i, it.col() ) = it.value();
       rpi[i] = M.r(pi[i])[i];
    }

    Ppi->finalize();
    return Ppi;
}

/* Main function */
template <class MDP_TYPE>
pt_agent policy_iteration(const MDP_TYPE &M, const Real gamma,
                          vecn *pi = NULL, const Real epsilon = const_epsilon_pi, const int max_it = -1)
{
    pt_agent agt(new agent(M.num_states(), M.num_actions()));

    bool destroy = false;
    if (pi == NULL)
    {
        pi = new vecn(M.num_states());
        for (int i = 0; i < pi->size(); ++i) (*pi)[i] = random_Natural(0, M.num_actions()-1);
        destroy = true;
    }

    policy pi_old(*pi);
    Natural max_iter = (max_it > 0) ? max_it : M.num_states(); /// think
    Natural it = 0;
    vec v;
    do
    {

        vec rpi(M.num_states());
        typename MDP_TYPE::mat_type *Ppi = generate_Markov_process(M, *pi, rpi); //using pointer because the assignment operator (=) seems to be buggy

        compute_v_iteratively(*Ppi, rpi, v, gamma, epsilon);
        delete Ppi;

        // Policy update
        pi_old = *pi;
        compute_policy(M, v, *pi, gamma, agt->Q());

    }
    while (*pi != pi_old && ++it < max_iter);

    if (destroy) delete pi;

    return agt;
}

/* currently not used:
inline void compute_v_exactly(const mat &Ppi, const vec &rpi, vec &v, const Real gamma, const Real epsilon)
{
    v = (mat::Identity(Ppi.rows(), Ppi.cols()) - gamma * Ppi).lu().solve(rpi);
}
*/

}

#endif // POLICY_ITERATION_H
