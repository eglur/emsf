#ifndef VALUE_ITERATION_H
#define VALUE_ITERATION_H

#include "mdp.h"
#include "agent.h"

namespace dp
{

    const Real const_epsilon_vi = 1e-6;


    template <class MDP_TYPE>
    pt_agent value_iteration(const MDP_TYPE &M, const Real gamma, const int max_it = -1, const Real epsilon = const_epsilon_vi)
    {
        vector<vec> Q(M.num_actions());
        for (Natural a = 0; a < M.num_actions(); ++a) Q[a] = M.r(a);
        value_iteration(Q, M, gamma, max_it, epsilon);
        pt_agent agt(new agent(M.num_states(), M.num_actions()));
        for (Natural a = 0; a < M.num_actions(); ++a) agt->Q().col(a) = Q[a];
        return agt;
    }


    inline void compute_v(vec &v, const vector<vec> &Q)
    {
        for (Natural i = 0; i < v.size(); ++i)
        {
            Real q_max = Q[0][i];
            for (size_t a = 1; a < Q.size(); ++a) if (Q[a][i] > q_max) q_max = Q[a][i];
            v[i] = q_max;
        }
    }

    template <class MDP_TYPE>
    void value_iteration(vector<vec> &Q, const MDP_TYPE &M, const Real gamma, const int max_it = -1, const Real epsilon = const_epsilon_vi)
    {
        const Natural max_it2 = (max_it > 0) ? max_it : M.num_states();  /// think
        /* This stop criterion is based on Puterman, Proposition 6.6.5, Eq. (6.6.11). It guarantees that
        ||v_true - v_approximate|| < epsilon */
        Real tolerance = epsilon * (1 - gamma) / (gamma);
        Natural it = 0;

        vec v(M.num_states());
        compute_v(v,Q);
        vec v_tmp;
        do
        {
            v_tmp = v;
            for (Natural a = 0; a < M.num_actions(); ++a) Q[a] = M.r(a) + gamma * M.P(a) * v;
            compute_v(v,Q);
            v_tmp = v - v_tmp;
            if (max_it > 0) ++it;
        } while (v_tmp.maxCoeff() - v_tmp.minCoeff() > tolerance && it < max_it2);
    }



}

#endif // VALUE_ITERATION_H
