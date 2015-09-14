#include "pisf.h"


namespace dp
{


    void pisf(vec &vb, const vector<mat> &Da, const mat &K, const vec &rb, const Real gamma, const int max_it, const Real epsilon)
    {
        Natural n  = K.cols();
        Natural m = K.rows();
        Natural num_actions = Da.size();

        mat Q(n, num_actions);

        policy pi(n);
        for (int i = 0; i < n; ++i) pi[i] = random_Natural(0, num_actions-1);

        policy pi_old(pi);
        int max_it2 = max_it;
        if (max_it2 < 0) max_it2 = n;
        Natural it = 0;
        do
        {
            mat Dpi = mat(int(n), int(m));

            // Policy evaluation
            for (int i = 0; i < Dpi.rows(); ++i) Dpi.row(i) = Da[pi[i]].row(i); // Não precisava da matrix Dpi...
            compute_v_iteratively(K * Dpi, rb, vb, gamma, epsilon);

            // Policy update
            pi_old = pi;

//            compute_policy(M, Dpi * vb, pi, gamma, (*agt)()); // Slower, but more precise (like in the paper)
            compute_policy(Da, vb, pi, gamma, Q);

            ++it;
        }
        while (pi != pi_old && it < max_it2);
    }

    pt_agent pisf(const vector<mat> &Da, const mat &K, const vec &rb, const Real gamma, const int max_it, const Real epsilon)
    {
        pt_agent agt(new agent(K.cols(), Da.size()));
        vec vb;
        pisf(vb, Da, K, rb, gamma, max_it, epsilon);
        for (size_t a = 0; a < Da.size(); ++a) agt->Q().col(a) = Da[a] * vb;
        return agt;
    }



// Sparse version
    void pisf(vec &vb, const vector<smat> &Da, const smat &K, const vec &rb, const Real gamma, const int max_it, const Real epsilon)
    {
        Natural n  = K.cols();
        Natural m = K.rows();
        Natural num_actions = Da.size();

        mat Q(n, num_actions);

        policy pi(n);
        for (int i = 0; i < n; ++i) pi[i] = random_Natural(0, num_actions-1);

        policy pi_old(pi);
        int max_it2 = max_it;
        if (max_it2 < 0) max_it2 = n;
        Natural it = 0;
        do
        {

            #ifdef VERBOSE
            cout << "PISF: iteration " << it << endl;
            #endif

            smat Ppi(n, m);

            // Policy evaluation
            #ifdef VERBOSE
            cout << "PISF: policy evaluation..." << endl;
            #endif

            for (int i = 0; i < Ppi.rows(); ++i)
            {
                for (smat::InnerIterator it(Da[pi[i]],i); it; ++it) Ppi.insert(i, it.col() ) = it.value();
            }

            Ppi.finalize();

            Ppi = K * Ppi;
            compute_v_iteratively(Ppi, rb, vb, gamma, epsilon);

            // Policy update
            #ifdef VERBOSE
            cout << "PISF: policy update..." << endl;
            #endif

            pi_old = pi;

//            compute_policy(M, Dpi * vb, pi, gamma, (*agt)()); // Slower, but more precise (like in the paper)
            compute_policy(Da, vb, pi, gamma, Q);

            ++it;
        }
        while (pi != pi_old && it < max_it2);
    }

    pt_agent pisf(const vector<smat> &Da, const smat &K, const vec &rb, const Real gamma, const int max_it, const Real epsilon)
    {
        pt_agent agt(new agent(K.cols(), Da.size()));
        vec vb;
        pisf(vb, Da, K, rb, gamma, max_it, epsilon);
        for (size_t a = 0; a < Da.size(); ++a) agt->Q().col(a) = Da[a] * vb;
        return agt;
    }

}

