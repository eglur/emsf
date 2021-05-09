#include <iostream>
#include "policy_iteration.h"

using namespace std;
using namespace dp;

inline void normalize(mat &A)
{
    for (Natural i = 0; i < A.rows(); ++i)
    {
        Real sum = A.row(i).sum();
        if (sum > 0) for (Natural j = 0; j < A.cols(); ++j) A(i,j) /= sum;
    }
}

inline void randomize(mat &A)
{
    for (Natural i = 0; i < A.rows(); ++i)
    {
       for (Natural j = 0; j < A.cols(); ++j) A(i,j) = util::random_Real();
    }
}


int main()
{
    const Natural n = 100;   // Number of regular states
    const Natural m = 10;    // Number of artificial states
    const Natural c = 10;    // Number of episodes
    const Natural T = 1000;  // Number of transitions per episode
    const Natural na = 4;    // Number of actions
    const Real gamma = 0.99; // Discount factor

    // MDP
    mdp M(n, na);

    // Create a factorizable MDP
    for (Natural a = 0; a < na; ++a)
    {
        mat D(n ,m);
        randomize(D);
        normalize(D);

        mat K(m, n);
        randomize(K);
        normalize(K);

        M.P(a) = D * K;
        for (Natural a = 0; a < na; ++a)
            for (Natural i = 0; i < n; ++i) M.r(a)[i] = util::random_Real();
    }

    // Solve the MDP using policy iteration
    pt_agent agt = policy_iteration(M, gamma);

    // Print value function
    cout << *agt->V() << endl;

    exit(EXIT_SUCCESS);

}

