#ifndef MDP_H
#define MDP_H

#include <vector>
#include <string>
#include <fstream>
#include "dp.h"

namespace dp{

class generic_mdp
{
    public:
        generic_mdp(const Natural num_states, const Natural num_actions) : m_ns(num_states), m_na(num_actions), m_r(num_actions)
        {
            for (Natural a = 0; a < num_actions; ++a) m_r[a] = vec::Zero(num_states);
        }

        virtual ~generic_mdp(void) { }

        Natural num_actions(void) const { return m_na; }
        Natural num_states(void) const { return m_ns; }
        vec &r(const Natural index) { return m_r[index]; }
        const vec &r(const Natural index) const { return m_r[index]; }
        Real r(const Natural a, const Natural i) const { return m_r[a][i]; }

        virtual void increase_size(const Natural num_states);

    protected:
        Natural m_ns, m_na;
        std::vector<vec> m_r;
};



class mdp : public generic_mdp
{
    public:
        mdp(const Natural num_states, const Natural num_actions) : generic_mdp(num_states, num_actions), m_P(num_actions)
        {
            for (Natural a = 0; a < num_actions; ++a) m_P[a] = mat::Zero(num_states, num_states);
        }

//        mdp(const string filename) { load_from_file(filename); }

        mat &P(const Natural index) { return m_P[index]; }
        const mat &P(const Natural index) const { return m_P[index]; }
        Real p(const Natural a, const Natural i, const Natural j) const { return m_P[a](i,j); }

        bool load_from_file(const string filename);
        bool save_to_file(const string filename, const Natural precision = 6);

        virtual void increase_size(const Natural num_states);

        typedef mat mat_type;

    protected:
        std::vector<mat_type> m_P;
        void resize_matrices(const Natural num_states, const Natural num_actions);
};


class smdp : public generic_mdp //sparse MDP
{
    public:

        smdp(const Natural num_states, const Natural num_actions): generic_mdp(num_states, num_actions), m_P(num_actions)
        {
            for (Natural a = 0; a < m_na; ++a) m_P[a] = new smat(num_states, num_states);
        }

        virtual ~smdp(void)
        {
            for (Natural a = 0; a < m_na; ++a) delete m_P[a];
        }

        smat &P(const Natural index) { return *(m_P[index]); }
        const smat &P(const Natural index) const { return *(m_P[index]); }

        virtual void increase_size(const Natural num_states);
        typedef smat mat_type;

    protected:
        std::vector<mat_type*> m_P; // using pointers because the assignment operator (=) seems to be buggy
};


    typedef std::auto_ptr<mdp> pt_mdp;
    typedef std::auto_ptr<smdp> pt_smdp;

    pt_mdp normal_mdp(const Natural num_states, const Natural num_actions, const Natural num_centers, const Real sd = 1, const Real noise = 0, Real threshold = 1e-3);

}

#endif // MDP_H
