#ifndef AGENT_H
#define AGENT_H

#include "dp.h"

namespace dp
{


    class agent  // agent with a lookup-table action value function
    {

    public:
        agent(const Natural num_states, const Natural num_actions): m_Q(num_states, num_actions) {  }

        mat &Q(void) { return m_Q; }
//        vec &Q(const Natural a) { return m_Q.col(a); this doesn't work }
        Real &Q(const Natural index, const Natural action) { return m_Q(index, action); }
        pt_vec Q(const Natural a)
        {
            pt_vec pv(new vec(m_Q.rows()));
            *pv = m_Q.col(a);
            return pv;
        }

        pt_vec V(void) const
        {
            pt_vec pv(new vec(m_Q.rows()));
            *pv = m_Q.rowwise().maxCoeff();
            return  pv;
        }

        Real V(const Natural index) const { return m_Q.row(index).maxCoeff(); }

        pt_vecn pi(void) const
        {
            pt_policy pi(new policy(num_states()));
            for (Natural i = 0; i < num_states(); ++i) m_Q.row(i).maxCoeff(&((*pi)[i]));
            return pi;
        }

        Natural pi(const Natural index) const
        {
            Natural a;
            m_Q.row(index).maxCoeff(&a);
            return a;
        }



        bool load_from_file(const string filename) { return m_Q.loadFromFile(filename + to_string("_Q")); }
        bool save_to_file(const string filename)   { return m_Q.saveToFile(filename + to_string("_Q")); }

        Natural num_states(void) const { return m_Q.rows(); }
        Natural num_actions(void) const { return m_Q.cols(); }

    protected:
        mat m_Q;
    };



    class agent_continuous_generic // generic agent agent to work in continuous state spaces
    {

    public:
        virtual Natural num_actions(void) const  = 0;
        virtual Natural dim(void) const  = 0;
        virtual ~agent_continuous_generic(void) { }

        virtual Natural pi(const vec &s) const = 0;

        Natural e_pi(const vec &s, const Real epsilon) const
        {
            if (random_Real() > epsilon) return pi(s);
            else return random_Natural(0, num_actions() - 1);
        }


    };


    class agent_continuous : public agent_continuous_generic // agent to work in a continuous state space using state-action value function
    {

    public:
        virtual Real Q(const vec &s, const Natural a) const = 0;
        Real   V(const vec &s) const
        {
            return Q(s, pi(s));
        }

        virtual Natural pi(const vec &s) const
        {
            Real max_Q = Q(s, 0);
            Natural max_a = 0;
            for (Natural a = 1; a < num_actions(); ++a)
            {
                Real value = Q(s,a);
                if (value > max_Q)
                {
                    max_Q = value;
                    max_a = a;
                }
            }
            return max_a;
        }


    };



    class random_agent_continuous : public agent_continuous_generic
    {
    public:
        random_agent_continuous(const Natural dim, const Natural num_actions)
            : m_dim(dim), m_num_actions(num_actions) { }
        virtual Natural num_actions(void) const  {return m_num_actions; }
        virtual Natural dim(void) const  { return m_dim; }
        Natural pi(const vec &s) const { return util::random_Natural(0, m_num_actions - 1); } // this is for efficiency only

        virtual ~random_agent_continuous(void) { }

    private:
        Natural m_dim;
        Natural m_num_actions;

    };


    typedef auto_ptr<agent> pt_agent;
    typedef auto_ptr<agent_continuous_generic> pt_agent_continuous_generic;
    typedef auto_ptr<agent_continuous> pt_agent_continuous;
}

#endif // AGENT_H

