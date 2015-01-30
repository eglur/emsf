#include "mdp.h"

namespace dp
{

    void generic_mdp::increase_size(const Natural num_states)
    {

        for (Natural a = 0; a < num_actions(); ++a)
        {
            m_r[a].conservativeResize(num_states);
           // Fill out the new rewards with zeros
            m_r[a].segment(m_ns, num_states - m_ns) = vec::Zero(num_states - m_ns);
        }

        m_ns = num_states;
    }



    inline void mdp::resize_matrices(const Natural num_states, const Natural num_actions)
    {
        m_P.resize(num_actions);
        m_r.resize(num_actions);
        for (Natural a = 0; a < num_actions; ++a)
        {
            m_P[a] = mat::Zero(num_states, num_states);
            m_r[a] = vec::Zero(num_states);
        }
    }


    void mdp::increase_size(const Natural num_states)
    {

        generic_mdp::increase_size(num_states);
        for (Natural a = 0; a < num_actions(); ++a)
        {
            m_P[a].conservativeResize(num_states, num_states);
            // Fill out the transitions to the new states with zeros
            m_P[a].block(0, m_ns, m_ns, num_states - m_ns) = mat::Zero(m_ns, num_states - m_ns);
            // Fill out the transitions from the new states with zeros
            m_P[a].block(m_ns, 0, num_states - m_ns, num_states) = mat::Zero(num_states - m_ns, num_states);
            // Fill out the new rewards with zeros
        }

    }


    bool mdp::load_from_file(const string filename)
    {
        ifstream ifs((filename + "_info").c_str());

        if (!ifs.good()) return false;

        ifs >> m_ns;
        ifs >> m_na;

        ifs.close();

        resize_matrices(m_ns, m_na);

        for (Natural a=0; a < m_na; ++ a)
        {
            bool ok = m_P[a].loadFromFile((filename+"_Pa"+ util::to_string(a)).c_str(), m_ns, m_ns);
            if (ok) ok = m_r[a].loadFromFile((filename+"_ra"+ util::to_string(a)).c_str(), m_ns, 1);
            if (!ok) return false;
        }
        return true;
    }


    bool mdp::save_to_file(const string filename, const Natural precision)
    {
        ofstream ofs((filename + "_info").c_str());

        if (!ofs.good()) return false;

        ofs << m_ns << " " << m_na << endl;
        ofs.close();

        for (Natural a=0; a < m_na; ++ a)
        {
            bool ok = m_P[a].saveToFile((filename+"_Pa"+ to_string(a)).c_str(), precision);
            if (ok) ok = m_r[a].saveToFile((filename+"_ra"+ to_string(a)).c_str(), precision);
            if (!ok) return false;
        }
        return true;
    }


    void smdp::increase_size(const Natural num_states)
    {
        generic_mdp::increase_size(num_states);
        for (Natural a = 0; a < m_na; ++a)
        {
            smat *T = new smat(num_states, num_states);
            T->reserve(m_P[a]->nonZeros());
            for (Natural i = 0; i < m_P[a]->rows(); ++i)
            {
                for (smat::InnerIterator it(*(m_P[a]),i); it; ++it) T->insert(i, it.col() ) = it.value();

            }
            T->finalize();
            delete m_P[a];
            m_P[a] = T;

        }
    }



    pt_mdp normal_mdp(const Natural num_states, const Natural num_actions, const Natural num_centers, const Real sd, const Real noise, Real threshold)
    {
        pt_mdp pm = pt_mdp(new mdp(num_states, num_actions));


        for (Natural a = 0; a < num_actions; ++a)
        {
            // Generate rows associated with centers
            pt_vecn pcenters = util::sample(0, num_states - 1, num_centers, false);


            for (Natural i = 0; i < num_centers; ++i)
            {
                for (Natural j = 0; j < num_states; ++j) pm->P(a)(i, j) = util::normal(j, (*pcenters)[i], sd) + util::random_Real(0, noise);
                // Normalize
                Real sum = pm->P(a).row(i).sum();
                for (Natural j = 0; j < num_states; ++j) pm->P(a)(i, j) /= sum;

            }

            // Generate remaining rows (using the same centers)
            pt_vecn prem = util::sample(0, num_centers - 1, num_states - num_centers, true);

            for (Natural i = 0; i < prem->size(); ++i)
            {
                Natural c = (*prem)[i], r = num_centers + i;

                for (Natural j = 0; j < num_states; ++j) pm->P(a)(r, j) = util::normal(j, (*pcenters)[c], sd) + util::random_Real(0, noise);
                // Normalize
                Real sum = pm->P(a).row(r).sum();
                for (Natural j = 0; j < num_states; ++j) pm->P(a)(r, j) /= sum;
            }


            // Generate rewards
            for (Natural i = 0; i < num_states; ++i) pm->r(a)[i] = random_Real();
        }

        /*
        Natural goal = random_Natural(0, num_states - 1);
        for (Natural a = 0; a < num_actions; ++a) pm->r(a)[goal] = 100.0;
        */

        return pm;
    }


}
