#ifndef UTIL_H
#define UTIL_H

/*
#ifndef VERBOSE
#define VERBOSE
#endif
*/

#ifndef NDEBUG
#define NDEBUG //faster
#endif

#include <string>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <memory>
#include <iostream>
#include <limits>
#include <sys/time.h>


#include "matrix_ext.h"


namespace util
{
    typedef int Natural; // the indices of Eigen are ints...
    typedef float Real; // change to double for RL-Glue?

    typedef Eigen::Matrix<util::Real, Eigen::Dynamic, Eigen::Dynamic> mat;
    typedef Eigen::Matrix<util::Real, Eigen::Dynamic, 1> vec;
    typedef Eigen::Matrix<util::Real, 1, Eigen::Dynamic> row_vec;
    typedef Eigen::Matrix<util::Natural, Eigen::Dynamic, Eigen::Dynamic> matn;
    typedef Eigen::Matrix<util::Natural, 1, Eigen::Dynamic> vecn;
    typedef std::auto_ptr<vec> pt_vec;
    typedef std::auto_ptr<mat> pt_mat;
    typedef std::auto_ptr<vecn> pt_vecn;
    typedef Eigen::SparseMatrix<util::Real,Eigen::RowMajor> smat;
    typedef std::auto_ptr<smat> pt_smat;
    typedef Eigen::DynamicSparseMatrix<util::Real,Eigen::RowMajor> dsmat;
    typedef Eigen::SparseMatrix<util::Real,Eigen::ColMajor> smat_col;
    typedef Eigen::SparseVector<util::Real,Eigen::ColMajor> svec;
    typedef Eigen::SparseVector<util::Natural,Eigen::ColMajor> svecn;
    typedef Eigen::SparseVector<util::Real,Eigen::RowMajor> srow_vec;
    const Real PI = 3.141592;
    const Real inf = std::numeric_limits<Real>::max(); // change this?
    const Real minf = -std::numeric_limits<Real>::infinity(); // can be diff -inf ?

    inline void random_Seed(unsigned long int s)
    {
        std::srand(s);
    }


    inline Natural random_Natural(const Natural lower_bound = 0, const Natural upper_bound = 1)
    {
        return lower_bound + static_cast<Natural>(rand() * (upper_bound - lower_bound + 1.0) / (RAND_MAX + 1.0) );
    }

    inline Real random_Real(const Real lower_bound = 0.0, const Real upper_bound = 1.0)
    {
        return (lower_bound + Real(std::rand())/Real(RAND_MAX) * (upper_bound - lower_bound));
    }

    inline Real random_Normal(const Real mean = 0, const Real sd = 1)
    {
        Real x1, x2, w;

        do
        {
            x1 = 2.0 * random_Real() - 1.0;
            x2 = 2.0 * random_Real() - 1.0;
            w = x1 * x1 + x2 * x2;
        }
        while ( w >= 1.0 );

        return mean + sd * (x1 * sqrt( (-2.0 * log( w ) ) / w ));
    }


    inline Real normal(const Real x, const Real mu = 0, const Real sd = 1)
    {
        return 1.0/(sqrt(2* PI)* sd) * exp(-(pow(x - mu, 2) / (2 * pow(sd,2))));
    }


    template<class T> inline std::string to_string(const T& t)
    {
        try
        {
            std::stringstream ss;
            ss << t;
            return ss.str();
        }
        catch ( ... )
        {
            return "";
        }
    }

    template<class T> inline bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&) = std::dec)
    {
        std::istringstream iss(s);
        return !(iss >> f >> t).fail();
    }


    inline pt_vec seq_length(const Real lower_bound, const Real upper_bound, const Natural length)
    {
        pt_vec s( new vec(length) );
        for (Natural i = 0; i < s->size(); ++i) (*s)[i] = lower_bound + i * (upper_bound - lower_bound) / (length - 1);
        return s;
    }

    inline pt_vecn seq_length(const Natural lower_bound, const Natural upper_bound, const Natural length)
    {
        pt_vecn s( new vecn(length) );
        for (Natural i = 0; i < s->size(); ++i) (*s)[i] = lower_bound + i * (upper_bound - lower_bound) / (length - 1);
        return s;
    }

    inline pt_vecn sample(const Natural lower_bound, const Natural upper_bound, const Natural n, const bool replacement = false)
    { /// use list for inds??
        pt_vecn pv( new vecn(n) );
        vecn inds(upper_bound - lower_bound + 1);
        for (Natural i = 0; i < inds.size(); ++i) inds[i] = lower_bound + i;
        Natural size = inds.size();
        for (Natural i = 0; i < n; ++i)
        {
            Natural j = random_Natural(0, size-1);
            (*pv)[i] = inds[j];
            if (!replacement)
            {
                for (Natural k = j; k < size - 1; ++k) inds[k] = inds[k + 1];
                --size;
            }
        }
        return pv;
    }


    /* std::isnan and std::isinf are not portable */
    template<typename T> inline bool is_nan(T const& x)
    {
        return (x != x);
    }

    template <typename T> inline bool is_inf(T const& x)
    {
        return std::numeric_limits<T>::has_infinity && (x == std::numeric_limits<T>::infinity() or x == -std::numeric_limits<T>::infinity());
    }

    inline bool is_valid(const Real v) { return !(is_inf(v) || is_nan(v)); }

    inline bool is_valid(const vec &v)
    {
        bool valid = true;
        Natural i = 0;
        while (i < v.size() && valid) valid = is_valid(v[i++]);
        return valid;
    }


    inline Real mean(const vec &v)
    {
        return v.sum() / v.size();
    }

    inline Real mean(const vec &v, const vecn &inds)
    {
        Real sum = v[inds[0]];
        for (Natural i = 1; i < inds.size(); ++i) sum += v[inds[i]];
        return sum / inds.size();
    }

    inline Real var(const vec &v)
    {
        Real m = mean(v);
        Real var = 0;
        for (Natural i = 0; i < v.size(); ++i) var += pow(v[i] - m, 2);
        Real den = v.size() > 1 ? v.size() - 1 : 1;
        return var / den;
    }


    inline Real var(const vec &v, const vecn &inds)
    {
        Real m = mean(v, inds);
        Real var = 0;
        for (Natural i = 0; i < inds.size(); ++i) var += pow(v[inds[i]] - m, 2);
        Real den = inds.size() > 1 ? inds.size() - 1 : 1;
        return var / den;
    }


    inline Real max(const vec &v, const vecn &inds)
    {
        Real vmax = v[inds[0]];
        for (Natural i = 1; i < inds.size(); ++i)
        {
            Real value = v[inds[i]];
            if (value > vmax) vmax = value;
        }
        return vmax;
    }

    inline Real min(const vec &v, const vecn &inds)
    {
        Real vmin = v[inds[0]];
        for (Natural i = 1; i < inds.size(); ++i)
        {
            Real value = v[inds[i]];
            if (value < vmin) vmin = value;
        }
        return vmin;
    }


    inline Real sd(const vec &v)
    {
        return sqrt(var(v));
    }

    inline Real sd(const vec &v, const vecn &inds)
    {
        return sqrt(var(v, inds));
    }

    inline pt_vec normalize(const vec &v)
    {
        pt_vec pv(new vec(v.size()));
        *pv = (v.array() - mean(v)) / sd(v);
        return pv;
    }

    inline void swap(vec &v, Natural i, Natural j)
    {
        Real tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
    }


}

#endif

