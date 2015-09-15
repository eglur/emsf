/********************************
	Class Tools
	Author: Eduardo Krempser 

	Define basic structures to 
	handle and generate stochastic matrices

	Created: 2015 - July, 30
*********************************/

#ifndef TOOLS.H
#define TOOLS.H

#include "iostream"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/Core"
#include "stdlib.h"
#include <vector>
#include <ctime>
#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>

#include "omp.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "time.h"

using namespace Eigen;
using namespace std;

class Tools{
 private:
  //Variables to handle a Dirichlet distribution
  static mt19937 rng;
  static gsl_rng *rngDirichlet;
  static float alphaParameter;

 public:
  //Initialize random variables
  static void init();

  //Define a value to random variable seed
  static void setSeed(int seed);

  //Define a alpha value to Dirichlet Distribution
  static void setAlphaParameter(float alpha);

  //Generate stochastic matrix 
  static SparseMatrix<float, RowMajor> generateStochasticMatrixDirichlet(int nrows, int ncols);

  //Sort a number - normal distribution	
  static int sample(int start, int end);
  static float sample();

  //Select a value from a dist
  static float sampleFromDist(SparseMatrix<float, RowMajor> D);

  //Calc errors (different norms)
  static float frobenius(SparseMatrix<float, RowMajor> P, SparseMatrix<float, RowMajor> Pt);

  static float frobeniusWeighted(SparseMatrix<float, RowMajor> P, SparseMatrix<float, RowMajor> Pt, VectorXf w);

  static float KL(SparseMatrix<float, RowMajor> P, SparseMatrix<float, RowMajor> Pt);

  static float KLWeighted(SparseMatrix<float, RowMajor> P, SparseMatrix<float, RowMajor> Pt, VectorXf w);
};

#endif
