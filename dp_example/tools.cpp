#include "Tools.h"

float Tools::alphaParameter = 0.5;
mt19937 Tools::rng;
gsl_rng* Tools::rngDirichlet;

//Init seed
void Tools::init(){
  rng.seed(static_cast<unsigned> (std::time(0)));

  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rngDirichlet = gsl_rng_alloc(T);
  gsl_rng_set(rngDirichlet, std::time(0));
}

//Define initial seed
void Tools::setSeed(int seed){
  rng.seed(static_cast<unsigned> (seed));
  gsl_rng_set(rngDirichlet, seed);

  gsl_rng_set(rngDirichlet, seed);
}

//Define a Dirichlet Alpha Parameter
void Tools::setAlphaParameter(float alpha){
  alphaParameter = alpha;
}

//Generate a stochastic matrix by a Dirichlet Distribution
SparseMatrix<float, RowMajor> Tools::generateStochasticMatrixDirichlet(int nrows, int ncols){
  SparseMatrix<float, RowMajor> A(nrows, ncols);

  double *alpha = new double[ncols];
  double *theta = new double[ncols];

  for(int j = 0; j < ncols; j++){
    alpha[j] = alphaParameter;
  }


  //Dirichlet
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      theta[j] = 1;
    }
    gsl_ran_dirichlet(rngDirichlet, ncols, alpha, theta);
    for(int j = 0; j < ncols; j++){
      if (theta[j] > 1e-6){
        A.coeffRef(i, j) = (float)theta[j];
      }
    }
    A.row(i) = A.row(i) / A.row(i).sum();
  }

  //free memory
  delete alpha;
  delete theta;

  return A;
}

//Sort a number (uniform distribution) between (start, end)
int Tools::sample(int start, int end){
  float v = sample();
  return (int)(start + (end-start)*v);
}

//Sort a number (uniform distribution) between (0,1)
float Tools::sample(){
  std::uniform_real_distribution<float> distribution(0.0,1.0);

  return distribution(rng);
}

//Select a value from a dist
float Tools::sampleFromDist(SparseMatrix<float, RowMajor> D){
  float v = sample();

  float vv = 0.0;
  int ind = 0;

  for(ind = 0; ind < D.cols() && vv < v; ind++){
    vv = vv + D.coeff(0, ind);
  }

  ind--;
  return ind;
}

//Calc erro by different norms
float Tools::frobenius(SparseMatrix<float, RowMajor> P, SparseMatrix<float, RowMajor> Pt){
  float sse = 0.0;

  //for (a in 1:dim(P)[3]) sse <- sse + sum((P[,,a] - Pt[,,a])^2)
  SparseMatrix<float, RowMajor> mAux = P - Pt;
  float sum = 0;
  for (int k=0; k<mAux.outerSize(); ++k){
    for (SparseMatrix<float, RowMajor>::InnerIterator it(mAux,k); it; ++it){
      sum += it.value() * it.value();
    }
  }
  sse += sum;


  return sse;
}

float Tools::frobeniusWeighted(SparseMatrix<float, RowMajor> P, SparseMatrix<float, RowMajor> Pt, VectorXf w){
  float sse = 0.0;

  //for (a in 1:dim(P)[3]) sse <- sse + sum((P[,,a] - Pt[,,a])^2)
  SparseMatrix<float, RowMajor> mAux = P - Pt;
  float sum = 0;
  for (int k=0; k<mAux.outerSize(); ++k){
    for (SparseMatrix<float, RowMajor>::InnerIterator it(mAux,k); it; ++it){
      mAux.coeffRef(it.row(), it.col()) = it.value() * it.value();
    }
  }

  for(int i = 0; i < mAux.rows(); i++){
    mAux.row(i) = w.coeff(i) * mAux.row(i);
    sum += mAux.row(i).sum();
  }
  sse += sum;


  return sse;
}

float Tools::KL(SparseMatrix<float, RowMajor> P, SparseMatrix<float, RowMajor> Pt){
  float sse = 0.0;


  SparseMatrix<float, RowMajor> mAux(P.rows(), P.cols());
  for(int i = 0; i < mAux.rows(); i++){
    for(int j = 0; j < mAux.cols(); j++){
      mAux.coeffRef(i,j) = P.coeff(i,j) * log(P.coeff(i,j)/Pt.coeff(i,j)) - P.coeff(i,j) + Pt.coeff(i,j);
    }
  }

  float sum = 0;
  for (int k=0; k<mAux.outerSize(); ++k){
    for (SparseMatrix<float, RowMajor>::InnerIterator it(mAux,k); it; ++it){
      sum += it.value();
    }
  }
  sse += sum;


  return sse;
}

float Tools::KLWeighted(SparseMatrix<float, RowMajor> P, SparseMatrix<float, RowMajor> Pt, VectorXf w){
  float sse = 0.0;

  SparseMatrix<float, RowMajor> mAux(P.rows(), P.cols());
  for(int i = 0; i < mAux.rows(); i++){
    for(int j = 0; j < mAux.cols(); j++){
      mAux.coeffRef(i,j) = P.coeff(i,j) * log(P.coeff(i,j)/Pt.coeff(i,j)) - P.coeff(i,j) + Pt.coeff(i,j);
    }
  }

  float sum = 0;
  /*for (int k=0; k<mAux.outerSize(); ++k){
    for (SparseMatrix<float, RowMajor>::InnerIterator it(mAux,k); it; ++it){
    sum += it.value();
    }
    }*/
  //cout << "w: " << w << endl;
  for(int i = 0; i < mAux.rows(); i++){
    mAux.row(i) = w.coeff(i) * mAux.row(i);
    sum += mAux.row(i).sum();
  }

  sse += sum;


  return sse;
}
