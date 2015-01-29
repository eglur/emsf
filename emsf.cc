#include <iostream>
#include "eigen3/Eigen/Dense"

using namespace std;
using namespace Eigen;

MatrixXd generate_stochastic_matrix(int nrows, int ncols)
{
  int i;
  
  MatrixXd A = MatrixXd::Random(nrows,ncols).cwiseAbs();
  MatrixXd result = A.rowwise().sum();
  for(i = 0; i < ncols; i++) {
    A.col(i) = A.col(i).array() / result.array();
  }
  return A;
}


MatrixXd * generate_stochastic_matrices(int nrows, int ncols, int na)
{
  MatrixXd *matrices = malloc(na * sizeof(MatrixXd));
  return matrices;
}


int main()
{
  int nrows = 3;
  int ncols = 3;

  MatrixXd m = generate_stochastic_matrix(nrows, ncols);

  MatrixXd *matrices = generate_stochastic_matrices(3,3,3);

  cout << "m =" << endl << m << endl;
  return 0;
}
