#include <iostream>
#include "eigen3/Eigen/Dense"

using namespace std;

// MatrixXd generate_stochastic_matrix(nrows, ncols)
// {
//   A = NULL;
//   return A;
// }


int main()
{
  MatrixXd m = MatrixXd::Random(3,3);
  cout << "m=" << endl << m << endl;
  return 0;
}
