#ifndef MATRIX_EXT_H
#define MATRIX_EXT_H

#include <fstream>
#include <string>
#include <math.h>

#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

// REMOVE; THIS IS ONLY FOR DEBUGGING
//#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
//#define EIGEN_DONT_ALIGN

// This is relative to where the file Matrix.h is
#define EIGEN_MATRIX_PLUGIN "../../../../EigenMatrixExtension.h"  //melhorar

#include "eigen3/Eigen/Dense"
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET // This is no longer necessary in the new version of Eigen
#include "eigen3/Eigen/Sparse"

#endif
