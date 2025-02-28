//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "FactorizedRankTwoTensor.h"
#include "RaccoonUtils.h"
#include "ADReal.h"
#include "MooseUtils.h"
#include "RankTwoTensorForward.h"
#include "RankTwoTensorImplementation.h"
#include "metaphysicl/raw_type.h"

namespace RaccoonUtils
{

ADReal
Macaulay(const ADReal x, const bool deriv)
{
  if (deriv)
    return x > 0 ? 1 : 0;
  return 0.5 * (x + std::abs(x));
}

std::vector<ADReal>
Macaulay(const std::vector<ADReal> & v, const bool deriv)
{
  std::vector<ADReal> m = v;
  for (auto & x : m)
    x = Macaulay(x, deriv);
  return m;
}

ADRankTwoTensor
spectralDecomposition(const ADRankTwoTensor & r2t)
{
  ADRankTwoTensor eigvecs;
  std::vector<ADReal> eigvals(LIBMESH_DIM);
  r2t.symmetricEigenvaluesEigenvectors(eigvals, eigvecs);

  ADRankTwoTensor eigvals_pos;
  eigvals_pos.fillFromInputVector(Macaulay(eigvals));
  return eigvecs * eigvals_pos * eigvecs.transpose();
}

ADRankTwoTensor
log(const ADRankTwoTensor & r2t)
{
  FactorizedRankTwoTensor A = MetaPhysicL::raw_value(r2t);
  ADRankTwoTensor I2;
  I2.setToIdentity();
  int s = 0;

  const double tol = 0.5;
  while (MetaPhysicL::raw_value(A.get() - I2).norm() > tol)
  {
    A = MathUtils::sqrt(A);
    s++;
  }

  ADRankTwoTensor X = A.get() - I2;

  ADRankTwoTensor logA = X;
  ADRankTwoTensor term = X;
  const double series_tol = 1e-8;
  int k = 2;

  while (MetaPhysicL::raw_value(term).norm() > series_tol && k < 100)
  {
    term = term * X;

    double coeff = ((k % 2 == 0) ? -1 : 1) / double(k);

    logA += coeff * term;

    k++;
  }

  ADRankTwoTensor logF = logA * std::pow(2, s);

  return logF;
  // FactorizedRankTwoTensor A = MetaPhysicL::raw_value(r2t);
  // A = MathUtils::log(A);
  // A = MathUtils::sqrt(const FactorizedRankTwoTensorTempl<T> &A)
  // ADRankTwoTensor B;
  // B = A.get();
  // return B;

  // std::vector<ADReal> d;
  // ADRankTwoTensor V, D;
  // r2t.symmetricEigenvaluesEigenvectors(d, V);
  // for (auto & di : d)
  //   di = std::log(di);
  // D.fillFromInputVector(d);
  // return V * D * V.transpose();
}

ADRankTwoTensor
exp(const ADRankTwoTensor & r2t)
{

  // FactorizedRankTwoTensor A = MetaPhysicL::raw_value(r2t);
  // A = MathUtils::exp(A);
  // ADRankTwoTensor B;
  // B = A.get();
  // return B;
  // int accuracy = 10;

  // scaling
  // int N = 4;

  // // M_small = M/(2^N);
  // ADRankTwoTensor A_small = r2t / std::pow(2, N);

  // ADRankTwoTensor expA;
  // ADRankTwoTensor term;
  // expA.setToIdentity();
  // term.setToIdentity();

  // ADRankTwoTensor A_power = A_small;

  // double factorial = 1;

  // const int max_iterations = 100;
  // double tol = 1e-12;
  // auto termnorm = MetaPhysicL::raw_value(term).norm();
  // auto expnorm = MetaPhysicL::raw_value(expA).norm();

  // for (int i = 1; i < max_iterations; i++)
  // {
  //   factorial = factorial * i;

  //   // term = A_small^i / i!
  //   term = (A_power) / factorial;
  //   expA += term;
  //   termnorm = MetaPhysicL::raw_value(term).norm();
  //   expnorm = MetaPhysicL::raw_value(expA).norm();

  //   if (termnorm < tol * expnorm || termnorm < tol)
  //     break;
  //   // Check for convergence using the tensor norm (if available)
  //   // if (MetaPhysicL::raw_value(term).norm() < 1e-12)
  //   //   ;

  //   A_power = A_power * A_small;
  // }

  // for (int i = 0; i < N; i++)
  // {
  //   expA = expA * expA;
  // }

  // return expA;
  std::vector<ADReal> d;
  ADRankTwoTensor V, D;
  r2t.symmetricEigenvaluesEigenvectors(d, V);
  for (auto & di : d)
    di = std::exp(di);
  D.fillFromInputVector(d);
  return V * D * V.transpose();
}

} // end namespace MooseUtils
