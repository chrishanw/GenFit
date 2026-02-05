/* Copyright 2013, Technische Universitaet Muenchen,
   Authors: Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "HMatrixPhi.h"

#include "IO.h"

#include <cassert>
#include <alloca.h>
#include <math.h>

namespace genfit {


// 0, 0, 0, cos(phi), sin(phi)


HMatrixPhi::HMatrixPhi(double phi) :
  phi_(phi),
  cosPhi_(cos(phi)),
  sinPhi_(sin(phi))
{
  ;
}

const SMatrix15& HMatrixPhi::getMatrix() const {
  static const double HMatrixContent[5] = {0, 0, 0, cosPhi_, sinPhi_};

  static const SMatrix15 HMatrix(HMatrixContent);

  return HMatrix;
}


SVector1 HMatrixPhi::Hv(const SVector5& v) const {
  double* retValArray =(double *)alloca(sizeof(double) * 1);

  retValArray[0] = cosPhi_*v(3) + sinPhi_*v(4);

  return SVector1(retValArray);
}


SMatrix51 HMatrixPhi::MHt(const SMatrixSym5& M) const {
  double* retValArray =(double *)alloca(sizeof(double) * 5);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<5; ++i) {
    retValArray[i] = cosPhi_*MatArray[i*5 + 3] + sinPhi_*MatArray[i*5 + 4];
  }

  return SMatrix51(retValArray);
}


template<unsigned int nRows>
ROOT::Math::SMatrix<double, nRows, 1> HMatrixPhi::MHt(const ROOT::Math::SMatrix<double, nRows, 5>& M) const {

  double* retValArray =(double *)alloca(sizeof(double) * nRows);
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < nRows; ++i) {
    retValArray[i] = cosPhi_*MatArray[i*5 + 3] + sinPhi_*MatArray[i*5 + 4];
  }

  return  ROOT::Math::SMatrix<double, nRows, 1>(retValArray);
}


SMatrixSym1 HMatrixPhi::HMHt(SMatrixSym5& M) const {
  const double retVal =   cosPhi_ * (cosPhi_*M(3,3) + sinPhi_*M(3,4))
                        + sinPhi_ * (cosPhi_*M(4,3) + sinPhi_*M(4,4));

  return SMatrixSym1(retVal);
}


bool HMatrixPhi::isEqual(const AbsHMatrix& other) const {
  if (dynamic_cast<const HMatrixPhi*>(&other) == nullptr)
    return false;

  return (phi_ == static_cast<const HMatrixPhi*>(&other)->phi_);
}

void HMatrixPhi::Print(const Option_t*) const
{
  printOut << "phi = " << phi_ << std::endl;
}

} /* End of namespace genfit */
