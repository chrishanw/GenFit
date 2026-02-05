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

#include "HMatrixU.h"

#include "IO.h"

#include <cassert>
#include <alloca.h>


namespace genfit {


// 0, 0, 0, 1, 0

const SMatrix15& HMatrixU::getMatrix() const {
  static const double HMatrixContent[5] = {0, 0, 0, 1, 0};

  static const SMatrix15 HMatrix(HMatrixContent);

  return HMatrix;
}


SVector1 HMatrixU::Hv(const SVector5& v) const {
  double* retValArray =(double *)alloca(sizeof(double) * 1);

  retValArray[0] = v(3); // u

  return SVector1(retValArray);
}


SMatrix51 HMatrixU::MHt(const SMatrixSym5& M) const {
  double* retValArray =(double *)alloca(sizeof(double) * 5);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<5; ++i) {
    retValArray[i] = MatArray[i*5 + 3];
  }

  return SMatrix51(retValArray);
}


template<unsigned int nRows>
ROOT::Math::SMatrix<double, nRows, 1> HMatrixU::MHt(const ROOT::Math::SMatrix<double, nRows, 5>& M) const {
  double* retValArray =(double *)alloca(sizeof(double) * M.GetNrows());
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < M.GetNrows(); ++i) {
    retValArray[i] = MatArray[i*5 + 3];
  }

  return ROOT::Math::SMatrix<double, nRows, 1>(retValArray);
}

SMatrixSym1 HMatrixU::HMHt(SMatrixSym5& M) const {
  return SMatrixSym1(M(3,3)):
}


void HMatrixU::Print(const Option_t*) const {
  printOut << "U" << std::endl;
}


} /* End of namespace genfit */
