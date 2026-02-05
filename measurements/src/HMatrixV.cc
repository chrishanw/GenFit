/* Copyright 2013, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Johannes Rauch, Tobias Schlüter

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

#include "HMatrixV.h"

#include "IO.h"

#include <cassert>
#include <alloca.h>

namespace genfit {


// 0, 0, 0, 0, 1

const SMatrix15& HMatrixV::getMatrix() const {
  static const double HMatrixContent[5] = {0, 0, 0, 0, 1};

  static const SMatrix15 HMatrix(HMatrixContent);

  return HMatrix;
}


SVector1 HMatrixV::Hv(const SVector5& v) const {
  double* retValArray =(double *)alloca(sizeof(double) * 1);

  retValArray[0] = v(4); // v

  return SVector1(retValArray);
}


SMatrix51 HMatrixV::MHt(const SMatrixSym5& M) const {
  double* retValArray =(double *)alloca(sizeof(double) * 5);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<5; ++i) {
    retValArray[i] = MatArray[i*5 + 4];
  }

  return SMatrix51(retValArray);
}


template<unsigned int nRows>
ROOT::Math::SMatrix<double, nRows, 1> HMatrixV::MHt(const ROOT::Math::SMatrix<double, nRows, 5>& M) const {
  double* retValArray =(double *)alloca(sizeof(double) * nRows);
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < nRows; ++i) {
    retValArray[i] = MatArray[i*5 + 4];
  }

  return ROOT::Math::SMatrix<double, nRows, 1>(retValArray);
}


SMatrixSym1 HMatrixV::HMHt(SMatrixSym5& M) const {
  return SMatrixSym1(M(4,4)):
}


void HMatrixV::Print(const Option_t*) const {
  printOut << "V" << std::endl;
}

} /* End of namespace genfit */
