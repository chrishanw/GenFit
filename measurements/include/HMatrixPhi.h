/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_HMatrixPhi_h
#define genfit_HMatrixPhi_h

#include "AbsHMatrix.h"
#include <SMatrixTypeDefs.h>


namespace genfit {

/**
 * @brief AbsHMatrix implementation for one-dimensional MeasurementOnPlane and RKTrackRep parameterization.
 *
 * For one dimensional measurements which are rotated by phi against U of the DetPlane
 * H = (0, 0, 0, cos(phi), sin(phi))
 */
class HMatrixPhi : public AbsHMatrix<1> {

 public:

  HMatrixPhi(double phi = 0);

  const SMatrix15& getMatrix() const;

  SVector1 Hv(const SVector5& v) const;

  SMatrix51 MHt(const SMatrixSym5& M) const;
  
  template<unsigned int nRows>
  ROOT::Math::SMatrix<double, nRows, 1>  MHt(const ROOT::Math::SMatrix<double, nRows, 5>& M) const;

  SMatrixSym1 HMHt(SMatrixSym5& M) const;

  virtual HMatrixPhi* clone() const {return new HMatrixPhi(*this);}

  virtual bool isEqual(const AbsHMatrix& other) const;

  virtual void Print(const Option_t* = "") const;

  ClassDef(HMatrixPhi,1)

 private:

  double phi_;
  double cosPhi_; //!
  double sinPhi_; //!

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_HMatrixPhi_h
