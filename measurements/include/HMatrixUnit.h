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

#ifndef genfit_HMatrixUnit_h
#define genfit_HMatrixUnit_h

#include "AbsHMatrix.h"
#include <SMatrixTypeDefs.h>


namespace genfit {

/**
 * @brief AbsHMatrix implementation for 5-dimensional MeasurementOnPlane and RKTrackRep parameterization.
 *
 * H = (1, 0, 0, 0, 0)
 *     (0, 1, 0, 0, 0)
 *     (0, 0, 1, 0, 0)
 *     (0, 0, 0, 1, 0)
 *     (0, 0, 0, 0, 1)
 */
class HMatrixUnit : public AbsHMatrix<5> {

 public:

  HMatrixUnit() {;}

  const SMatrix55& getMatrix() const;

  SVector5 Hv(const SVector5& v) const {return v;}

  SMatrix55 MHt(const SMatrixSym5& M) const {return SMatrix55(M);}

  template<unsigned int nRows>
  ROOT::Math::SMatrix<double, nRows, 5> MHt(const ROOT::Math::SMatrix<double, nRows, 5>& M) const {return M;}

  SMatrixSym5 HMHt(SMatrixSym5& M) const {return M;}

  virtual HMatrixUnit* clone() const {return new HMatrixUnit(*this);}

  virtual bool isEqual(const AbsHMatrix& other) const {return (dynamic_cast<const HMatrixUnit*>(&other) != nullptr);}

  virtual void Print(const Option_t* = "") const;

  ClassDef(HMatrixUnit,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_HMatrixUnit_h
