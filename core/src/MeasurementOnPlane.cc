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

#include "MeasurementOnPlane.h"

#include "IO.h"

#include <TClass.h>

namespace genfit {

MeasurementOnPlane::MeasurementOnPlane(const MeasurementOnPlane& other) :
    MeasuredStateOnPlane(other),
    weight_(other.weight_)
{
  hMatrix_.reset(other.hMatrix_->clone());
}


MeasurementOnPlane& MeasurementOnPlane::operator=(MeasurementOnPlane other) {
  swap(other);
  return *this;
}


void MeasurementOnPlane::swap(MeasurementOnPlane& other) {
  MeasuredStateOnPlane::swap(other);
  this->hMatrix_.swap(other.hMatrix_);
  std::swap(this->weight_, other.weight_);
}


void MeasurementOnPlane::Print(Option_t*) const
{
  printOut << "genfit::MeasurementOnPlane, weight = " << weight_ << "\n";
  printOut << " state vector: "; state_.Print();
  printOut << " covariance matrix: "; cov_.Print();
  if (sharedPlane_ != nullptr) {
      printOut << " defined in plane ";
      sharedPlane_->Print();
  }
  printOut << " hMatrix: "; hMatrix_->Print();

}

} /* End of namespace genfit */
