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

#include "AbsMeasurement.h"
#include "IO.h"

#include <cassert>


namespace genfit {

template<unsigned int dimMeas>
AbsMeasurement<dimMeas>::AbsMeasurement(const SVectorCoord& rawHitCoords, const SMatrixSymCoord& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : rawHitCoords_(rawHitCoords), rawHitCov_(rawHitCov), detId_(detId), hitId_(hitId), trackPoint_(trackPoint)
{}


template<unsigned int dimMeas>
AbsMeasurement<dimMeas>::AbsMeasurement(const AbsMeasurement<dimMeas>& o)
  : rawHitCoords_(o.rawHitCoords_),
    rawHitCov_(o.rawHitCov_),
    detId_(o.detId_),
    hitId_(o.hitId_),
    trackPoint_(o.trackPoint_)
{
  ;
}


template<unsigned int dimMeas>
AbsMeasurement<dimMeas>& AbsMeasurement<dimMeas>::operator=(const AbsMeasurement<dimMeas>&) {
  fputs ("must not call AbsMeasurement::operator=\n",stderr);
  abort();
  return *this;
}


template<unsigned int dimMeas>
void AbsMeasurement<dimMeas>::Print(const Option_t*) const {
  printOut << "genfit::AbsMeasurement, detId = " << detId_ << ". hitId = " << hitId_ << "\n";
  printOut << "Raw hit coordinates: "; rawHitCoords_.Print();
  printOut << "Raw hit covariance: "; rawHitCov_.Print();
}


} /* End of namespace genfit */
