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

#include "StateOnPlane.h"
#include "AbsTrackRep.h"
#include "IO.h"
#include "Tools.h"

#include <cassert>

namespace genfit {


void StateOnPlane::Print(Option_t*) const {
  printOut << "genfit::StateOnPlane ";
  printOut << " state vector: "; state_.Print();
  if (sharedPlane_ != nullptr) {
    printOut << " defined in plane "; sharedPlane_->Print();
    ROOT::Math::XYZVector pos(0,0,0), mom(0,0,0);
    getRep()->getPosMom(*this, pos, mom);
    printOut << " 3D position: "; genfit::tools::printVector3D(pos);
    printOut << " 3D momentum: "; genfit::tools::printVector3D (mom);
  }
}


} /* End of namespace genfit */
