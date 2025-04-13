/**
 * Author: Christian Wessel
 * Creation date: 25.03.2025
 */

#pragma once

#include <Rtypes.h>

namespace genfit {

/**
 * @brief Abstract TrackPoint class.
 *
 */
class AbsTrackPoint {

 public:

  AbsTrackPoint(unsigned int dim = 1) : dim_(dim) {};

  unsigned int getDim() const { return dim_;}

 private:
  unsigned int dim_;

 public:

  ClassDef(AbsTrackPoint,1)

};

} /* End of namespace genfit */
/** @} */
