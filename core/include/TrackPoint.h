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

#ifndef genfit_TrackPoint_h
#define genfit_TrackPoint_h

#include "AbsMeasurement.fwd.h"
#include "AbsTrackPoint.h"
#include "AbsFitterInfo.h"
#include "ThinScatterer.h"
#include <Track.fwd.h>

#include <map>
#include <vector>
#include <memory>


namespace genfit {

/**
 * @brief Object containing AbsMeasurement and AbsFitterInfo objects.
 *
 */
template<unsigned int dimMeas>
class TrackPoint : public AbsTrackPoint {

  using Super = AbsTrackPoint;

 public:

  TrackPoint();
  explicit TrackPoint(Track* track);

  /**
   * @brief Contructor taking list of measurements.
   *
   * AbsMeasurement::setTrackPoint() of each measurement will be called.
   * TrackPoint takes ownership over rawMeasurements.
   */
  TrackPoint(const std::vector< genfit::AbsMeasurement<dimMeas>* >& rawMeasurements, Track* track);

  /**
   * @brief Contructor taking one measurement.
   *
   * AbsMeasurement::setTrackPoint() of the measurement will be called.
   * TrackPoint takes ownership over the rawMeasurement.
   */
  TrackPoint(genfit::AbsMeasurement<dimMeas>* rawMeasurement, Track* track);

  TrackPoint(const TrackPoint<dimMeas>&); // copy constructor
  TrackPoint<dimMeas>& operator=(TrackPoint<dimMeas>); // assignment operator
  void swap(TrackPoint<dimMeas>& other);

  /**
   * custom copy constructor where all TrackRep pointers are exchanged according to the map.
   * FitterInfos with a rep in repsToIgnore will NOT be copied.
   */
  TrackPoint(const TrackPoint<dimMeas>& rhs,
      const std::map<const genfit::AbsTrackRep*, genfit::AbsTrackRep*>& map,
      const std::vector<const genfit::AbsTrackRep*> * repsToIgnore = nullptr);

  virtual ~TrackPoint();

  const std::vector< genfit::AbsMeasurement<dimMeas>* >& getRawMeasurements() const {return rawMeasurements_;}
  AbsMeasurement<dimMeas>* getRawMeasurement(int i = 0) const;
  unsigned int getNumRawMeasurements() const {return rawMeasurements_.size();}
  bool hasRawMeasurements() const {return (rawMeasurements_.size() != 0);}


  //! Takes ownership and sets this as measurement's trackPoint
  void addRawMeasurement(genfit::AbsMeasurement<dimMeas>* rawMeasurement) {assert(rawMeasurement!=nullptr); rawMeasurement->setTrackPoint(this); rawMeasurements_.push_back(rawMeasurement);}
  void deleteRawMeasurements();

  void Print(const Option_t* = "") const;

 private:
  //! Can be more than one, e.g. multiple measurements in the same Si detector, left and right measurements of a wire detector etc.
  std::vector<AbsMeasurement<dimMeas>*> rawMeasurements_; // Ownership

 public:

  ClassDef(TrackPoint,2)
  // Version history:
  // ver 2: no longer derives from TObject

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_TrackPoint_h
