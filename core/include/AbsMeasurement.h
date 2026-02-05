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

#ifndef genfit_AbsMeasurement_h
#define genfit_AbsMeasurement_h

#include "MeasurementOnPlane.h"
#include "AbsHMatrix.h"
// #include "TrackPoint.fwd.h"
#include "AbsTrackPoint.h"

#include <Math/SVector.h>

namespace genfit {

class AbsTrackRep;

/**
 *  @brief Contains the measurement and covariance in raw detector coordinates.
 *
 *  Detector and hit ids can be used to point back to the original detector hits (clusters etc.).
 */
template<unsigned int dimMeas>
class AbsMeasurement {

  using SVectorCoord = ROOT::Math::SVector<double, dimMeas>;
  using SMatrixSymCoord = ROOT::Math::SMatrix<double, dimMeas, dimMeas, ROOT::Math::MatRepSym<double, dimMeas> >;

 public:

  AbsMeasurement() : rawHitCoords_(), rawHitCov_(), detId_(-1), hitId_(-1), trackPoint_(nullptr) {;}
  AbsMeasurement(int nDims) : rawHitCoords_(nDims), rawHitCov_(nDims), detId_(-1), hitId_(-1), trackPoint_(nullptr) {;}
  AbsMeasurement(const SVectorCoord& rawHitCoords, const SMatrixSymCoord& rawHitCov, int detId, int hitId, AbsTrackPoint* trackPoint);

  virtual ~AbsMeasurement() = default;

  //! Deep copy ctor for polymorphic class.
  virtual AbsMeasurement* clone() const = 0;

  AbsTrackPoint* getTrackPoint() const {return trackPoint_;}
  void setTrackPoint(AbsTrackPoint* tp) {trackPoint_ = tp;}

  const SVectorCoord& getRawHitCoords() const {return rawHitCoords_;}
  const SMatrixSymCoord& getRawHitCov() const {return rawHitCov_;}
  SVectorCoord& getRawHitCoords() {return rawHitCoords_;}
  SMatrixSymCoord& getRawHitCov() {return rawHitCov_;}
  int getDetId() const {return detId_;}
  int getHitId() const {return hitId_;}

  //! If the AbsMeasurement is a wire hit, the left/right resolution will be used.
  virtual bool isLeftRightMeasurement() const {return false;}
  virtual int getLeftRightResolution() const {return 0;}

  unsigned int getDim() const {return rawHitCoords_.GetNrows();}

  void setRawHitCoords(const SVectorCoord& coords) {rawHitCoords_ = coords;}
  void setRawHitCov(const SMatrixSymCoord& cov) {rawHitCov_ = cov;}
  void setDetId(int detId) {detId_ = detId;}
  void setHitId(int hitId) {hitId_ = hitId;}


  /**
   * Construct (virtual) detector plane (use state's AbsTrackRep).
   * It's possible to make corrections to the plane here.
   * The state should be defined somewhere near the measurement.
   * For virtual planes, the state will be extrapolated to the POCA to point (SpacepointMeasurement)
   * or line (WireMeasurement), and from this info the plane will be constructed.
   */
  template<unsigned int dim, unsigned int dimAux>
  SharedPlanePtr constructPlane(const StateOnPlane<dim, dimAux>& state) const;

  /**
   * Construct MeasurementOnPlane on plane of the state
   * and wrt the states TrackRep.
   * The state will usually be the prediction or reference state,
   * and has to be defined AT the measurement.
   * The AbsMeasurement will be projected onto the plane.
   * It's possible to make corrections to the coordinates here (e.g. by using the state coordinates).
   * Usually the vector will contain only one element. But in the case of e.g. a WireMeasurement, it will be 2 (left and right).
   */
  template<unsigned int dim, unsigned int dimAux>
  std::vector<genfit::MeasurementOnPlane<dim, dimAux>*> constructMeasurementsOnPlane(const StateOnPlane<dim, dimAux>& state) const;

  /**
   * Returns a new AbsHMatrix object. Caller must take ownership.
   */
  virtual const AbsHMatrix<dimMeas>* constructHMatrix(const AbsTrackRep*) const = 0;

  virtual void Print(const Option_t* = "") const;


 private:
  //! protect from calling assignment operator from outside the class. Use #clone() if you want a copy!
  AbsMeasurement& operator=(const AbsMeasurement&); // default cannot work because TVector and TMatrix = operators don't do resizing

 protected:
  //! protect from calling copy c'tor from outside the class. Use #clone() if you want a copy!
  AbsMeasurement(const AbsMeasurement&);

  SVectorCoord rawHitCoords_;
  SMatrixSymCoord rawHitCov_;
  int detId_; // detId id is -1 per default
  int hitId_; // hitId id is -1 per default

  //! Pointer to TrackPoint where the measurement belongs to
  AbsTrackPoint* trackPoint_; //! No ownership

 public:
  ClassDef(AbsMeasurement, 4)
  // Version history:
  // ver 4: no longer derives from TObject
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsMeasurement_h
