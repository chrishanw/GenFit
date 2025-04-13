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

#ifndef genfit_StateOnPlane_h
#define genfit_StateOnPlane_h

#include "SharedPlanePtr.h"
#include "AbsTrackRep.h"
#include <SMatrixTypeDefs.h>

#include <TVectorD.h>

#include <cassert>


namespace genfit {

/**
 * @brief A state with arbitrary dimension defined in a DetPlane.
 *
 * The dimension and meaning of the #state_ vector are defined by the track parameterization of the #rep_.
 * #sharedPlane_ is a shared_pointer, the ownership over that plane is shared between all StateOnPlane objects defined in that plane.
 * The definition of the state is bound to the TrackRep #rep_. Therefore, the StateOnPlane contains a pointer to a AbsTrackRep.
 * It will provide functionality to extrapolate it and translate the state it into cartesian coordinates.
 * Shortcuts to all functions of the AbsTrackRep which use this StateOnPlane are also provided here.
 */
template<unsigned int dim, unsigned int dimAux = 0>
class StateOnPlane {

  using SVectorState = ROOT::Math::SVector<double, dim>;
  using SVectorAux = ROOT::Math::SVector<double, dimAux>;

 public:

  StateOnPlane(const genfit::StateOnPlane<dim, dimAux>&) = default;

  StateOnPlane(const genfit::AbsTrackRep* rep = nullptr);
  //! state is defined by the TrackReps parameterization
  StateOnPlane(const SVectorState& state, const SharedPlanePtr& plane, const AbsTrackRep* rep);
  StateOnPlane(const SVectorState& state, const SharedPlanePtr& plane, const AbsTrackRep* rep, const SVectorAux& auxInfo);

  StateOnPlane& operator=(StateOnPlane other);
  void swap(StateOnPlane& other); // nothrow

  virtual ~StateOnPlane() {}
  virtual StateOnPlane* clone() const {return new StateOnPlane(*this);}

  const SVectorState& getState() const {return state_;}
  SVectorState& getState() {return state_;}
  const SVectorAux& getAuxInfo() const {return auxInfo_;}
  SVectorAux& getAuxInfo() {return auxInfo_;}
  const SharedPlanePtr& getPlane() const {return sharedPlane_;}
  const AbsTrackRep* getRep() const {return rep_;}

  void setState(const SVectorState& state) {state_ = state;}
  void setPlane(const SharedPlanePtr& plane) {sharedPlane_ = plane;}
  void setStatePlane(const SVectorState& state, const SharedPlanePtr& plane) {state_ = state; sharedPlane_ = plane;}
  void setAuxInfo(const SVectorAux& auxInfo) {auxInfo_ = auxInfo;}
  void setRep(const AbsTrackRep* rep) {rep_ = rep;}

  // Shortcuts to TrackRep functions
  double extrapolateToPlane(const SharedPlanePtr& plane,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToPlane(*this, plane, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToLine(const ROOT::Math::XYZVector& linePoint,
        const ROOT::Math::XYZVector& lineDirection,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToLine(*this, linePoint, lineDirection, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToPoint(const ROOT::Math::XYZVector& point,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToPoint(*this, point, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToPoint(const ROOT::Math::XYZVector& point,
        const TMatrixDSym& G, // weight matrix (metric)
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToPoint(*this, point, G, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToCylinder(double radius,
        const ROOT::Math::XYZVector& linePoint = ROOT::Math::XYZVector(0.,0.,0.),
        const ROOT::Math::XYZVector& lineDirection = ROOT::Math::XYZVector(0.,0.,1.),
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToCylinder(*this, radius, linePoint, lineDirection, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToCone(double openingAngle,
        const ROOT::Math::XYZVector& conePoint = ROOT::Math::XYZVector(0.,0.,0.),
        const ROOT::Math::XYZVector& coneDirection = ROOT::Math::XYZVector(0.,0.,1.),
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToCone(*this, openingAngle, conePoint, coneDirection, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToSphere(double radius,
        const ROOT::Math::XYZVector& point = ROOT::Math::XYZVector(0.,0.,0.),
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToSphere(*this, radius, point, stopAtBoundary, calcJacobianNoise);}
  double extrapolateBy(double step,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateBy(*this, step, stopAtBoundary, calcJacobianNoise);}
  template<unsigned int dimMeas>
  double extrapolateToMeasurement(const AbsMeasurement<dimMeas>* measurement,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToMeasurement(*this, measurement, stopAtBoundary, calcJacobianNoise);}


  ROOT::Math::XYZVector getPos() const {return rep_->getPos(*this);}
  ROOT::Math::XYZVector getMom() const {return rep_->getMom(*this);}
  ROOT::Math::XYZVector getDir() const {return rep_->getDir(*this);}
  void getPosMom(ROOT::Math::XYZVector& pos, ROOT::Math::XYZVector& mom) const {rep_->getPosMom(*this, pos, mom);}
  void getPosDir(ROOT::Math::XYZVector& pos, ROOT::Math::XYZVector& dir) const {rep_->getPosDir(*this, pos, dir);}
  SVector6 get6DState() const {return rep_->get6DState(*this);}
  double getMomMag() const {return rep_->getMomMag(*this);}
  int getPDG() const {return rep_->getPDG();}
  double getCharge() const {return rep_->getCharge(*this);}
  double getQop() const {return rep_->getQop(*this);}
  double getMass() const {return rep_->getMass(*this);}
  double getTime() const {return rep_->getTime(*this);}

  void setPosMom(const ROOT::Math::XYZVector& pos, const ROOT::Math::XYZVector& mom) {rep_->setPosMom(*this, pos, mom);}
  void setPosMom(const SVector6& state6) {rep_->setPosMom(*this, state6);}
  void setChargeSign(double charge) {rep_->setChargeSign(*this, charge);}
  void setQop(double qop) {rep_->setQop(*this, qop);}
  void setTime(double time) {rep_->setTime(*this, time);}


  virtual void Print(Option_t* option = "") const;

 protected:

  SVectorState state_; // state vector
  SVectorAux auxInfo_; // auxiliary information (e.g. charge, flight direction etc.)
  SharedPlanePtr sharedPlane_; //! Shared ownership.  '!' in order to silence ROOT, custom streamer writes and reads this.

 private:

  /** Pointer to TrackRep with respect to which StateOnPlane is defined
   */
  const AbsTrackRep* rep_; //! No ownership

 public:
  ClassDef(StateOnPlane,2)
  // Version history:
  // ver 2: no longer derives from TObject (the TObject parts were not 
  //        streamed, so no compatibility issues arise.)
};


template<unsigned int dim, unsigned int dimAux>
inline StateOnPlane<dim, dimAux>::StateOnPlane(const AbsTrackRep* rep) :
  state_(0), auxInfo_(0), sharedPlane_(), rep_(rep)
{
  if (rep != nullptr) {
    state_.ResizeTo(rep->getDim());
  }
}

template<unsigned int dim, unsigned int dimAux>
inline StateOnPlane<dim, dimAux>::StateOnPlane(const SVectorState& state,
                                               const SharedPlanePtr& plane,
                                               const AbsTrackRep* rep) :
  state_(state), auxInfo_(0), sharedPlane_(plane), rep_(rep)
{
  assert(rep != nullptr);
  assert(sharedPlane_.get() != nullptr);
}

template<unsigned int dim, unsigned int dimAux>
inline StateOnPlane<dim, dimAux>::StateOnPlane(const SVectorState& state,
                                               const SharedPlanePtr& plane,
                                               const AbsTrackRep* rep, 
                                               const SVectorAux& auxInfo) :
  state_(state), auxInfo_(auxInfo), sharedPlane_(plane), rep_(rep)
{
  assert(rep != nullptr);
  assert(sharedPlane_.get() != nullptr);
}

template<unsigned int dim, unsigned int dimAux>
inline StateOnPlane<dim, dimAux>& StateOnPlane<dim, dimAux>::operator=(StateOnPlane other) {
  swap(other);
  return *this;
}

template<unsigned int dim, unsigned int dimAux>
inline void StateOnPlane<dim, dimAux>::swap(StateOnPlane& other) {
  this->state_.ResizeTo(other.state_);
  std::swap(this->state_, other.state_);
  this->auxInfo_.ResizeTo(other.auxInfo_);
  std::swap(this->auxInfo_, other.auxInfo_);
  this->sharedPlane_.swap(other.sharedPlane_);
  std::swap(this->rep_, other.rep_);
}

} /* End of namespace genfit */
/** @} */

#endif // genfit_StateOnPlane_h
