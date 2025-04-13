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

#ifndef genfit_AbsTrackRep_h
#define genfit_AbsTrackRep_h

#include "SharedPlanePtr.h"
//#include "MaterialInfo.h"
#include "Material.h"
#include <StateOnPlane.fwd.h>
#include <MeasuredStateOnPlane.fwd.h>
#include <SMatrixTypeDefs.h>
#include <AbsMeasurement.fwd.h>

#include <Math/Vector3D.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>


namespace genfit {

/**
 * @brief Simple struct containing MaterialProperties and stepsize in the material.
 */
struct MatStep {
  Material material_;
  double stepSize_;

  MatStep() {
    stepSize_ = 0;
  }

};

// class StateOnPlane;
// class MeasuredStateOnPlane;
// class AbsMeasurement;

/**
 * @brief Abstract base class for a track representation
 *
 *  Provides functionality to extrapolate a StateOnPlane to another DetPlane,
 *  to the POCA to a line or a point, or a cylinder or sphere.
 *  Defines a set of parameters describing the track.
 *  StateOnPlane objects are always defined with a track parameterization of a specific AbsTrackRep.
 *  The AbsTrackRep provides functionality to translate from the internal representation of a state
 *  into cartesian position and momentum (and covariance) and vice versa.
 */
class AbsTrackRep {

 public:

  AbsTrackRep();
  AbsTrackRep(int pdgCode, char propDir = 0);

  virtual ~AbsTrackRep() {;}

  //! Clone the trackRep.
  virtual AbsTrackRep* clone() const = 0;

  /**
   * @brief Extrapolates the state to plane, and returns the extrapolation length
   *        and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateToPlane(
      StateOnPlane<dim, dimAux>& state,
      const genfit::SharedPlanePtr& plane,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  /**
   * @brief Extrapolates the state to the POCA to a line, and returns the extrapolation length
   *        and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateToLine(StateOnPlane<dim, dimAux>& state,
      const ROOT::Math::XYZVector& linePoint,
      const ROOT::Math::XYZVector& lineDirection,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  /**
   * @brief Resembles the interface of GFAbsTrackRep in old versions of genfit
   *
   * This interface to extrapolateToLine is intended to resemble the
   * interface of GFAbsTrackRep in old versions of genfit and is
   * implemented by default via the preceding function.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateToLine(StateOnPlane<dim, dimAux>& state,
      const ROOT::Math::XYZVector& point1,
      const ROOT::Math::XYZVector& point2,
      ROOT::Math::XYZVector& poca,
      ROOT::Math::XYZVector& dirInPoca,
      ROOT::Math::XYZVector& poca_onwire,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const {
    ROOT::Math::XYZVector wireDir(point2 - point1);
    wireDir = wireDir.Unit();
    double retval = this->extrapolateToLine(state, point1, wireDir, stopAtBoundary, calcJacobianNoise);
    poca = this->getPos(state);
    dirInPoca = this->getMom(state);
    dirInPoca = dirInPoca.Unit();

    poca_onwire = point1 + wireDir*((poca - point1).Dot(wireDir));
    
    return retval;
  }

  /**
   * @brief Extrapolates the state to the POCA to a point, and returns the extrapolation length
   *        and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateToPoint(StateOnPlane<dim, dimAux>& state,
      const ROOT::Math::XYZVector& point,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  /**
   * @brief Extrapolates the state to the POCA to a point in the metric of G, and returns the extrapolation length
   *        and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateToPoint(StateOnPlane<dim, dimAux>& state,
      const ROOT::Math::XYZVector& point,
      const TMatrixDSym& G, // weight matrix (metric)
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  /**
   * @brief Extrapolates the state to the cylinder surface, and returns the extrapolation length
   *       and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateToCylinder(StateOnPlane<dim, dimAux>& state,
      double radius,
      const ROOT::Math::XYZVector& linePoint = ROOT::Math::XYZVector(0.,0.,0.),
      const ROOT::Math::XYZVector& lineDirection = ROOT::Math::XYZVector(0.,0.,1.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  /**
   * @brief Extrapolates the state to the cone surface, and returns the extrapolation length
   *       and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateToCone(StateOnPlane<dim, dimAux>& state,
      double radius,
      const ROOT::Math::XYZVector& linePoint = ROOT::Math::XYZVector(0.,0.,0.),
      const ROOT::Math::XYZVector& lineDirection = ROOT::Math::XYZVector(0.,0.,1.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  /**
   * @brief Extrapolates the state to the sphere surface, and returns the extrapolation length
   *       and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateToSphere(StateOnPlane<dim, dimAux>& state,
      double radius,
      const ROOT::Math::XYZVector& point = ROOT::Math::XYZVector(0.,0.,0.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  /**
   * @brief Extrapolates the state by step (cm) and returns the extrapolation length
   *       and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  template<unsigned int dim, unsigned int dimAux>
  double extrapolateBy(StateOnPlane<dim, dimAux>& state,
      double step,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  //! extrapolate to an AbsMeasurement
  template<unsigned int dim, unsigned int dimAux, unsigned int dimMeas>
  double extrapolateToMeasurement(StateOnPlane<dim, dimAux>& state,
      const AbsMeasurement<dimMeas>* measurement,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  //! Get the dimension of the state vector used by the track representation.
  virtual unsigned int getDim() const = 0;

  //! Get the cartesian position of a state.
  template<unsigned int dim, unsigned int dimAux>
  ROOT::Math::XYZVector getPos(const StateOnPlane<dim, dimAux>& state) const;

  //! Get the cartesian momentum vector of a state.
  template<unsigned int dim, unsigned int dimAux>
  ROOT::Math::XYZVector getMom(const StateOnPlane<dim, dimAux>& state) const;

  //! Get the direction vector of a state.
  template<unsigned int dim, unsigned int dimAux>
  ROOT::Math::XYZVector getDir(const StateOnPlane<dim, dimAux>& state) const {return getMom(state).Unit();}

  //! Get cartesian position and momentum vector of a state.
  template<unsigned int dim, unsigned int dimAux>
  void getPosMom(const StateOnPlane<dim, dimAux>& state, ROOT::Math::XYZVector& pos, ROOT::Math::XYZVector& mom) const;

  //! Get cartesian position and direction vector of a state.
  template<unsigned int dim, unsigned int dimAux>
  void getPosDir(const StateOnPlane<dim, dimAux>& state, ROOT::Math::XYZVector& pos, ROOT::Math::XYZVector& dir) const {getPosMom(state, pos, dir); dir *= 1. / dir.R();}

  //! Get the 6D state vector (x, y, z, p_x, p_y, p_z).
  template<unsigned int dim, unsigned int dimAux>
  SVector6 get6DState(const StateOnPlane<dim, dimAux>& state) const;

  //! Get the 6D covariance.
  template<unsigned int dim, unsigned int dimAux>
  SMatrixSym6 get6DCov(const MeasuredStateOnPlane<dim, dimAux>& state) const;

  //! Translates MeasuredStateOnPlane into 3D position, momentum and 6x6 covariance.
  template<unsigned int dim, unsigned int dimAux>
  void getPosMomCov(const MeasuredStateOnPlane<dim, dimAux>& state, ROOT::Math::XYZVector& pos, ROOT::Math::XYZVector& mom, SMatrixSym6& cov) const;

  //! Translates MeasuredStateOnPlane into 6D state vector (x, y, z, p_x, p_y, p_z) and 6x6 covariance.
  template<unsigned int dim, unsigned int dimAux>
  void get6DStateCov(const MeasuredStateOnPlane<dim, dimAux>& state, SVector6& stateVec, SMatrixSym6& cov) const;

  //! get the magnitude of the momentum in GeV.
  template<unsigned int dim, unsigned int dimAux>
  double getMomMag(const StateOnPlane<dim, dimAux>& state) const;
  //! get the variance of the absolute value of the momentum .
  template<unsigned int dim, unsigned int dimAux>
  double getMomVar(const MeasuredStateOnPlane<dim, dimAux>& state) const;

  //! Get the pdg code.
  int getPDG() const {return pdgCode_;}

  //! Get the charge of the particle of the pdg code
  double getPDGCharge() const;

  /**
   * @brief Get the (fitted) charge of a state.
   * This is not always equal the pdg charge (e.g. if the charge sign was flipped during the fit).
   */
  template<unsigned int dim, unsigned int dimAux>
  double getCharge(const StateOnPlane<dim, dimAux>& state) const;
  //! Get charge over momentum.
  template<unsigned int dim, unsigned int dimAux>
  double getQop(const StateOnPlane<dim, dimAux>& state) const;
  //! Get tha particle mass in GeV/c^2
  template<unsigned int dim, unsigned int dimAux>
  double getMass(const StateOnPlane<dim, dimAux>& state) const;

  //! Get propagation direction. (-1, 0, 1) -> (backward, auto, forward).
  char getPropDir() const {return propDir_;}

  //! Get the jacobian and noise matrix of the last extrapolation.
  virtual void getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const = 0;

  //! Get the jacobian and noise matrix of the last extrapolation if it would have been done in opposite direction.
  virtual void getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const = 0;

  //! Get stepsizes and material properties of crossed materials of the last extrapolation.
  virtual std::vector<genfit::MatStep> getSteps() const = 0;

  //! Get the accumulated X/X0 (path / radiation length) of the material crossed in the last extrapolation.
  virtual double getRadiationLenght() const = 0;

  //! Get the time corresponding to the StateOnPlane.  Extrapolation
  // should keep this up to date with the time of flight.
  template<unsigned int dim, unsigned int dimAux>
  double getTime(const StateOnPlane<dim, dimAux>&) const;

  /**
   * @brief Calculate Jacobian of transportation numerically.
   * Slow but accurate. Can be used to validate (semi)analytic calculations.
   */
  template<unsigned int dim, unsigned int dimAux>
  void calcJacobianNumerically(const genfit::StateOnPlane<dim, dimAux>& origState,
                                   const genfit::SharedPlanePtr destPlane,
                                   TMatrixD& jacobian) const;

  //! try to multiply pdg code with -1. (Switch from particle to anti-particle and vice versa).
  bool switchPDGSign();

  //! Set position and momentum of state.
  template<unsigned int dim, unsigned int dimAux>
  void setPosMom(StateOnPlane<dim, dimAux>& state, const ROOT::Math::XYZVector& pos, const ROOT::Math::XYZVector& mom) const;
  //! Set position and momentum of state.
  template<unsigned int dim, unsigned int dimAux>
  void setPosMom(StateOnPlane<dim, dimAux>& state, const SVector6& state6) const;
  //! Set position and momentum and error of state.
  template<unsigned int dim, unsigned int dimAux>
  void setPosMomErr(MeasuredStateOnPlane<dim, dimAux>& state, const ROOT::Math::XYZVector& pos, const ROOT::Math::XYZVector& mom, const ROOT::Math::XYZVector& posErr, const ROOT::Math::XYZVector& momErr) const;
  //! Set position, momentum and covariance of state.
  template<unsigned int dim, unsigned int dimAux>
  void setPosMomCov(MeasuredStateOnPlane<dim, dimAux>& state, const ROOT::Math::XYZVector& pos, const ROOT::Math::XYZVector& mom, const SMatrixSym6& cov6x6) const;
  //! Set position, momentum and covariance of state.
  template<unsigned int dim, unsigned int dimAux>
  void setPosMomCov(MeasuredStateOnPlane<dim, dimAux>& state, const SVector6& state6, const SMatrixSym6& cov6x6) const;

  //! Set the sign of the charge according to charge.
  template<unsigned int dim, unsigned int dimAux>
  void setChargeSign(StateOnPlane<dim, dimAux>& state, double charge) const;
  //! Set charge/momentum.
  template<unsigned int dim, unsigned int dimAux>
  void setQop(StateOnPlane<dim, dimAux>& state, double qop) const;
  //! Set time at which the state was defined
  template<unsigned int dim, unsigned int dimAux>
  void setTime(StateOnPlane<dim, dimAux>& state, double time) const;

  //! Set propagation direction. (-1, 0, 1) -> (backward, auto, forward).
  void setPropDir(int dir) {
    if (dir>0) propDir_ = 1;
    else if (dir<0) propDir_ = -1;
    else propDir_ = 0;
  };

  //! Switch propagation direction. Has no effect if propDir_ is set to 0.
  void switchPropDir(){propDir_ = -1*propDir_;}

  //! check if other is of same type (e.g. RKTrackRep).
  virtual bool isSameType(const AbsTrackRep* other) = 0;

  //! check if other is of same type (e.g. RKTrackRep) and has same pdg code.
  virtual bool isSame(const AbsTrackRep* other) = 0;

  virtual void setDebugLvl(unsigned int lvl = 1) {debugLvl_ = lvl;}

  virtual void Print(const Option_t* = "") const;

 protected:

  //! protect from calling copy c'tor from outside the class. Use #clone() if you want a copy!
  AbsTrackRep(const AbsTrackRep&);
  //! protect from calling assignment operator from outside the class. Use #clone() instead!
  AbsTrackRep& operator=(const AbsTrackRep&);


  //! Particle code
  int pdgCode_;
  //! propagation direction (-1, 0, 1) -> (backward, auto, forward)
  char propDir_;

  unsigned int debugLvl_;

 public:
  ClassDef(AbsTrackRep,2)
  // Version history:
  // ver 2: no longer derives from TObject

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsTrackRep_h
