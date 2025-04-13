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

#ifndef genfit_MeasuredStateOnPlane_h
#define genfit_MeasuredStateOnPlane_h

#include "StateOnPlane.h"
#include "AbsTrackRep.h"
#include <SMatrixTypeDefs.h>

#include <TMatrixDSym.h>

#include <cassert>

namespace genfit {

/**
 *  @brief #StateOnPlane with additional covariance matrix.
 */
template<unsigned int dim, unsigned int dimAux = 0>
class MeasuredStateOnPlane : public StateOnPlane<dim, dimAux> {

  using Super = StateOnPlane<dim, dimAux>;
  using SVectorState = ROOT::Math::SVector<double, dim>;
  using SMatrixCov = ROOT::Math::SMatrix<double, dim, dim, ROOT::Math::MatRepSym<double, dim> >;
  using SVectorAux = ROOT::Math::SVector<double, dimAux>;

 public:

  MeasuredStateOnPlane(const AbsTrackRep* rep = nullptr);
  MeasuredStateOnPlane(const SVectorState& state, const SMatrixCov& cov, const genfit::SharedPlanePtr& plane, const AbsTrackRep* rep);
  MeasuredStateOnPlane(const SVectorState& state, const SMatrixCov& cov, const genfit::SharedPlanePtr& plane, const AbsTrackRep* rep, const SVectorAux& auxInfo);
  MeasuredStateOnPlane(const MeasuredStateOnPlane& o);
  MeasuredStateOnPlane(const Super& state, const SMatrixCov& cov);

  MeasuredStateOnPlane& operator=(MeasuredStateOnPlane other);
  void swap(MeasuredStateOnPlane& other); // nothrow

  virtual ~MeasuredStateOnPlane() {}
  virtual MeasuredStateOnPlane* clone() const override {return new MeasuredStateOnPlane(*this);}


  const SMatrixCov& getCov() const {return cov_;}
  SMatrixCov& getCov() {return cov_;}

  //! Blow up covariance matrix with blowUpFac. Per default, off diagonals are reset to 0 and the maximum values are limited to maxVal.
  void blowUpCov(double blowUpFac, bool resetOffDiagonals = true, double maxVal = -1.);

  void setStateCov(const SVectorState& state, const SMatrixCov& cov) {setState(state); setCov(cov);}
  void setStateCovPlane(const SVectorState& state, const SMatrixCov& cov, const SharedPlanePtr& plane) {setStatePlane(state, plane); setCov(cov);}
  void setCov(const SMatrixCov& cov) {if(cov_.GetNrows() == 0) cov_.ResizeTo(cov); cov_ = cov;}

  // Shortcuts to TrackRep functions
  SMatrixSym6 get6DCov() const {return Super::getRep()->get6DCov(*this);};
  void getPosMomCov(ROOT::Math::XYZVector& pos, ROOT::Math::XYZVector& mom, SMatrixSym6& cov) const {Super::getRep()->getPosMomCov(*this, pos, mom, cov);}
  void get6DStateCov(SVector6& stateVec, SMatrixSym6& cov) const {Super::getRep()->get6DStateCov(*this, stateVec, cov);}
  double getMomVar() const {return Super::getRep()->getMomVar(*this);}

  void setPosMomErr(const ROOT::Math::XYZVector& pos, const ROOT::Math::XYZVector& mom, const ROOT::Math::XYZVector& posErr, const ROOT::Math::XYZVector& momErr) {Super::getRep()->setPosMomErr(*this, pos, mom, posErr, momErr);}
  void setPosMomCov(const ROOT::Math::XYZVector& pos, const ROOT::Math::XYZVector& mom, const SMatrixSym6& cov6x6) {Super::getRep()->setPosMomCov(*this, pos, mom, cov6x6);}
  void setPosMomCov(const SVector6& state6, const SMatrixSym6& cov6x6) {Super::getRep()->setPosMomCov(*this, state6, cov6x6);}


  virtual void Print(Option_t* option = "") const override;

 protected:

  SMatrixCov cov_;

 public:
  ClassDefOverride(MeasuredStateOnPlane,1)

};


/**
 * @brief Calculate weighted average between two MeasuredStateOnPlanes
 */
template<unsigned int dim, unsigned int dimAux>
MeasuredStateOnPlane<dim, dimAux> calcAverageState(const MeasuredStateOnPlane<dim, dimAux>& forwardState, const MeasuredStateOnPlane<dim, dimAux>& backwardState);


template<unsigned int dim, unsigned int dimAux>
inline void MeasuredStateOnPlane<dim, dimAux>::swap(MeasuredStateOnPlane<dim, dimAux>& other) {
  Super::swap(other);
  this->cov_.ResizeTo(other.cov_);
  std::swap(this->cov_, other.cov_);
}

template<unsigned int dim, unsigned int dimAux>
inline MeasuredStateOnPlane<dim, dimAux>::MeasuredStateOnPlane(const AbsTrackRep* rep) :
  Super(rep), cov_(0,0)
{
  if (rep != nullptr) {
    cov_.ResizeTo(rep->getDim(), rep->getDim());
  }
}

template<unsigned int dim, unsigned int dimAux>
inline MeasuredStateOnPlane<dim, dimAux>::MeasuredStateOnPlane(const SVectorState& state,
                                                  const SMatrixCov& cov,
                                                  const SharedPlanePtr& plane,
                                                  const AbsTrackRep* rep) :
  Super(state, plane, rep), cov_(cov)
{
  assert(rep != nullptr);
  //assert(cov_.GetNcols() == (signed)rep->getDim());
}

template<unsigned int dim, unsigned int dimAux>
inline MeasuredStateOnPlane<dim, dimAux>::MeasuredStateOnPlane(const SVectorState& state,
                                                  const SMatrixCov& cov,
                                                  const SharedPlanePtr& plane,
                                                  const AbsTrackRep* rep,
                                                  const SVectorAux& auxInfo) :
  Super(state, plane, rep, auxInfo), cov_(cov)
{
  assert(rep != nullptr);
  //assert(cov_.GetNcols() == (signed)rep->getDim());
}

template<unsigned int dim, unsigned int dimAux>
inline MeasuredStateOnPlane<dim, dimAux>::MeasuredStateOnPlane(const MeasuredStateOnPlane<dim, dimAux>& o) :
  Super(o), cov_(o.cov_)
{
}

template<unsigned int dim, unsigned int dimAux>
inline MeasuredStateOnPlane<dim, dimAux>::MeasuredStateOnPlane(const Super& state,
                                                  const SMatrixCov& cov) :
  Super(state), cov_(cov)
{
  //assert(cov_.GetNcols() == (signed)getRep()->getDim());
}

template<unsigned int dim, unsigned int dimAux>
inline MeasuredStateOnPlane<dim, dimAux>& MeasuredStateOnPlane<dim, dimAux>::operator=(MeasuredStateOnPlane<dim, dimAux> other) {
  swap(other);
  return *this;
}


} /* End of namespace genfit */
/** @} */

#endif // genfit_MeasuredStateOnPlane_h
