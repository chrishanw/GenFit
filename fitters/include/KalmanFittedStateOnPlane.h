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

#ifndef genfit_KalmanFittedStateOnPlane_h
#define genfit_KalmanFittedStateOnPlane_h

#include "MeasuredStateOnPlane.h"


namespace genfit {


/**
 *  @brief #MeasuredStateOnPlane with additional info produced by a Kalman filter or DAF.
 */
template<unsigned int dim, unsigned int dimAux>
class KalmanFittedStateOnPlane : public MeasuredStateOnPlane<dim, dimAux> {

  using Super = MeasuredStateOnPlane<dim, dimAux>;
  using SVectorState = ROOT::Math::SVector<double, dim>;
  using SMatrixCov = ROOT::Math::SMatrix<double, dim, dim, ROOT::Math::MatRepSym<double, dim> >;
  using SVectorAux = ROOT::Math::SVector<double, dimAux>;

 public:

  KalmanFittedStateOnPlane();
  KalmanFittedStateOnPlane(const KalmanFittedStateOnPlane&) = default;
  KalmanFittedStateOnPlane(const SVectorState& state,
                           const SMatrixCov& cov,
                           const SharedPlanePtr& plane,
                           const AbsTrackRep* rep,
                           double chiSquareIncrement,
                           double ndf);
  KalmanFittedStateOnPlane(const SVectorState& state,
                           const SMatrixCov& cov,
                           const SharedPlanePtr& plane,
                           const AbsTrackRep* rep,
                           const SVectorAux& auxInfo,
                           double chiSquareIncrement,
                           double ndf);
  KalmanFittedStateOnPlane(const Super& state, double chiSquareIncrement, double ndf);

  KalmanFittedStateOnPlane& operator=(KalmanFittedStateOnPlane other);
  void swap(KalmanFittedStateOnPlane& other); // nothrow

  virtual ~KalmanFittedStateOnPlane() {}

  double getChiSquareIncrement() const {return chiSquareIncrement_;}
  double getNdf() const {return ndf_;}

  void setChiSquareIncrement(double chiSquareIncrement) {chiSquareIncrement_ = chiSquareIncrement;}
  void setNdf(double ndf) {ndf_ = ndf;}


 protected:

  double chiSquareIncrement_;
  
  //! Degrees of freedom. Needs to be a double because of DAF.
  double ndf_;


 public:

  ClassDef(KalmanFittedStateOnPlane,1)

};


template<unsigned int dim, unsigned int dimAux>
inline KalmanFittedStateOnPlane<dim, dimAux>::KalmanFittedStateOnPlane() :
  Super(), chiSquareIncrement_(0), ndf_(0)
{
  ;
}

template<unsigned int dim, unsigned int dimAux>
inline KalmanFittedStateOnPlane<dim, dimAux>::KalmanFittedStateOnPlane(const SVectorState& state,
                                                          const SMatrixCov& cov,
                                                          const SharedPlanePtr& plane,
                                                          const AbsTrackRep* rep,
                                                          double chiSquareIncrement,
                                                          double ndf) :
  Super(state, cov, plane, rep), chiSquareIncrement_(chiSquareIncrement), ndf_(ndf)
{
  ;
}

template<unsigned int dim, unsigned int dimAux>
inline KalmanFittedStateOnPlane<dim, dimAux>::KalmanFittedStateOnPlane(const SVectorState& state,
                                                          const SMatrixCov& cov,
                                                          const SharedPlanePtr& plane,
                                                          const AbsTrackRep* rep,
                                                          const SVectorAux& auxInfo,
                                                          double chiSquareIncrement,
                                                          double ndf) :
  Super(state, cov, plane, rep, auxInfo), chiSquareIncrement_(chiSquareIncrement), ndf_(ndf)
{
  ;
}

template<unsigned int dim, unsigned int dimAux>
inline KalmanFittedStateOnPlane<dim, dimAux>::KalmanFittedStateOnPlane(const Super& state, double chiSquareIncrement, double ndf) :
  Super(state), chiSquareIncrement_(chiSquareIncrement), ndf_(ndf)
{
  ;
}

template<unsigned int dim, unsigned int dimAux>
inline KalmanFittedStateOnPlane<dim, dimAux>& KalmanFittedStateOnPlane<dim, dimAux>::operator=(KalmanFittedStateOnPlane<dim, dimAux> other) {
  swap(other);
  return *this;
}

template<unsigned int dim, unsigned int dimAux>
inline void KalmanFittedStateOnPlane<dim, dimAux>::swap(KalmanFittedStateOnPlane<dim, dimAux>& other) {
  Super::swap(other);
  std::swap(this->chiSquareIncrement_, other.chiSquareIncrement_);
  std::swap(this->ndf_, other.ndf_);
}

} /* End of namespace genfit */
/** @} */

#endif // genfit_KalmanFittedStateOnPlane_h
