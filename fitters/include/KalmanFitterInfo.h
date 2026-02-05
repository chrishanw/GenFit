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

#ifndef genfit_KalmanFitterInfo_h
#define genfit_KalmanFitterInfo_h

#include "AbsFitterInfo.h"
#include "KalmanFittedStateOnPlane.h"
#include "MeasuredStateOnPlane.fwd.h"
#include "MeasurementOnPlane.h"
#include "ReferenceStateOnPlane.h"
#include "StateOnPlane.fwd.h"

#include <vector>

#include <memory>


namespace genfit {


/**
 *  @brief Collects information needed and produced by a AbsKalmanFitter implementations and is specific to one AbsTrackRep of the Track.
 */
template<unsigned int dim, unsigned int dimAux>
class KalmanFitterInfo : public AbsFitterInfo {

 public:

  KalmanFitterInfo();
  KalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);
  virtual ~KalmanFitterInfo();

  virtual KalmanFitterInfo* clone() const override;

  ReferenceStateOnPlane<dim, dimAux>* getReferenceState() const {return referenceState_.get();}
  MeasuredStateOnPlane<dim, dimAux>* getForwardPrediction() const {return forwardPrediction_.get();}
  MeasuredStateOnPlane<dim, dimAux>* getBackwardPrediction() const {return backwardPrediction_.get();}
  MeasuredStateOnPlane<dim, dimAux>* getPrediction(int direction) const {if (direction >=0) return forwardPrediction_.get(); return backwardPrediction_.get();}
  KalmanFittedStateOnPlane<dim, dimAux>* getForwardUpdate() const {return forwardUpdate_.get();}
  KalmanFittedStateOnPlane<dim, dimAux>* getBackwardUpdate() const {return backwardUpdate_.get();}
  KalmanFittedStateOnPlane<dim, dimAux>* getUpdate(int direction) const {if (direction >=0) return forwardUpdate_.get(); return backwardUpdate_.get();}
  const std::vector< genfit::MeasurementOnPlane<dim, dimAux>* >& getMeasurementsOnPlane() const {return measurementsOnPlane_;}
  MeasurementOnPlane<dim, dimAux>* getMeasurementOnPlane(int i = 0) const {if (i<0) i += measurementsOnPlane_.size(); return measurementsOnPlane_.at(i);}
  //! Get weighted mean of all measurements.
  //! @param ignoreWeights If set, the weights of the individual measurements will be ignored (they will be treated as if they all had weight 1)
  MeasurementOnPlane<dim, dimAux> getAvgWeightedMeasurementOnPlane(bool ignoreWeights = false) const;
  //! Get measurements which is closest to state.
  MeasurementOnPlane<dim, dimAux>* getClosestMeasurementOnPlane(const StateOnPlane<dim, dimAux>*) const;
  unsigned int getNumMeasurements() const {return measurementsOnPlane_.size();}
  //! Get weights of measurements.
  std::vector<double> getWeights() const;
  //! Are the weights fixed?
  bool areWeightsFixed() const {return fixWeights_;}
  //! Get unbiased or biased (default) smoothed state.
  const MeasuredStateOnPlane<dim, dimAux>& getFittedState(bool biased = true) const;
  //! Get unbiased (default) or biased residual from ith measurement.
  MeasurementOnPlane<dim, dimAux> getResidual(unsigned int iMeasurement = 0, bool biased = false, bool onlyMeasurementErrors = true) const; // calculate residual, track and measurement errors are added if onlyMeasurementErrors is false
  double getSmoothedChi2(unsigned int iMeasurement = 0) const;

  bool hasMeasurements() const override {return getNumMeasurements() > 0;}
  bool hasReferenceState() const override {return (referenceState_.get() != nullptr);}
  bool hasForwardPrediction() const override {return (forwardPrediction_.get()  != nullptr);}
  bool hasBackwardPrediction() const override {return (backwardPrediction_.get() != nullptr);}
  bool hasForwardUpdate() const override {return (forwardUpdate_.get() != nullptr);}
  bool hasBackwardUpdate() const override {return (backwardUpdate_.get() != nullptr);}
  bool hasUpdate(int direction) const override {if (direction < 0) return hasBackwardUpdate(); return hasForwardUpdate();}
  bool hasPredictionsAndUpdates() const {return (hasForwardPrediction() && hasBackwardPrediction() && hasForwardUpdate() && hasBackwardUpdate());}

  void setReferenceState(ReferenceStateOnPlane<dim, dimAux>* referenceState);
  void setForwardPrediction(MeasuredStateOnPlane<dim, dimAux>* forwardPrediction);
  void setBackwardPrediction(MeasuredStateOnPlane<dim, dimAux>* backwardPrediction);
  void setPrediction(MeasuredStateOnPlane<dim, dimAux>* prediction, int direction)  {if (direction >=0) setForwardPrediction(prediction); else setBackwardPrediction(prediction);}
  void setForwardUpdate(KalmanFittedStateOnPlane<dim, dimAux>* forwardUpdate);
  void setBackwardUpdate(KalmanFittedStateOnPlane<dim, dimAux>* backwardUpdate);
  void setUpdate(KalmanFittedStateOnPlane<dim, dimAux>* update, int direction)  {if (direction >=0) setForwardUpdate(update); else setBackwardUpdate(update);}
  void setMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane<dim, dimAux>* >& measurementsOnPlane);
  void addMeasurementOnPlane(MeasurementOnPlane<dim, dimAux>* measurementOnPlane);
  void addMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane<dim, dimAux>* >& measurementsOnPlane);
  //! Set weights of measurements.
  void setWeights(const std::vector<double>&);
  void fixWeights(bool arg = true) {fixWeights_ = arg;}
  void setRep(const AbsTrackRep* rep) override;

  void deleteForwardInfo() override;
  void deleteBackwardInfo() override;
  void deletePredictions();
  void deleteReferenceInfo() {setReferenceState(nullptr);}
  void deleteMeasurementInfo() override;

  virtual void Print(const Option_t* = "") const override;

  virtual bool checkConsistency(const genfit::PruneFlags* = nullptr) const override;

 private:

  //! Reference state. Used by KalmanFitterRefTrack.
  std::unique_ptr<ReferenceStateOnPlane<dim, dimAux>> referenceState_; // Ownership
  std::unique_ptr<MeasuredStateOnPlane<dim, dimAux>> forwardPrediction_; // Ownership
  std::unique_ptr<KalmanFittedStateOnPlane<dim, dimAux>> forwardUpdate_; // Ownership
  std::unique_ptr<MeasuredStateOnPlane<dim, dimAux>> backwardPrediction_; // Ownership
  std::unique_ptr<KalmanFittedStateOnPlane<dim, dimAux>> backwardUpdate_; // Ownership
  mutable std::unique_ptr<MeasuredStateOnPlane<dim, dimAux>> fittedStateUnbiased_; //!  cache
  mutable std::unique_ptr<MeasuredStateOnPlane<dim, dimAux>> fittedStateBiased_; //!  cache

 //> TODO ! ptr implement: to the special ownership version
  /* class owned_pointer_vector : private std::vector<MeasuredStateOnPlane*> {
   public: 
    ~owned_pointer_vector() { for (size_t i = 0; i < this->size(); ++i)
                         delete this[i]; }
    size_t size() const { return this->size(); };
    void push_back(MeasuredStateOnPlane* measuredState) { this->push_back(measuredState); };
    const  MeasuredStateOnPlane* at(size_t i)  const { return this->at(i); }; 
	//owned_pointer_vector::iterator erase(owned_pointer_vector::iterator position) ;
	//owned_pointer_vector::iterator erase(owned_pointer_vector::iterator first, owned_pointer_vector::iterator last);
};
	*/

  // template<unsigned int dim, unsigned int dimAux>
  // std::vector<MeasurementOnPlane<dim, dimAux>*> measurementsOnPlane_; // Ownership
  std::vector<MeasurementOnPlane<dim, dimAux>*> measurementsOnPlane_; // Ownership
  bool fixWeights_; // weights should not be altered by fitters anymore

 public:

  ClassDefOverride(KalmanFitterInfo,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_KalmanFitterInfo_h
