/** @addtogroup genfit
 * @{
 */

#pragma once

#include "AbsFitterInfo.h"
#include "MeasuredStateOnPlane.h"
#include "MeasurementOnPlane.h"

#include <memory>
#include <vector>

namespace genfit {

/** A Gaussian component state and its normalized mixture weight. */
struct GSFComponent {
  GSFComponent() : state(), weight(0.) {}
  GSFComponent(const MeasuredStateOnPlane& componentState, double componentWeight)
      : state(componentState), weight(componentWeight) {}

  MeasuredStateOnPlane state;
  double weight;

  ClassDefNV(GSFComponent, 1)
};

/**
 * @brief Per-point fit information for the Gaussian-sum fitter.
 *
 * Stores moment-matched predictions and updates for standard GenFit fitted-
 * state access, as well as the full predicted and updated component mixtures.
 */
class GSFFitterInfo : public AbsFitterInfo {

 public:
  GSFFitterInfo();
  GSFFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);
  ~GSFFitterInfo() override = default;

  GSFFitterInfo* clone() const override;

  MeasuredStateOnPlane* getForwardPrediction() const { return forwardPrediction_.get(); }
  MeasuredStateOnPlane* getBackwardPrediction() const { return backwardPrediction_.get(); }
  MeasuredStateOnPlane* getPrediction(int direction) const { return direction >= 0 ? getForwardPrediction() : getBackwardPrediction(); }
  MeasuredStateOnPlane* getForwardUpdate() const { return forwardUpdate_.get(); }
  MeasuredStateOnPlane* getBackwardUpdate() const { return backwardUpdate_.get(); }
  MeasuredStateOnPlane* getUpdate(int direction) const { return direction >= 0 ? getForwardUpdate() : getBackwardUpdate(); }
  double getForwardUpdateChi2Increment() const { return forwardUpdateChi2Increment_; }
  double getBackwardUpdateChi2Increment() const { return backwardUpdateChi2Increment_; }
  double getForwardUpdateNdf() const { return forwardUpdateNdf_; }
  double getBackwardUpdateNdf() const { return backwardUpdateNdf_; }

  const std::vector<MeasurementOnPlane>& getMeasurementsOnPlane() const { return measurementsOnPlane_; }
  MeasurementOnPlane* getMeasurementOnPlane(int i = 0) { if (i < 0) i += measurementsOnPlane_.size(); return &measurementsOnPlane_.at(i); }
  const MeasurementOnPlane* getMeasurementOnPlane(int i = 0) const { if (i < 0) i += measurementsOnPlane_.size(); return &measurementsOnPlane_.at(i); }
  unsigned int getNumMeasurements() const { return measurementsOnPlane_.size(); }
  std::vector<double> getWeights() const override;
  bool areWeightsFixed() const { return fixWeights_; }

  void setForwardPrediction(MeasuredStateOnPlane* prediction);
  void setBackwardPrediction(MeasuredStateOnPlane* prediction);
  void setPrediction(MeasuredStateOnPlane* prediction, int direction);
  void setForwardUpdate(const MeasuredStateOnPlane& update, double chi2Increment = 0., double ndf = 0.);
  void setBackwardUpdate(const MeasuredStateOnPlane& update, double chi2Increment = 0., double ndf = 0.);
  void setUpdate(const MeasuredStateOnPlane& update, int direction, double chi2Increment = 0., double ndf = 0.);
  void addMeasurementOnPlane(const MeasurementOnPlane& measurement);
  void setWeights(const std::vector<double>& weights);
  void fixWeights(bool value = true) { fixWeights_ = value; }

  //! Get predicted (updated=false) or updated (updated=true) components.
  const std::vector<GSFComponent>& getComponents(int direction, bool updated = true) const;
  //! Store predicted (updated=false) or updated (updated=true) components.
  void setComponents(const std::vector<GSFComponent>& components, int direction, bool updated = true);

  void setRep(const AbsTrackRep* rep) override;
  bool hasMeasurements() const override { return !measurementsOnPlane_.empty(); }
  bool hasReferenceState() const override { return false; }

  bool hasForwardPrediction() const override {return (forwardPrediction_.get()  != nullptr);}
  bool hasBackwardPrediction() const override {return (backwardPrediction_.get() != nullptr);}
  bool hasForwardUpdate() const override {return (forwardUpdate_.get() != nullptr);}
  bool hasBackwardUpdate() const override {return (backwardUpdate_.get() != nullptr);}
  void deleteReferenceInfo() override {}
  void deleteForwardInfo() override;
  void deleteBackwardInfo() override;
  void deleteMeasurementInfo() override;
  const MeasuredStateOnPlane& getFittedState(bool biased = true) const override;
  MeasurementOnPlane getResidual(unsigned int iMeasurement = 0, bool biased = false, bool onlyMeasurementErrors = true) const override;
  bool checkConsistency(const PruneFlags* = nullptr) const override;

 private:
  std::unique_ptr<MeasuredStateOnPlane> forwardPrediction_;
  std::unique_ptr<MeasuredStateOnPlane> backwardPrediction_;
  std::unique_ptr<MeasuredStateOnPlane> forwardUpdate_;
  std::unique_ptr<MeasuredStateOnPlane> backwardUpdate_;
  double forwardUpdateChi2Increment_ = 0.;
  double backwardUpdateChi2Increment_ = 0.;
  double forwardUpdateNdf_ = 0.;
  double backwardUpdateNdf_ = 0.;
  bool fixWeights_ = false;
  std::vector<MeasurementOnPlane> measurementsOnPlane_;
  std::vector<GSFComponent> forwardPredictedComponents_;
  std::vector<GSFComponent> forwardUpdatedComponents_;
  std::vector<GSFComponent> backwardPredictedComponents_;
  std::vector<GSFComponent> backwardUpdatedComponents_;
  mutable std::unique_ptr<MeasuredStateOnPlane> fittedState_;

 public:
  ClassDefOverride(GSFFitterInfo, 1)
};

} /* End of namespace genfit */
/** @} */