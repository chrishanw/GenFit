#include "GSFFitterInfo.h"

#include "AbsTrackRep.h"
#include "Exception.h"
#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"

#include <cmath>

namespace genfit {

GSFFitterInfo::GSFFitterInfo() : AbsFitterInfo() {}

GSFFitterInfo::GSFFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep)
    : AbsFitterInfo(trackPoint, rep) {}

GSFFitterInfo* GSFFitterInfo::clone() const {
  GSFFitterInfo* result = new GSFFitterInfo(getTrackPoint(), getRep());
  if (hasForwardPrediction())
    result->setForwardPrediction(new MeasuredStateOnPlane(*getForwardPrediction()));
  if (hasForwardUpdate())
    result->setForwardUpdate(*getForwardUpdate(), getForwardUpdateChi2Increment(), getForwardUpdateNdf());
  if (hasBackwardPrediction())
    result->setBackwardPrediction(new MeasuredStateOnPlane(*getBackwardPrediction()));
  if (hasBackwardUpdate())
    result->setBackwardUpdate(*getBackwardUpdate(), getBackwardUpdateChi2Increment(), getBackwardUpdateNdf());

  for (const MeasurementOnPlane& measurement : measurementsOnPlane_)
    result->addMeasurementOnPlane(measurement);
  if (!measurementsOnPlane_.empty())
    result->setWeights(getWeights());
  result->fixWeights(fixWeights_);
  result->setPlane(getPlane());

  result->forwardPredictedComponents_ = forwardPredictedComponents_;
  result->forwardUpdatedComponents_ = forwardUpdatedComponents_;
  result->backwardPredictedComponents_ = backwardPredictedComponents_;
  result->backwardUpdatedComponents_ = backwardUpdatedComponents_;
  return result;
}

std::vector<double> GSFFitterInfo::getWeights() const {
  std::vector<double> weights;
  weights.reserve(measurementsOnPlane_.size());
  for (const MeasurementOnPlane& measurement : measurementsOnPlane_)
    weights.push_back(measurement.getWeight());
  return weights;
}

void GSFFitterInfo::setForwardPrediction(MeasuredStateOnPlane* prediction) {
  if (prediction != nullptr)
    setPlane(prediction->getPlane());
  forwardPrediction_.reset(prediction);
}

void GSFFitterInfo::setBackwardPrediction(MeasuredStateOnPlane* prediction) {
  if (prediction != nullptr)
    setPlane(prediction->getPlane());
  backwardPrediction_.reset(prediction);
}

void GSFFitterInfo::setPrediction(MeasuredStateOnPlane* prediction, int direction) {
  if (direction >= 0)
    setForwardPrediction(prediction);
  else
    setBackwardPrediction(prediction);
}

void GSFFitterInfo::setForwardUpdate(const MeasuredStateOnPlane& update, double chi2Increment, double ndf) {
  forwardUpdate_.reset(new MeasuredStateOnPlane(update));
  forwardUpdateChi2Increment_ = chi2Increment;
  forwardUpdateNdf_ = ndf;
  setPlane(update.getPlane());
}

void GSFFitterInfo::setBackwardUpdate(const MeasuredStateOnPlane& update, double chi2Increment, double ndf) {
  backwardUpdate_.reset(new MeasuredStateOnPlane(update));
  backwardUpdateChi2Increment_ = chi2Increment;
  backwardUpdateNdf_ = ndf;
  setPlane(update.getPlane());
}

void GSFFitterInfo::setUpdate(const MeasuredStateOnPlane& update, int direction, double chi2Increment, double ndf) {
  if (direction >= 0)
    setForwardUpdate(update, chi2Increment, ndf);
  else
    setBackwardUpdate(update, chi2Increment, ndf);
}

void GSFFitterInfo::addMeasurementOnPlane(const MeasurementOnPlane& measurement) {
  if (measurementsOnPlane_.empty())
    setPlane(measurement.getPlane());
  measurementsOnPlane_.push_back(measurement);
}

void GSFFitterInfo::setWeights(const std::vector<double>& weights) {
  if (weights.size() != measurementsOnPlane_.size()) {
    Exception exception("GSFFitterInfo::setWeights: weight count does not match measurement count.", __LINE__, __FILE__);
    throw exception;
  }
  for (std::size_t i = 0; i < weights.size(); ++i)
    measurementsOnPlane_[i].setWeight(weights[i]);
}

const std::vector<GSFComponent>& GSFFitterInfo::getComponents(int direction, bool updated) const {
  if (direction >= 0)
    return updated ? forwardUpdatedComponents_ : forwardPredictedComponents_;
  return updated ? backwardUpdatedComponents_ : backwardPredictedComponents_;
}

void GSFFitterInfo::setComponents(const std::vector<GSFComponent>& components, int direction, bool updated) {
  std::vector<GSFComponent>& destination = direction >= 0
      ? (updated ? forwardUpdatedComponents_ : forwardPredictedComponents_)
      : (updated ? backwardUpdatedComponents_ : backwardPredictedComponents_);
  destination = components;
}

void GSFFitterInfo::setRep(const AbsTrackRep* rep) {
  AbsFitterInfo::setRep(rep);
  if (forwardPrediction_)
    forwardPrediction_->setRep(rep);
  if (backwardPrediction_)
    backwardPrediction_->setRep(rep);
  if (forwardUpdate_)
    forwardUpdate_->setRep(rep);
  if (backwardUpdate_)
    backwardUpdate_->setRep(rep);
  for (MeasurementOnPlane& measurement : measurementsOnPlane_)
    measurement.setRep(rep);
  for (std::vector<GSFComponent>* components : {&forwardPredictedComponents_, &forwardUpdatedComponents_,
                                                &backwardPredictedComponents_, &backwardUpdatedComponents_}) {
    for (GSFComponent& component : *components)
      component.state.setRep(rep);
  }
}

void GSFFitterInfo::deleteForwardInfo() {
  forwardPrediction_.reset();
  forwardUpdate_.reset();
  forwardUpdateChi2Increment_ = 0.;
  forwardUpdateNdf_ = 0.;
  forwardPredictedComponents_.clear();
  forwardUpdatedComponents_.clear();
}

void GSFFitterInfo::deleteBackwardInfo() {
  backwardPrediction_.reset();
  backwardUpdate_.reset();
  backwardUpdateChi2Increment_ = 0.;
  backwardUpdateNdf_ = 0.;
  backwardPredictedComponents_.clear();
  backwardUpdatedComponents_.clear();
}

void GSFFitterInfo::deleteMeasurementInfo() {
  measurementsOnPlane_.clear();
}

const MeasuredStateOnPlane& GSFFitterInfo::getFittedState(bool biased) const {
  const TrackPoint* trackPoint = getTrackPoint();
  if (trackPoint == nullptr || trackPoint->getTrack() == nullptr) {
    Exception exception("GSFFitterInfo::getFittedState: no owning TrackPoint or Track.", __LINE__, __FILE__);
    throw exception;
  }

  const Track* track = trackPoint->getTrack();
  bool first = false;
  bool last = false;
  FitStatus* status = track->getFitStatus(getRep());
  if (status != nullptr && status->isTrackPruned() && track->getNumPoints() == 1) {
    first = status->getPruneFlags().hasFlags("F");
    last = status->getPruneFlags().hasFlags("L");
  } else {
    first = track->getPointWithFitterInfo(0, getRep()) == trackPoint;
    last = track->getPointWithFitterInfo(-1, getRep()) == trackPoint;
  }

  if (biased) {
    if (last && !hasBackwardPrediction() && hasForwardUpdate())
      return *forwardUpdate_;
    if (first && (!hasForwardUpdate() || (hasBackwardUpdate() && !hasForwardPrediction()))) {
      if (hasBackwardUpdate())
        return *backwardUpdate_;
    }
    if (hasForwardUpdate() && hasBackwardPrediction()) {
      fittedState_ = calcAverageState(*forwardUpdate_, *backwardPrediction_);
      return fittedState_;
    }
  } else {
    if (last && !hasBackwardPrediction() && hasForwardPrediction())
      return *forwardPrediction_;
    if (first && !hasForwardPrediction() && hasBackwardPrediction())
      return *backwardPrediction_;
    if (hasForwardPrediction() && hasBackwardPrediction()) {
      fittedState_ = calcAverageState(*forwardPrediction_, *backwardPrediction_);
      return fittedState_;
    }
  }

  Exception exception("GSFFitterInfo::getFittedState: required predictions or updates are unavailable.", __LINE__, __FILE__);
  throw exception;
}

MeasurementOnPlane GSFFitterInfo::getResidual(unsigned int iMeasurement, bool biased, bool onlyMeasurementErrors) const {

  const MeasuredStateOnPlane& fittedState = getFittedState(biased);
  const MeasurementOnPlane& measurement = measurementsOnPlane_.at(iMeasurement);
  const SharedPlanePtr& plane = measurement.getPlane();

  // check equality of planes and reps
  if(*(fittedState.getPlane()) != *plane) {
    Exception exception("GSFFitterInfo::getResidual: fittedState and measurement are not defined in the same plane.", __LINE__,__FILE__);
    throw exception;
  }
  if(fittedState.getRep() != measurement.getRep()) {
    Exception exception("GSFFitterInfo::getResidual: fittedState and measurement are not defined wrt the same TrackRep.", __LINE__,__FILE__);
    throw exception;
  }

  const AbsHMatrix* hMatrix = measurement.getHMatrix();

  TVectorD residual(hMatrix->Hv(fittedState.getState()) - measurement.getState());
  residual *= -1.;
  TMatrixDSym covariance(measurement.getCov());

  if (onlyMeasurementErrors) {
    return MeasurementOnPlane(residual, covariance, plane, fittedState.getRep(), hMatrix->clone(), measurement.getWeight());
  }

  TMatrixDSym stateCov(fittedState.getCov());
  hMatrix->HMHt(stateCov);
  covariance += stateCov;
  return MeasurementOnPlane(residual, covariance, measurement.getPlane(), fittedState.getRep(),
                             hMatrix->clone(), measurement.getWeight());
}

bool GSFFitterInfo::checkConsistency(const PruneFlags*) const {
  if (!measurementsOnPlane_.empty() && getPlane() != measurementsOnPlane_.front().getPlane())
    return false;
  for (const MeasurementOnPlane& measurement : measurementsOnPlane_) {
    if (measurement.getRep() != getRep() || measurement.getPlane() != getPlane())
      return false;
  }
  const auto validComponents = [](const std::vector<GSFComponent>& components) {
    double weightSum = 0.;
    for (const GSFComponent& component : components) {
      if (component.weight < 0. || !std::isfinite(component.weight))
        return false;
      weightSum += component.weight;
    }
    return components.empty() || std::fabs(weightSum - 1.) < 1e-8;
  };
  return validComponents(forwardPredictedComponents_) && validComponents(forwardUpdatedComponents_) &&
         validComponents(backwardPredictedComponents_) && validComponents(backwardUpdatedComponents_);
}

} /* End of namespace genfit */
