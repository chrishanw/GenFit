#include "GSFFitterInfo.h"

#include "AbsTrackRep.h"
#include "Exception.h"
#include "IO.h"
#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"

#include <cmath>

namespace genfit {

GSFFitterInfo::GSFFitterInfo() : AbsFitterInfo() {}

GSFFitterInfo::GSFFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep)
    : AbsFitterInfo(trackPoint, rep) {}

GSFFitterInfo::~GSFFitterInfo() {
  deleteMeasurementInfo();
}

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

  for (const MeasurementOnPlane* measurement : measurementsOnPlane_)
    result->addMeasurementOnPlane(*measurement);
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
  for (const MeasurementOnPlane* measurement : measurementsOnPlane_)
    weights.push_back(measurement->getWeight());
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
  measurementsOnPlane_.push_back(new MeasurementOnPlane(measurement));
}

void GSFFitterInfo::setWeights(const std::vector<double>& weights) {
  if (weights.size() != measurementsOnPlane_.size()) {
    Exception exception("GSFFitterInfo::setWeights: weight count does not match measurement count.", __LINE__, __FILE__);
    throw exception;
  }
  for (std::size_t i = 0; i < weights.size(); ++i)
    measurementsOnPlane_[i]->setWeight(weights[i]);
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
  for (MeasurementOnPlane* measurement : measurementsOnPlane_)
    measurement->setRep(rep);
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
  for (MeasurementOnPlane* measurement : measurementsOnPlane_)
    delete measurement;
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
      fittedState_.reset(new MeasuredStateOnPlane(calcAverageState(*forwardUpdate_, *backwardPrediction_)));
      return *fittedState_;
    }
  } else {
    if (last && !hasBackwardPrediction() && hasForwardPrediction())
      return *forwardPrediction_;
    if (first && !hasForwardPrediction() && hasBackwardPrediction())
      return *backwardPrediction_;
    if (hasForwardPrediction() && hasBackwardPrediction()) {
      fittedState_.reset(new MeasuredStateOnPlane(calcAverageState(*forwardPrediction_, *backwardPrediction_)));
      return *fittedState_;
    }
  }

  Exception exception("GSFFitterInfo::getFittedState: required predictions or updates are unavailable.", __LINE__, __FILE__);
  throw exception;
}

MeasurementOnPlane GSFFitterInfo::getResidual(unsigned int iMeasurement, bool biased, bool onlyMeasurementErrors) const {

  const MeasuredStateOnPlane& fittedState = getFittedState(biased);
  const MeasurementOnPlane& measurement = *measurementsOnPlane_.at(iMeasurement);
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

bool GSFFitterInfo::checkConsistency(const PruneFlags* flags) const {
  bool retVal(true);

  // check if in a TrackPoint
  if (trackPoint_ == nullptr) {
    errorOut << "GSFFitterInfo::checkConsistency(): trackPoint_ is nullptr" << std::endl;
    retVal = false;
  }

  // a rep is required for the plane/rep/dimension checks below
  if (getRep() == nullptr) {
    errorOut << "GSFFitterInfo::checkConsistency(): rep_ is nullptr" << std::endl;
    return false;
  }

  const SharedPlanePtr plane = getPlane();

  // An empty GSFFitterInfo is a valid state
  if (plane.get() == nullptr) {
    if (!(forwardPrediction_ || forwardUpdate_ || backwardPrediction_ || backwardUpdate_ ||
          !measurementsOnPlane_.empty() ||
          !forwardPredictedComponents_.empty() || !forwardUpdatedComponents_.empty() ||
          !backwardPredictedComponents_.empty() || !backwardUpdatedComponents_.empty()))
      return true;
    errorOut << "GSFFitterInfo::checkConsistency(): plane is nullptr" << std::endl;
    retVal = false;
  }
  else {
    // see if the plane object is actually there
    // cppcheck-suppress unreadVariable
    TVector3 oTest = plane->getO();
    // cppcheck-suppress unreadVariable
    oTest *= 47.11;
  }

  const int dim = getRep()->getDim();

  // check a stored state: non-null, in the correct plane and rep, right dimension
  const auto checkState = [&](const MeasuredStateOnPlane* state, const char* name) {
    if (state == nullptr)
      return;   // an absent state is not an inconsistency
    if (state->getPlane().get() == nullptr) {
      errorOut << "GSFFitterInfo::checkConsistency(): " << name << " has a nullptr plane" << std::endl;
      retVal = false;
    }
    else if (state->getPlane() != plane) {
      errorOut << "GSFFitterInfo::checkConsistency(): " << name
               << " is not defined with the correct plane "
               << state->getPlane().get() << " vs. " << plane.get() << std::endl;
      retVal = false;
    }
    if (state->getRep() != getRep()) {
      errorOut << "GSFFitterInfo::checkConsistency(): " << name
               << " is not defined with the correct TrackRep" << std::endl;
      retVal = false;
    }
    if (state->getState().GetNrows() != dim || state->getCov().GetNrows() != dim) {
      errorOut << "GSFFitterInfo::checkConsistency(): " << name
               << " does not have the right dimension!" << std::endl;
      retVal = false;
    }
  };

  checkState(forwardPrediction_.get(),  "forwardPrediction_");
  checkState(forwardUpdate_.get(),      "forwardUpdate_");
  checkState(backwardPrediction_.get(), "backwardPrediction_");
  checkState(backwardUpdate_.get(),     "backwardUpdate_");

  // measurements
  for (unsigned int i = 0; i < measurementsOnPlane_.size(); ++i) {
    const MeasurementOnPlane* measurement = measurementsOnPlane_[i];
    if (measurement->getPlane() != plane) {
      errorOut << "GSFFitterInfo::checkConsistency(): measurement " << i
               << " is not defined with the correct plane" << std::endl;
      retVal = false;
    }
    if (measurement->getRep() != getRep()) {
      errorOut << "GSFFitterInfo::checkConsistency(): measurement " << i
               << " is not defined with the correct TrackRep" << std::endl;
      retVal = false;
    }
    if (measurement->getState().GetNrows() == 0) {
      errorOut << "GSFFitterInfo::checkConsistency(): measurement " << i
               << " has dimension 0!" << std::endl;
      retVal = false;
    }
  }

  // if more than one measurement, check that they share a dimensionality
  if (measurementsOnPlane_.size() > 1) {
    const int mDim = measurementsOnPlane_[0]->getState().GetNrows();
    for (unsigned int i = 1; i < measurementsOnPlane_.size(); ++i) {
      if (measurementsOnPlane_[i]->getState().GetNrows() != mDim) {
        errorOut << "GSFFitterInfo::checkConsistency(): measurementsOnPlane_ do not all have the same dimensionality" << std::endl;
        retVal = false;
      }
    }
    if (mDim == 0) {
      errorOut << "GSFFitterInfo::checkConsistency(): measurementsOnPlane_ have dimension 0" << std::endl;
      retVal = false;
    }
  }

  // component mixtures: state plane/rep/dimension and normalized weights
  const auto validComponents = [&](const std::vector<GSFComponent>& components, const char* name) {
    if (components.empty())
      return;
    double weightSum = 0.;
    for (unsigned int i = 0; i < components.size(); ++i) {
      const GSFComponent& component = components[i];
      if (component.weight < 0. || !std::isfinite(component.weight)) {
        errorOut << "GSFFitterInfo::checkConsistency(): " << name << "[" << i
                 << "] has an invalid weight" << std::endl;
        retVal = false;
      }
      if (component.state.getPlane().get() == nullptr) {
        errorOut << "GSFFitterInfo::checkConsistency(): " << name << "[" << i
                 << "] has a nullptr plane" << std::endl;
        retVal = false;
      }
      else if (component.state.getPlane() != plane) {
        errorOut << "GSFFitterInfo::checkConsistency(): " << name << "[" << i
                 << "] is not defined with the correct plane" << std::endl;
        retVal = false;
      }
      if (component.state.getRep() != getRep()) {
        errorOut << "GSFFitterInfo::checkConsistency(): " << name << "[" << i
                 << "] is not defined with the correct TrackRep" << std::endl;
        retVal = false;
      }
      if (component.state.getState().GetNrows() != dim || component.state.getCov().GetNrows() != dim) {
        errorOut << "GSFFitterInfo::checkConsistency(): " << name << "[" << i
                 << "] does not have the right dimension!" << std::endl;
        retVal = false;
      }
      weightSum += component.weight;
    }
    if (std::fabs(weightSum - 1.) > 1e-8) {
      errorOut << "GSFFitterInfo::checkConsistency(): " << name
               << " weights do not sum to 1" << std::endl;
      retVal = false;
    }
  };

  validComponents(forwardPredictedComponents_,  "forwardPredictedComponents_");
  validComponents(forwardUpdatedComponents_,    "forwardUpdatedComponents_");
  validComponents(backwardPredictedComponents_, "backwardPredictedComponents_");
  validComponents(backwardUpdatedComponents_,   "backwardUpdatedComponents_");

  // update w/o prediction / measurement, unless the corresponding info was pruned
  if (flags == nullptr || !flags->hasFlags("U")) {
    if (forwardUpdate_ && !forwardPrediction_) {
      errorOut << "GSFFitterInfo::checkConsistency(): forwardUpdate_ w/o forwardPrediction_" << std::endl;
      retVal = false;
    }
    if (backwardUpdate_ && !backwardPrediction_) {
      errorOut << "GSFFitterInfo::checkConsistency(): backwardUpdate_ w/o backwardPrediction_" << std::endl;
      retVal = false;
    }
    if (flags == nullptr || !flags->hasFlags("M")) {
      if (forwardUpdate_ && measurementsOnPlane_.empty()) {
        errorOut << "GSFFitterInfo::checkConsistency(): forwardUpdate_ w/o measurement" << std::endl;
        retVal = false;
      }
      if (backwardUpdate_ && measurementsOnPlane_.empty()) {
        errorOut << "GSFFitterInfo::checkConsistency(): backwardUpdate_ w/o measurement" << std::endl;
        retVal = false;
      }
    }
  }

  return retVal;
}

} /* End of namespace genfit */
