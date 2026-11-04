#include "GSF.h"

#include "AbsMeasurement.h"
#include "AbsTrackRep.h"
#include "Exception.h"
#include "FitStatus.h"
#include "GSFFitStatus.h"
#include "GSFFitterInfo.h"
#include "IO.h"
#include "MeasurementOnPlane.h"
#include "MeasuredStateOnPlane.h"
#include "Track.h"
#include "TrackPoint.h"
#include "Tools.h"

#include <TDecompChol.h>
#include <TMatrixDSymEigen.h>
#include <TMath.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace genfit;

namespace {

using GaussianComponent = GSFComponent;
using Components = std::vector<GaussianComponent>;

struct ComponentUpdate {
  GaussianComponent component;
  double logLikelihood;
  double chi2;
  unsigned int measurementDimension;

  ComponentUpdate(const MeasuredStateOnPlane& state, double weight, double logLikelihoodIn,
                  double chi2In, unsigned int dimension)
      : component(state, weight), logLikelihood(logLikelihoodIn), chi2(chi2In),
        measurementDimension(dimension) {}
};

Components makeSeedMixture(const MeasuredStateOnPlane& seed, unsigned int maxComponents) {
  const TMatrixDSym& covariance = seed.getCov();
  const int dimension = covariance.GetNrows();
  TMatrixDSymEigen eigen(covariance);
  const TVectorD eigenvalues = eigen.GetEigenValues();
  const TMatrixD eigenvectors = eigen.GetEigenVectors();

  std::vector<int> directions;
  double largestEigenvalue = 0.;
  for (int i = 0; i < dimension; ++i)
    largestEigenvalue = std::max(largestEigenvalue, eigenvalues(i));
  for (int i = 0; i < dimension; ++i) {
    if (eigenvalues(i) > std::max(1e-14, largestEigenvalue * 1e-12))
      directions.push_back(i);
  }

  const unsigned int directionCount = std::min<unsigned int>(
      static_cast<unsigned int>(directions.size()), maxComponents / 2);
  if (directionCount == 0)
    return Components(1, GaussianComponent(seed, 1.));

  std::sort(directions.begin(), directions.end(), [&eigenvalues](int lhs, int rhs) {
    return eigenvalues(lhs) > eigenvalues(rhs);
  });
  directions.resize(directionCount);

  // Symmetric offsets carry the selected eigenmode variance; the shared
  // component covariance carries the remainder. Together they preserve the
  // seed's first two moments, including any modes omitted by the cap.
  TMatrixDSym componentCovariance(covariance);
  std::vector<TVectorD> offsets;
  offsets.reserve(directionCount);
  for (unsigned int i = 0; i < directionCount; ++i) {
    const int axis = directions[i];
    TVectorD offset(dimension);
    const double amplitude = std::sqrt(directionCount * eigenvalues(axis) / 2.);
    for (int row = 0; row < dimension; ++row)
      offset(row) = amplitude * eigenvectors(row, axis);
    offsets.push_back(offset);

    for (int row = 0; row < dimension; ++row) {
      for (int column = 0; column < dimension; ++column) {
        componentCovariance(row, column) -=
            0.5 * eigenvalues(axis) * eigenvectors(row, axis) * eigenvectors(column, axis);
      }
    }
  }

  Components components;
  components.reserve(2 * directionCount);
  const double weight = 1. / (2. * directionCount);
  for (const TVectorD& offset : offsets) {
    for (int sign : {-1, 1}) {
      MeasuredStateOnPlane state(seed);
      TVectorD componentState(seed.getState());
      componentState += static_cast<double>(sign) * offset;
      state.setStateCovPlane(componentState, componentCovariance, seed.getPlane());
      components.emplace_back(state, weight);
    }
  }
  return components;
}

MeasuredStateOnPlane momentMatch(const Components& components) {
  const MeasuredStateOnPlane& reference = components.front().state;
  const int dimension = reference.getState().GetNrows();
  TVectorD mean(dimension);
  mean.Zero();
  double totalWeight = 0.;
  for (const GaussianComponent& component : components) {
    mean += component.weight * component.state.getState();
    totalWeight += component.weight;
  }
  mean *= 1. / totalWeight;

  // Apply the law of total covariance: average each component covariance and
  // add the covariance of the component means.
  TMatrixDSym covariance(dimension);
  covariance.Zero();
  for (const GaussianComponent& component : components) {
    const TVectorD delta(component.state.getState() - mean);
    const double weight = component.weight / totalWeight;
    for (int row = 0; row < dimension; ++row) {
      for (int column = 0; column < dimension; ++column) {
        covariance(row, column) += weight *
            (component.state.getCov()(row, column) + delta(row) * delta(column));
      }
    }
  }

  MeasuredStateOnPlane result(reference);
  result.setStateCovPlane(mean, covariance, reference.getPlane());
  return result;
}

double componentDistance(const GaussianComponent& lhs, const GaussianComponent& rhs) {
  // A diagonal Mahalanobis distance is inexpensive and avoids another matrix
  // inversion in the repeated pair-selection loop.
  const TVectorD difference(lhs.state.getState() - rhs.state.getState());
  double distance = 0.;
  for (int i = 0; i < difference.GetNrows(); ++i) {
    const double variance = lhs.state.getCov()(i, i) + rhs.state.getCov()(i, i);
    distance += difference(i) * difference(i) / std::max(variance, 1e-15);
  }
  return distance;
}

GaussianComponent mergeComponents(const GaussianComponent& lhs, const GaussianComponent& rhs) {
  const double weight = lhs.weight + rhs.weight;
  const double lhsFraction = lhs.weight / weight;
  const double rhsFraction = rhs.weight / weight;
  const TVectorD mean(lhsFraction * lhs.state.getState() + rhsFraction * rhs.state.getState());
  const TVectorD lhsDelta(lhs.state.getState() - mean);
  const TVectorD rhsDelta(rhs.state.getState() - mean);
  TMatrixDSym covariance(lhs.state.getCov());
  covariance *= lhsFraction;
  TMatrixDSym rhsCovariance(rhs.state.getCov());
  rhsCovariance *= rhsFraction;
  covariance += rhsCovariance;
  for (int row = 0; row < mean.GetNrows(); ++row) {
    for (int column = 0; column < mean.GetNrows(); ++column) {
      covariance(row, column) += lhsFraction * lhsDelta(row) * lhsDelta(column) +
          rhsFraction * rhsDelta(row) * rhsDelta(column);
    }
  }
  MeasuredStateOnPlane state(lhs.state);
  state.setStateCovPlane(mean, covariance, lhs.state.getPlane());
  return GaussianComponent(state, weight);
}

void compressComponents(Components& components, unsigned int maxComponents, double pruningThreshold) {
  const auto largestComponent = std::max_element(components.begin(), components.end(),
      [](const GaussianComponent& lhs, const GaussianComponent& rhs) {
        return lhs.weight < rhs.weight;
      });
  GaussianComponent bestComponent = *largestComponent;
  components.erase(std::remove_if(components.begin(), components.end(), [pruningThreshold](const GaussianComponent& component) {
    return component.weight < pruningThreshold;
  }), components.end());
  if (components.empty())
    components.push_back(bestComponent);

  // Greedily merge the closest pair, preserving its total weight, mean, and
  // covariance. This bounds later branching while retaining those moments.
  while (components.size() > maxComponents) {
    std::size_t bestLeft = 0;
    std::size_t bestRight = 1;
    double bestDistance = componentDistance(components[bestLeft], components[bestRight]);
    for (std::size_t left = 0; left < components.size(); ++left) {
      for (std::size_t right = left + 1; right < components.size(); ++right) {
        const double distance = componentDistance(components[left], components[right]);
        if (distance < bestDistance) {
          bestDistance = distance;
          bestLeft = left;
          bestRight = right;
        }
      }
    }
    components[bestLeft] = mergeComponents(components[bestLeft], components[bestRight]);
    components.erase(components.begin() + bestRight);
  }

  double totalWeight = 0.;
  for (const GaussianComponent& component : components)
    totalWeight += component.weight;
  for (GaussianComponent& component : components)
    component.weight /= totalWeight;
}

ComponentUpdate updateComponent(const GaussianComponent& component, const MeasurementOnPlane& measurement) {
  const double measurementWeight = measurement.getWeight();
  if (measurementWeight <= 1e-10)
    throw std::runtime_error("Measurement has negligible weight");

  const TVectorD& measurementState = measurement.getState();
  const AbsHMatrix* hMatrix = measurement.getHMatrix();
  TMatrixDSym measurementCovariance(measurement.getCov());
  if (measurementWeight < 0.99999)
    measurementCovariance *= 1. / measurementWeight;

  const TVectorD& stateVector = component.state.getState();
  const TMatrixDSym& covariance = component.state.getCov();
  const TVectorD residual(measurementState - hMatrix->Hv(stateVector));

  // Apply a linear-Gaussian update to each component. The innovation
  // determinant and chi-square provide its measurement likelihood.
  TMatrixDSym innovationCovariance(covariance);
  hMatrix->HMHt(innovationCovariance);
  innovationCovariance += measurementCovariance;
  TDecompChol decomposition(innovationCovariance);
  if (!decomposition.Decompose())
    throw std::runtime_error("Could not decompose innovation covariance");
  const TMatrixD& upper = decomposition.GetU();
  double logDeterminant = 0.;
  for (int i = 0; i < upper.GetNrows(); ++i)
    logDeterminant += 2. * std::log(upper(i, i));

  tools::invertMatrix(innovationCovariance);
  const double chi2 = tools::similarity(residual, innovationCovariance);
  const TMatrixD covarianceTimesHt(hMatrix->MHt(covariance));
  const TVectorD update(TMatrixD(covarianceTimesHt, TMatrixD::kMult, innovationCovariance) * residual);
  TVectorD updatedState(stateVector + update);

  TMatrixDSym covarianceReduction(covariance);
  hMatrix->HMHt(covarianceReduction);
  covarianceReduction += measurementCovariance;
  tools::invertMatrix(covarianceReduction);
  tools::similarity(covarianceTimesHt, covarianceReduction);
  TMatrixDSym updatedCovariance(covariance);
  updatedCovariance -= covarianceReduction;

  MeasuredStateOnPlane state(component.state);
  state.setStateCovPlane(updatedState, updatedCovariance, component.state.getPlane());
  const double logLikelihood = -0.5 * (chi2 + logDeterminant +
      measurementState.GetNrows() * std::log(2. * TMath::Pi())) + std::log(measurementWeight);
  return ComponentUpdate(state, component.weight, logLikelihood, chi2,
                         static_cast<unsigned int>(measurementState.GetNrows()));
}

void normalizeUpdates(std::vector<ComponentUpdate>& updates) {
  // Normalize in log space so unlikely hypotheses do not underflow before
  // their relative probabilities are computed.
  double largestLogWeight = -std::numeric_limits<double>::infinity();
  for (const ComponentUpdate& update : updates)
    largestLogWeight = std::max(largestLogWeight, std::log(update.component.weight) + update.logLikelihood);
  double normalization = 0.;
  for (const ComponentUpdate& update : updates)
    normalization += std::exp(std::log(update.component.weight) + update.logLikelihood - largestLogWeight);
  const double logNormalization = largestLogWeight + std::log(normalization);
  for (ComponentUpdate& update : updates)
    update.component.weight = std::exp(std::log(update.component.weight) + update.logLikelihood - logNormalization);
}

} // namespace

GSF::GSF(unsigned int maxComponents, double pruningThreshold)
    : maxComponents_(maxComponents), pruningThreshold_(pruningThreshold) {
  setMaxComponents(maxComponents);
  setPruningThreshold(pruningThreshold);
}

void GSF::setMaxComponents(unsigned int value) {
  if (value == 0)
    throw std::invalid_argument("GSF requires at least one Gaussian component");
  maxComponents_ = value;
}

void GSF::setPruningThreshold(double value) {
  if (value < 0. || value >= 1.)
    throw std::invalid_argument("GSF pruning threshold must be in [0, 1)");
  pruningThreshold_ = value;
}

void GSF::processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool) {
  if (tr->hasFitStatus(rep) && tr->getFitStatus(rep)->isTrackPruned()) {
    Exception exc("GSF::processTrackWithRep: Cannot process pruned track!", __LINE__, __FILE__);
    throw exc;
  }

  const unsigned int nPoints = tr->getNumPointsWithMeasurement();
  GSFFitStatus* status = new GSFFitStatus();
  status->setMaxComponents(maxComponents_);
  tr->setFitStatus(status, rep);
  if (nPoints == 0) {
    status->setIsFitted(false);
    return;
  }

  MeasuredStateOnPlane seed(rep);
  rep->setTime(seed, tr->getTimeSeed());
  rep->setPosMomCov(seed, tr->getStateSeed(), tr->getCovSeed());
  Components components = makeSeedMixture(seed, maxComponents_);
  status->setMaxComponentsUsed(static_cast<unsigned int>(components.size()));

  double forwardChi2 = 0.;
  double forwardNdf = -static_cast<double>(rep->getDim());
  double backwardChi2 = 0.;
  double backwardNdf = -static_cast<double>(rep->getDim());
  int failedPointsForward = 0;
  int failedPointsBackward = 0;
  std::vector<bool> forwardPointValid(nPoints, false);

  // Both passes use the same component update. The reverse pass reuses the
  // planes and measurements built during the forward pass.
  auto runPass = [&](int direction, double& chi2, double& ndf, int& failedPoints) {
    for (int pointIndex = direction > 0 ? 0 : static_cast<int>(nPoints) - 1;
         pointIndex >= 0 && pointIndex < static_cast<int>(nPoints); pointIndex += direction) {
      TrackPoint* trackPoint = tr->getPointWithMeasurement(pointIndex);
      if (!trackPoint->hasRawMeasurements())
        continue;

      if (direction < 0 && !forwardPointValid[pointIndex])
        continue;

      GSFFitterInfo* fitterInfo = nullptr;
      if (trackPoint->hasFitterInfo(rep))
        fitterInfo = dynamic_cast<GSFFitterInfo*>(trackPoint->getFitterInfo(rep));
      if (fitterInfo == nullptr) {
        if (trackPoint->hasFitterInfo(rep))
          trackPoint->deleteFitterInfo(rep);
        fitterInfo = new GSFFitterInfo(trackPoint, rep);
        trackPoint->setFitterInfo(fitterInfo);
      }

      Components componentsBeforePoint(components);
      try {
        SharedPlanePtr plane;
        if (direction > 0) {
          fitterInfo->deleteForwardInfo();
          fitterInfo->deleteBackwardInfo();
          fitterInfo->deleteMeasurementInfo();
          plane = trackPoint->getRawMeasurements().front()->constructPlane(components.front().state);
        }
        else {
          plane = fitterInfo->getPlane();
        }

        if (!plane)
          throw std::runtime_error("TrackPoint has no detector plane");

        for (GaussianComponent& component : components)
          rep->extrapolateToPlane(component.state, plane);

        fitterInfo->setComponents(components, direction, false);
        MeasuredStateOnPlane prediction = momentMatch(components);
        fitterInfo->setPrediction(new MeasuredStateOnPlane(prediction), direction);

        if (direction > 0) {
          for (AbsMeasurement* rawMeasurement : trackPoint->getRawMeasurements()) {
            for (MeasurementOnPlane* measurement : rawMeasurement->constructMeasurementsOnPlane(prediction)) {
              fitterInfo->addMeasurementOnPlane(*measurement);
              delete measurement;
            }
          }
        }

        // Form the component x measurement-hypothesis branches, then normalize
        // their posterior weights before pruning or moment matching.
        std::vector<ComponentUpdate> updates;
        for (const GaussianComponent& component : components) {
          for (const MeasurementOnPlane& measurement : fitterInfo->getMeasurementsOnPlane()) {
            try {
              updates.push_back(updateComponent(component, measurement));
            }
            catch (const std::runtime_error&) {
            }
          }
        }
        if (updates.empty())
          throw std::runtime_error("No valid component measurement updates");

        normalizeUpdates(updates);
        status->setMaxComponentsUsed(std::max(status->getMaxComponentsUsed(),
          static_cast<unsigned int>(updates.size())));
        double pointChi2 = 0.;
        double pointNdf = 0.;
        components.clear();
        components.reserve(updates.size());
        for (const ComponentUpdate& update : updates) {
          pointChi2 += update.component.weight * update.chi2;
          pointNdf += update.component.weight * update.measurementDimension;
          components.push_back(update.component);
        }
        compressComponents(components, maxComponents_, pruningThreshold_);
        if (components.empty())
          throw std::runtime_error("All Gaussian components were pruned");

        MeasuredStateOnPlane updatedState = momentMatch(components);
        fitterInfo->setComponents(components, direction);
        fitterInfo->setUpdate(updatedState, direction, pointChi2, pointNdf);
        chi2 += pointChi2;
        ndf += pointNdf;
        if (direction > 0)
          forwardPointValid[pointIndex] = true;
      }
      catch (Exception& exception) {
        components = componentsBeforePoint;
        errorOut << exception.what();
        ++failedPoints;
      }
      catch (const std::runtime_error& exception) {
        components = componentsBeforePoint;
        errorOut << "GSF failed at TrackPoint " << pointIndex << ": "
                 << exception.what() << std::endl;
        ++failedPoints;
      }
    }
  };

  runPass(1, forwardChi2, forwardNdf, failedPointsForward);
  if (failedPointsForward == static_cast<int>(nPoints)) {
    status->setIsFitted(false);
    status->setNFailedPoints(failedPointsForward);
    return;
  }

  // Inflate the forward covariance before the reverse pass so the backward
  // update is not constrained by the forward result.
  for (GaussianComponent& component : components)
    component.state.blowUpCov(1e3, true, 1e6);
  runPass(-1, backwardChi2, backwardNdf, failedPointsBackward);
  const int failedPoints = std::max(failedPointsForward, failedPointsBackward);

  status->setIsFitted(true);
  status->setIsFitConvergedFully(failedPoints == 0);
  status->setIsFitConvergedPartially(true);
  status->setNFailedPoints(failedPoints);
  status->setChi2(backwardChi2);
  status->setNdf(std::max(0., backwardNdf));
  status->setForwardChi2(forwardChi2);
  status->setForwardNdf(std::max(0., forwardNdf));
  status->setNumIterations(1);
  status->setCharge(rep->getCharge(components.front().state));
}