/** @addtogroup genfit
 * @{
 */

#pragma once

#include "FitStatus.h"

#include <Math/ProbFunc.h>

namespace genfit {

/** Fit summary with diagnostics specific to the Gaussian-sum fitter. */
class GSFFitStatus : public FitStatus {

 public:
  GSFFitStatus() : FitStatus(), numIterations_(0), forwardChi2_(-1e99), forwardNdf_(-1e99),
      maxComponents_(0), maxComponentsUsed_(0) {}

  FitStatus* clone() const override { return new GSFFitStatus(*this); }

  unsigned int getNumIterations() const { return numIterations_; }
  double getForwardChi2() const { return forwardChi2_; }
  double getBackwardChi2() const { return FitStatus::getChi2(); }
  double getForwardNdf() const { return forwardNdf_; }
  double getBackwardNdf() const { return FitStatus::getNdf(); }
  double getForwardPVal() const { return std::max(0., ROOT::Math::chisquared_cdf_c(forwardChi2_, forwardNdf_)); }
  double getBackwardPVal() const { return FitStatus::getPVal(); }

  //! Configured upper bound on the number of mixture components.
  unsigned int getMaxComponents() const { return maxComponents_; }
  //! Largest mixture size encountered during this fit, before compression.
  unsigned int getMaxComponentsUsed() const { return maxComponentsUsed_; }

  void setNumIterations(unsigned int value) { numIterations_ = value; }
  void setForwardChi2(double value) { forwardChi2_ = value; }
  void setBackwardChi2(double value) { FitStatus::setChi2(value); }
  void setForwardNdf(double value) { forwardNdf_ = value; }
  void setBackwardNdf(double value) { FitStatus::setNdf(value); }
  void setMaxComponents(unsigned int value) { maxComponents_ = value; }
  void setMaxComponentsUsed(unsigned int value) { maxComponentsUsed_ = value; }

 private:
  unsigned int numIterations_;
  double forwardChi2_;
  double forwardNdf_;
  unsigned int maxComponents_;
  unsigned int maxComponentsUsed_;

 public:
  ClassDefOverride(GSFFitStatus, 1)
};

} /* End of namespace genfit */
/** @} */