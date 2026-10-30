/** @addtogroup genfit
 * @{
 */

#pragma once

#include "AbsFitter.h"

namespace genfit {

/**
 * @brief Gaussian-sum fitter with bounded, weighted Gaussian components.
 *
 * The seed Gaussian is represented by a symmetric set of components along
 * its leading covariance eigenmodes. Each component is propagated and updated
 * independently. At a TrackPoint with multiple raw measurements, every
 * component is branched against each measurement hypothesis and the posterior
 * weights include that hypothesis's measurement weight and Gaussian
 * likelihood. Low-weight components are pruned and the closest remaining
 * pairs are moment-merged until the configured component limit is met.
 *
 * The fitter makes one forward pass followed by one reverse pass; it does not
 * perform iterative refitting or hit resorting. Moment-matched directional
 * predictions and updates, together with the component mixtures, are stored
 * in GSFFitterInfo, which derives from AbsFitterInfo. Fit summaries use
 * GSFFitStatus, which derives from FitStatus.
 *
 * Material transport and noise are supplied by the AbsTrackRep. This class
 * does not introduce a separate Bethe-Heitler or other non-Gaussian material
 * interaction model; the explicit mixture hypotheses here come from the seed
 * approximation and multiple measurements at a TrackPoint.
 */
class GSF : public AbsFitter {

 public:
  /**
   * @param maxComponents Maximum number of components retained after
   *        measurement branching and compression.
   * @param pruningThreshold Remove components whose normalized weight is
   *        below this value before merging. If all components fall below the
   *        threshold, the highest-weight component is retained.
   */
  explicit GSF(unsigned int maxComponents = 12, double pruningThreshold = 1e-6);
  virtual ~GSF() = default;

  //! Fit one representation of a track. Hit resorting is not supported.
  void processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool resortHits = false) override;

  //! Get the maximum number of mixture components retained per pass step.
  unsigned int getMaxComponents() const { return maxComponents_; }
  //! Get the normalized component-weight threshold used for pruning.
  double getPruningThreshold() const { return pruningThreshold_; }
  //! Set the component limit; zero is invalid.
  void setMaxComponents(unsigned int value);
  //! Set a pruning threshold in [0, 1).
  void setPruningThreshold(double value);

 private:
  unsigned int maxComponents_;
  double pruningThreshold_;
};

} /* End of namespace genfit */
/** @} */