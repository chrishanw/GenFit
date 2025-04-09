 
#include <Rtypes.h>

namespace genfit {
  class Track;
  class AbsTrackRep;
  class KalmanFitStatus;

  //! Check if track has a KalmanFitStatus for given AbsTrackRep. Per default, check for cardinal rep.
  static bool hasTrackKalmanFitStatus(const Track* track, const AbsTrackRep* rep);

  //! If FitStatus is a KalmanFitStatus, return it. Otherwise return nullptr
  static KalmanFitStatus* getTrackKalmanFitStatus(const Track* track, const AbsTrackRep* rep = nullptr);

  //! Check track consistency, in this case only for KalmanFitterInfo
  static void checkTrackConistency(const Track* track);

  //! Helper function: For all KalmanFitterInfos belonging to rep (if nullptr, for all reps),
  //! call the fixWeights() function, so that e.g. the DAF will not alter weights anymore.
  static void fixTrackWeights(const Track* track, AbsTrackRep* rep = nullptr, int startId = 0, int endId = -1);

  /**
   * @brief Delete unneeded information from the Track.
   *
   * Possible options: (see also PruneFlags defined in FitStatus.h)
   * U:  if fitterInfo is a KalmanFitterInfo, prune predictions and keep updates
   */
  static void pruneTrack(const Track* track, const Option_t* option = "U");

}