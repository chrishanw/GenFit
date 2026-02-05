 
#include <KalmanHelpers.h>
#include <Track.h>
#include <TrackPoint.h>
#include <AbsTrackRep.h>
#include <KalmanFitterInfo.h>

namespace genfit {

  static bool hasTrackKalmanFitStatus(const Track* track, const AbsTrackRep* rep) {
    if (rep == nullptr) {
      rep = track->getCardinalRep();
    }
    
    const auto& fitStatuses = track->getFitStatuses();
  
    if (fitStatuses.find(rep) == fitStatuses.end())
      return false;
  
    return (dynamic_cast<KalmanFitStatus*>(fitStatuses.at(rep)) != nullptr);
  }

  static KalmanFitStatus* getTrackKalmanFitStatus(const Track* track, const AbsTrackRep* rep) const {
    return dynamic_cast<KalmanFitStatus*>(track->getFitStatus(rep));
  }

  static void checkTrackConistency(const Track* track) {

    // cppcheck-suppress unreadVariable
    bool consistent = true;
    std::stringstream failures;
    const auto& trackPoints = track->getPoints();

    std::map<const AbsTrackRep*, const KalmanFitterInfo*> prevFis;

    // check TrackPoints
    for (std::vector<TrackPoint*>::const_iterator tp = trackPoints.begin(); tp != trackPoints.end(); ++tp) {

      // check fitterInfos
      const std::vector<AbsFitterInfo*>& fitterInfos = (*tp)->getFitterInfos();
  
      for (std::vector<AbsFitterInfo*>::const_iterator fi = fitterInfos.begin(); fi != fitterInfos.end(); ++fi) {

        if (dynamic_cast<KalmanFitterInfo*>(*fi) != nullptr) {
          if (prevFis[(*fi)->getRep()] != nullptr &&
              static_cast<KalmanFitterInfo*>(*fi)->hasReferenceState() &&
              prevFis[(*fi)->getRep()]->hasReferenceState() ) {
            const double len = static_cast<KalmanFitterInfo*>(*fi)->getReferenceState()->getForwardSegmentLength();
            const double prevLen = prevFis[(*fi)->getRep()]->getReferenceState()->getBackwardSegmentLength();
            if (fabs(prevLen + len) > 1E-10 ) {
              failures << "Track::checkConsistency(): segment lengths of reference states for rep " << (*fi)->getRep() << " (id " << getIdForRep((*fi)->getRep()) << ") at TrackPoint " << (*tp) << " don't match" << std::endl;
              failures << prevLen << " + " << len << " = " << prevLen + len << std::endl;
              failures << "TrackPoint " << *tp << ", FitterInfo " << *fi << ", rep " << getIdForRep((*fi)->getRep()) << std::endl;
        // cppcheck-suppress unreadVariable
              consistent = false;
            }
          }

          prevFis[(*fi)->getRep()] = static_cast<KalmanFitterInfo*>(*fi);
        }
        else
          prevFis[(*fi)->getRep()] = nullptr;

      } // end loop over FitterInfos

    } // end loop over TrackPoints

    if (not consistent) {
      throw genfit::Exception(failures.str(), __LINE__, __FILE__);
    }
  }


  static void fixTrackWeights(const Track* track, AbsTrackRep* rep, int startId, int endId) {

    const auto& trackPoints = track->getPoints();

    if (startId < 0)
      startId += trackPoints.size();
    if (endId < 0)
      endId += trackPoints.size();
  
    assert(startId >= 0);
    assert(startId <= endId);
    assert(endId <= (int)trackPoints.size());
  
    std::vector< AbsFitterInfo* > fis;
  
    for (std::vector<TrackPoint*>::iterator tp = trackPoints.begin() + startId; tp != trackPoints.begin() + endId; ++tp) {
      fis.clear();
      if (rep == nullptr) {
        fis = (*tp)->getFitterInfos();
      }
      else if ((*tp)->hasFitterInfo(rep)) {
        fis.push_back((*tp)->getFitterInfo(rep));
      }
  
      for (std::vector< AbsFitterInfo* >::iterator fi = fis.begin(); fi != fis.end(); ++fi) {
        KalmanFitterInfo* kfi = dynamic_cast<KalmanFitterInfo*>(*fi);
        if (kfi == nullptr)
          continue;
  
        kfi->fixWeights();
      }
    }
  }
  
  
  static void pruneTrack(const Track* track, const Option_t* option) {
  
    const auto& trackPoints = track->getPoints();

    PruneFlags f;
    f.setFlags(option);
  
    // from remaining trackPoints: prune measurementsOnPlane, unneeded fitterInfoStuff
    for (unsigned int i = 0; i<trackPoints.size(); ++i) {
      if (f.hasFlags("W"))
      trackPoints[i]->deleteRawMeasurements();
  
      std::vector< AbsFitterInfo* > fis =  trackPoints[i]->getFitterInfos();
      for (unsigned int j = 0; j<fis.size(); ++j) {
  
        if (f.hasFlags("U") && dynamic_cast<KalmanFitterInfo*>(fis[j]) != nullptr) {
          static_cast<KalmanFitterInfo*>(fis[j])->deletePredictions();
        }
      }
    }
  
    track->fillPointsWithMeasurement();
  
    #ifdef DEBUG
    debugOut << "pruned Track: "; Print();
    #endif
  
  }

  template<unsigned int dim, unsigned int dimAux>
  static KalmanFitterInfo<dim, dimAux>* getTrackPointKalmanFitterInfo(const TrackPoint* tp, const AbsTrackRep* rep) {
    return dynamic_cast<KalmanFitterInfo<dim, dimAux>*>(tp->getFitterInfo(rep));
  }

}