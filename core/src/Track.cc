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

#include "Track.h"

#include "Exception.h"
#include "IO.h"
#include "PlanarMeasurement.h"
#include "AbsMeasurement.h"

#include "WireTrackCandHit.h"

#include <algorithm>
#include <map>

#include <TDatabasePDG.h>
#include <TMath.h>

//#include <glog/logging.h>

//#define DEBUG


namespace genfit {

template<unsigned int dimMeas>
Track<dimMeas>::Track() :
  cardinalRep_(0), fitStatuses_(), mcTrackId_(-1), timeSeed_(0), stateSeed_(6), covSeed_(6)
{
  ;
}


template<unsigned int dimMeas>
Track<dimMeas>::Track(const TrackCand& trackCand, const MeasurementFactory<AbsMeasurement<dimMeas>>& factory, AbsTrackRep* rep) :
  cardinalRep_(0), fitStatuses_(), stateSeed_(SVector6()), covSeed_(SMatrixSym6())
{

  if (rep != nullptr)
    addTrackRep(rep);

  createMeasurements(trackCand, factory);

  // Copy seed information from candidate
  timeSeed_ = trackCand.getTimeSeed();
  stateSeed_ = trackCand.getStateSeed();
  covSeed_ = trackCand.getCovSeed();

  mcTrackId_ = trackCand.getMcTrackId();

  // fill cache
  fillPointsWithMeasurement();

  checkConsistency();
}

template<unsigned int dimMeas>
void
Track<dimMeas>::createMeasurements(const TrackCand& trackCand, const MeasurementFactory<AbsMeasurement<dimMeas>>& factory)
{
  // create the measurements using the factory.
  const std::vector <AbsMeasurement<dimMeas>*>& factoryHits = factory.createMany(trackCand);

  if (factoryHits.size() != trackCand.getNHits()) {
    Exception exc("Track::Track ==> factoryHits.size() != trackCand->getNHits()",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  // create TrackPoints
  for (unsigned int i=0; i<factoryHits.size(); ++i){
    TrackPoint<dimMeas>* tp = new TrackPoint(factoryHits[i], this);
    tp->setSortingParameter(trackCand.getHit(i)->getSortingParameter());
    insertPoint(tp);
  }
}


template<unsigned int dimMeas>
Track<dimMeas>::Track(AbsTrackRep* trackRep, const SVector6& stateSeed) :
  cardinalRep_(0), fitStatuses_(), mcTrackId_(-1), timeSeed_(0), stateSeed_(stateSeed),
  covSeed_(ROOT::Math::SMatrixIdentity())
{
  addTrackRep(trackRep);
}


template<unsigned int dimMeas>
Track<dimMeas>::Track(AbsTrackRep* trackRep, const ROOT::Math::XYZVector& posSeed, const ROOT::Math::XYZVector& momSeed) :
  cardinalRep_(0), fitStatuses_(), mcTrackId_(-1), timeSeed_(0), stateSeed_(SVector6()),
  covSeed_(ROOT::Math::SMatrixIdentity())
{
  addTrackRep(trackRep);
  setStateSeed(posSeed, momSeed);
}


template<unsigned int dimMeas>
Track<dimMeas>::Track(AbsTrackRep* trackRep, const SVector6& stateSeed, const SMatrixSym6& covSeed) :
  cardinalRep_(0), fitStatuses_(), mcTrackId_(-1), timeSeed_(0), stateSeed_(stateSeed),
  covSeed_(covSeed)
{
  addTrackRep(trackRep);
}


template<unsigned int dimMeas>
Track<dimMeas>::Track(const Track& rhs) :
  cardinalRep_(rhs.cardinalRep_), mcTrackId_(rhs.mcTrackId_), timeSeed_(rhs.timeSeed_),
  stateSeed_(rhs.stateSeed_), covSeed_(rhs.covSeed_)
{
  rhs.checkConsistency();

  std::map<const AbsTrackRep*, AbsTrackRep*> oldRepNewRep;

  for (std::vector<AbsTrackRep*>::const_iterator it=rhs.trackReps_.begin(); it!=rhs.trackReps_.end(); ++it) {
    AbsTrackRep* newRep = (*it)->clone();
    addTrackRep(newRep);
    oldRepNewRep[(*it)] = newRep;
  }

  trackPoints_.reserve(rhs.trackPoints_.size());
  for (std::vector<TrackPoint<dimMeas>*>::const_iterator it=rhs.trackPoints_.begin(); it!=rhs.trackPoints_.end(); ++it) {
    trackPoints_.push_back(new TrackPoint(**it, oldRepNewRep));
    trackPoints_.back()->setTrack(this);
  }

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=rhs.fitStatuses_.begin(); it!=rhs.fitStatuses_.end(); ++it) {
    setFitStatus(it->second->clone(), oldRepNewRep[it->first]);
  }

  fillPointsWithMeasurement();

  checkConsistency();
}

template<unsigned int dimMeas>
Track& Track<dimMeas>::operator=(Track other) {
  swap(other);

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator it=trackPoints_.begin(); it!=trackPoints_.end(); ++it) {
    (*it)->setTrack(this);
  }

  fillPointsWithMeasurement();

  checkConsistency();

  return *this;
}

template<unsigned int dimMeas>
void Track<dimMeas>::swap(Track& other) {
  std::swap(this->trackReps_, other.trackReps_);
  std::swap(this->cardinalRep_, other.cardinalRep_);
  std::swap(this->trackPoints_, other.trackPoints_);
  std::swap(this->trackPointsWithMeasurement_, other.trackPointsWithMeasurement_);
  std::swap(this->fitStatuses_, other.fitStatuses_);
  std::swap(this->mcTrackId_, other.mcTrackId_);
  std::swap(this->timeSeed_, other.timeSeed_);
  std::swap(this->stateSeed_, other.stateSeed_);
  std::swap(this->covSeed_, other.covSeed_);

}

template<unsigned int dimMeas>
Track<dimMeas>::~Track() {
  this->Clear();
}

template<unsigned int dimMeas>
void Track<dimMeas>::Clear(Option_t*)
{
  // This function is needed for TClonesArray embedding.
  // FIXME: smarter containers or pointers needed ...
  for (size_t i = 0; i < trackPoints_.size(); ++i)
    delete trackPoints_[i];

  trackPoints_.clear();
  trackPointsWithMeasurement_.clear();

  for (std::map< const AbsTrackRep*, FitStatus* >::iterator it = fitStatuses_.begin(); it!= fitStatuses_.end(); ++it)
    delete it->second;
  fitStatuses_.clear();

  for (size_t i = 0; i < trackReps_.size(); ++i)
    delete trackReps_[i];
  trackReps_.clear();

  cardinalRep_ = 0;

  mcTrackId_ = -1;

  timeSeed_ = 0;
  stateSeed_ = SVector6();
  covSeed_ = SMatrixSym6();
}


template<unsigned int dimMeas>
TrackPoint<dimMeas>* Track<dimMeas>::getPoint(int id) const {
  if (id < 0)
    id += trackPoints_.size();

  return trackPoints_.at(id);
}


template<unsigned int dimMeas>
TrackPoint<dimMeas>* Track<dimMeas>::getPointWithMeasurement(int id) const {
  if (id < 0)
    id += trackPointsWithMeasurement_.size();

  return trackPointsWithMeasurement_.at(id);
}


template<unsigned int dimMeas>
TrackPoint<dimMeas>* Track<dimMeas>::getPointWithMeasurementAndFitterInfo(int id, const AbsTrackRep* rep) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  if (id >= 0) {
    int i = 0;
    for (std::vector<TrackPoint<dimMeas>*>::const_iterator it = trackPointsWithMeasurement_.begin(); it != trackPointsWithMeasurement_.end(); ++it) {
      if ((*it)->hasFitterInfo(rep)) {
        if (id == i)
          return (*it);
        ++i;
      }
    }
  } else {
    // Search backwards.
    int i = -1;
    for (std::vector<TrackPoint<dimMeas>*>::const_reverse_iterator it = trackPointsWithMeasurement_.rbegin(); it != trackPointsWithMeasurement_.rend(); ++it) {
      if ((*it)->hasFitterInfo(rep)) {
        if (id == i)
          return (*it);
        --i;
      }
    }
  }

  // Not found, i.e. abs(id) > number of fitted TrackPoints
  return 0;
}


template<unsigned int dimMeas>
TrackPoint<dimMeas>* Track<dimMeas>::getPointWithFitterInfo(int id, const AbsTrackRep* rep) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  if (id >= 0) {
    int i = 0;
    for (std::vector<TrackPoint<dimMeas>*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
      if ((*it)->hasFitterInfo(rep)) {
        if (id == i)
          return (*it);
        ++i;
      }
    }
  } else {
    // Search backwards.
    int i = -1;
    for (std::vector<TrackPoint<dimMeas>*>::const_reverse_iterator it = trackPoints_.rbegin(); it != trackPoints_.rend(); ++it) {
      if ((*it)->hasFitterInfo(rep)) {
        if (id == i)
          return (*it);
        --i;
      }
    }
  }

  // Not found, i.e. abs(id) > number of fitted TrackPoints
  return 0;
}


template<unsigned int dimMeas>
const MeasuredStateOnPlane& Track<dimMeas>::getFittedState(int id, const AbsTrackRep* rep, bool biased) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  TrackPoint<dimMeas>* point = getPointWithFitterInfo(id, rep);
  if (point == nullptr) {
    Exception exc("Track::getFittedState ==> no trackPoint with fitterInfo for rep",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  return point->getFitterInfo(rep)->getFittedState(biased);
}


template<unsigned int dimMeas>
int Track<dimMeas>::getIdForRep(const AbsTrackRep* rep) const
{
  for (size_t i = 0; i < trackReps_.size(); ++i)
    if (trackReps_[i] == rep)
      return i;

  Exception exc("Track::getIdForRep ==> cannot find TrackRep in Track",__LINE__,__FILE__);
  exc.setFatal();
  throw exc;
}


template<unsigned int dimMeas>
bool Track<dimMeas>::hasFitStatus(const AbsTrackRep* rep) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  if (fitStatuses_.find(rep) == fitStatuses_.end())
    return false;

  return (fitStatuses_.at(rep) != nullptr);
}


template<unsigned int dimMeas>
void Track<dimMeas>::setFitStatus(FitStatus* fitStatus, const AbsTrackRep* rep) {
  if (fitStatuses_.find(rep) != fitStatuses_.end())
    delete fitStatuses_.at(rep);

  fitStatuses_[rep] = fitStatus;
}


template<unsigned int dimMeas>
void Track<dimMeas>::setStateSeed(const ROOT::Math::XYZVector& pos, const ROOT::Math::XYZVector& mom) {
  stateSeed_(0) = pos.X();
  stateSeed_(1) = pos.Y();
  stateSeed_(2) = pos.Z();

  stateSeed_(3) = mom.X();
  stateSeed_(4) = mom.Y();
  stateSeed_(5) = mom.Z();
}



template<unsigned int dimMeas>
void Track<dimMeas>::insertPoint(TrackPoint<dimMeas>* point, int id) {

  point->setTrack(this);

  #ifdef DEBUG
  debugOut << "Track::insertPoint at position " << id  << "\n";
  #endif
  assert(point!=nullptr);
  trackHasChanged();

  point->setTrack(this);

  if (trackPoints_.size() == 0) {
    trackPoints_.push_back(point);

    if (point->hasRawMeasurements())
      trackPointsWithMeasurement_.push_back(point);

    return;
  }

  if (id == -1 || id == (int)trackPoints_.size()) {
    trackPoints_.push_back(point);

    if (point->hasRawMeasurements())
      trackPointsWithMeasurement_.push_back(point);

    deleteReferenceInfo(std::max(0, (int)trackPoints_.size()-2), (int)trackPoints_.size()-1);

    // delete fitter infos if inserted point has a measurement
    if (point->hasRawMeasurements()) {
      deleteForwardInfo(-1, -1);
      deleteBackwardInfo(0, -2);
    }

    return;
  }

  // [-size, size-1] is the allowed range
  assert(id < (ssize_t)trackPoints_.size() || -id-1 <= (ssize_t)trackPoints_.size());

  if (id < 0)
    id += trackPoints_.size() + 1;

  // insert
  trackPoints_.insert(trackPoints_.begin() + id, point);  // insert inserts BEFORE

  // delete fitter infos if inserted point has a measurement
  if (point->hasRawMeasurements()) {
    deleteForwardInfo(id, -1);
    deleteBackwardInfo(0, id);
  }

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id+1));

  fillPointsWithMeasurement();
}


template<unsigned int dimMeas>
void Track<dimMeas>::insertPoints(std::vector<TrackPoint<dimMeas>*> points, int id) {

  int nBefore = getNumPoints();
  int n = points.size();

  if (n == 0)
    return;
  if (n == 1) {
    insertPoint(points[0], id);
    return;
  }

  for (std::vector<TrackPoint<dimMeas>*>::iterator p = points.begin(); p != points.end(); ++p)
    (*p)->setTrack(this);

  if (id == -1 || id == (int)trackPoints_.size()) {
    trackPoints_.insert(trackPoints_.end(), points.begin(), points.end());

    deleteReferenceInfo(std::max(0, nBefore-1), nBefore);

    deleteForwardInfo(nBefore, -1);
    deleteBackwardInfo(0, std::max(0, nBefore-1));

    fillPointsWithMeasurement();

    return;
  }


  assert(id < (ssize_t)trackPoints_.size() || -id-1 <= (ssize_t)trackPoints_.size());

  if (id < 0)
    id += trackPoints_.size() + 1;


  // insert
  trackPoints_.insert(trackPoints_.begin() + id, points.begin(), points.end());  // insert inserts BEFORE

  // delete fitter infos if inserted point has a measurement
  deleteForwardInfo(id, -1);
  deleteBackwardInfo(0, id+n);

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id));
  deleteReferenceInfo(std::max(0, id+n-1), std::min((int)trackPoints_.size()-1, id+n));

  fillPointsWithMeasurement();
}


template<unsigned int dimMeas>
void Track<dimMeas>::deletePoint(int id) {

  #ifdef DEBUG
  debugOut << "Track::deletePoint at position " << id  << "\n";
  #endif

  trackHasChanged();

  if (id < 0)
    id += trackPoints_.size();
  assert(id>0);


  // delete forwardInfo after point (backwardInfo before point) if deleted point has a measurement
  if (trackPoints_[id]->hasRawMeasurements()) {
    deleteForwardInfo(id, -1);
    deleteBackwardInfo(0, id-1);
  }

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id+1));

  // delete point
  std::vector<TrackPoint<dimMeas>*>::iterator it = std::find(trackPointsWithMeasurement_.begin(), trackPointsWithMeasurement_.end(), trackPoints_[id]);
  if (it != trackPointsWithMeasurement_.end())
    trackPointsWithMeasurement_.erase(it);

  delete trackPoints_[id];
  trackPoints_.erase(trackPoints_.begin()+id);

  fillPointsWithMeasurement();

}


template<unsigned int dimMeas>
void Track<dimMeas>::insertMeasurement(AbsMeasurement<dimMeas>* measurement, int id) {
  insertPoint(new TrackPoint(measurement, this), id);
}
  
void Track<dimMeas>::deleteFittedState(const genfit::AbsTrackRep* rep) {
  if(hasFitStatus(rep)) {
    delete fitStatuses_.at(rep);
    fitStatuses_.erase(rep);
  }

  // delete FitterInfos related to the deleted TrackRep
  for (const auto& trackPoint : trackPoints_) {
    if(trackPoint->hasFitterInfo(rep)) {
      trackPoint->deleteFitterInfo(rep);
    }
  }
}


template<unsigned int dimMeas>
void Track<dimMeas>::mergeTrack(const Track* other, int id) {

  #ifdef DEBUG
  debugOut << "Track::mergeTrack\n";
  #endif

  if (other->getNumPoints() == 0)
    return;

  std::map<const AbsTrackRep*, AbsTrackRep*> otherRepThisRep;
  std::vector<const AbsTrackRep*> otherRepsToRemove;
  otherRepsToRemove.reserve(other->trackReps_.size());

  for (std::vector<AbsTrackRep*>::const_iterator otherRep=other->trackReps_.begin(); otherRep!=other->trackReps_.end(); ++otherRep) {
    bool found(false);
    for (std::vector<AbsTrackRep*>::const_iterator thisRep=trackReps_.begin(); thisRep!=trackReps_.end(); ++thisRep) {
      if ((*thisRep)->isSame(*otherRep)) {
        otherRepThisRep[*otherRep] = *thisRep;
        #ifdef DEBUG
        debugOut << " map other rep " << *otherRep << " to " << (*thisRep) << "\n";
        #endif
        if (found) {
          Exception exc("Track::mergeTrack ==> more than one matching rep.",__LINE__,__FILE__);
          exc.setFatal();
          throw exc;
        }
        found = true;
        break;
      }
    }
    if (!found) {
      otherRepsToRemove.push_back(*otherRep);
      #ifdef DEBUG
      debugOut << " remove other rep " << *otherRep << "\n";
      #endif
    }
  }


  std::vector<TrackPoint<dimMeas>*> points;
  points.reserve(other->getNumPoints());

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator otherTp=other->trackPoints_.begin(); otherTp!=other->trackPoints_.end(); ++otherTp) {
    points.push_back(new TrackPoint(**otherTp, otherRepThisRep, &otherRepsToRemove));
  }

  insertPoints(points, id);
}


template<unsigned int dimMeas>
void Track<dimMeas>::addTrackRep(AbsTrackRep* trackRep) {
  trackReps_.push_back(trackRep);
  fitStatuses_[trackRep] = new FitStatus();
}


template<unsigned int dimMeas>
void Track<dimMeas>::deleteTrackRep(int id) {
  if (id < 0)
    id += trackReps_.size();

  AbsTrackRep* rep = trackReps_.at(id);

  // update cardinalRep_
  if (int(cardinalRep_) == id)
    cardinalRep_ = 0; // reset
  else if (int(cardinalRep_) > id)
    --cardinalRep_; // make cardinalRep_ point to the same TrackRep before and after deletion

  // delete FitterInfos related to the deleted TrackRep
  for (std::vector<TrackPoint<dimMeas>*>::const_iterator pointIt = trackPoints_.begin(); pointIt != trackPoints_.end(); ++pointIt) {
    (*pointIt)->deleteFitterInfo(rep);
  }

  // delete fitStatus
  delete fitStatuses_.at(rep);
  fitStatuses_.erase(rep);

  // delete rep
  delete rep;
  trackReps_.erase(trackReps_.begin()+id);
}


template<unsigned int dimMeas>
void Track<dimMeas>::setCardinalRep(int id) {

  if (id < 0)
    id += trackReps_.size();

  if (id >= 0 && (unsigned int)id < trackReps_.size())
    cardinalRep_ = id;
  else {
    cardinalRep_ = 0;
    errorOut << "Track::setCardinalRep: Attempted to set cardinalRep_ to a value out of bounds. Resetting  cardinalRep_ to 0." << std::endl;
  }
}


template<unsigned int dimMeas>
void Track<dimMeas>::determineCardinalRep() {

  // Todo: test

  if (trackReps_.size() <= 1)
    return;

  double minChi2(9.E99);
  const AbsTrackRep* bestRep(nullptr);

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it = fitStatuses_.begin(); it != fitStatuses_.end(); ++it) {
    if (it->second->isFitConverged()) {
      if (it->second->getChi2() < minChi2) {
        minChi2 = it->second->getChi2();
        bestRep = it->first;
      }
    }
  }

  if (bestRep != nullptr) {
    setCardinalRep(getIdForRep(bestRep));
  }
}


template<unsigned int dimMeas>
bool Track<dimMeas>::sort() {
  #ifdef DEBUG
  debugOut << "Track::sort \n";
  #endif

  int nPoints(trackPoints_.size());
  // original order
  const std::vector<TrackPoint<dimMeas>*> pointsBefore(trackPoints_);

  // sort
  std::stable_sort(trackPoints_.begin(), trackPoints_.end(), TrackPointComparator());

  // see where order changed
  int equalUntil(-1), equalFrom(nPoints);
  for (int i = 0; i<nPoints; ++i) {
    if (pointsBefore[i] == trackPoints_[i])
      equalUntil = i;
    else
      break;
  }

  if (equalUntil == nPoints-1)
    return false; // sorting did not change anything


  trackHasChanged();

  for (int i = nPoints-1; i>equalUntil; --i) {
    if (pointsBefore[i] == trackPoints_[i])
      equalFrom = i;
    else
      break;
  }

  #ifdef DEBUG
  debugOut << "Track::sort. Equal up to (including) hit " << equalUntil << " and from (including) hit " << equalFrom << " \n";
  #endif

  deleteForwardInfo(equalUntil+1, -1);
  deleteBackwardInfo(0, equalFrom-1);

  deleteReferenceInfo(std::max(0, equalUntil+1), std::min((int)trackPoints_.size()-1, equalFrom-1));

  fillPointsWithMeasurement();

  return true;
}


template<unsigned int dimMeas>
bool Track<dimMeas>::udpateSeed(int id, AbsTrackRep* rep, bool biased) {
  try {
    const MeasuredStateOnPlane& fittedState = getFittedState(id, rep, biased);
    setTimeSeed(fittedState.getTime());
    setStateSeed(fittedState.get6DState());
    setCovSeed(fittedState.get6DCov());

    double fittedCharge = fittedState.getCharge();

    for (unsigned int i = 0; i<trackReps_.size(); ++i) {
      if (trackReps_[i]->getPDGCharge() * fittedCharge < 0) {
        trackReps_[i]->switchPDGSign();
      }
    }
  }
  catch (Exception& e) {
    // in this case the original track seed will be used
    return false;
  }
  return true;
}


template<unsigned int dimMeas>
void Track<dimMeas>::reverseTrackPoints() {

  std::reverse(trackPoints_.begin(),trackPoints_.end());

  deleteForwardInfo(0, -1);
  deleteBackwardInfo(0, -1);
  deleteReferenceInfo(0, -1);

  fillPointsWithMeasurement();
}


template<unsigned int dimMeas>
void Track<dimMeas>::switchPDGSigns(AbsTrackRep* rep) {
  if (rep != nullptr) {
    rep->switchPDGSign();
    return;
  }

  for (unsigned int i = 0; i<trackReps_.size(); ++i) {
    trackReps_[i]->switchPDGSign();
  }
}


template<unsigned int dimMeas>
void Track<dimMeas>::reverseTrack() {
  udpateSeed(-1); // set fitted state of last hit as new seed
  reverseMomSeed(); // flip momentum direction
  switchPDGSigns();
  reverseTrackPoints(); // also deletes all fitterInfos
}


template<unsigned int dimMeas>
void Track<dimMeas>::deleteForwardInfo(int startId, int endId, const AbsTrackRep* rep) {
  #ifdef DEBUG
  debugOut << "Track::deleteForwardInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteForwardInfo();
    }
    else {
      const std::vector<AbsFitterInfo*>& fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteForwardInfo();
      }
    }
  }
}

template<unsigned int dimMeas>
void Track<dimMeas>::deleteBackwardInfo(int startId, int endId, const AbsTrackRep* rep) {

  #ifdef DEBUG
  debugOut << "Track::deleteBackwardInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);


  for (std::vector<TrackPoint<dimMeas>*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteBackwardInfo();
    }
    else {
      const std::vector<AbsFitterInfo*>& fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteBackwardInfo();
      }
    }
  }
}

template<unsigned int dimMeas>
void Track<dimMeas>::deleteReferenceInfo(int startId, int endId, const AbsTrackRep* rep) {

  #ifdef DEBUG
  debugOut << "Track::deleteReferenceInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteReferenceInfo();
    }
    else {
      const std::vector<AbsFitterInfo*>& fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteReferenceInfo();
      }
    }
  }
}

template<unsigned int dimMeas>
void Track<dimMeas>::deleteMeasurementInfo(int startId, int endId, const AbsTrackRep* rep) {

  #ifdef DEBUG
  debugOut << "Track::deleteMeasurementInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteMeasurementInfo();
    }
    else {
      const std::vector<AbsFitterInfo*>& fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteMeasurementInfo();
      }
    }
  }
}

template<unsigned int dimMeas>
void Track<dimMeas>::deleteFitterInfo(int startId, int endId, const AbsTrackRep* rep) {

  #ifdef DEBUG
  debugOut << "Track::deleteFitterInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->deleteFitterInfo(rep);
    }
    else {
      for (std::vector<AbsTrackRep*>::const_iterator repIt = trackReps_.begin(); repIt != trackReps_.end(); ++repIt) {
        if ((*pointIt)->hasFitterInfo(*repIt))
          (*pointIt)->deleteFitterInfo(*repIt);
      }
    }
  }
}


template<unsigned int dimMeas>
double Track<dimMeas>::getTrackLen(AbsTrackRep* rep, int startId, int endId) const {

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();

  bool backwards(false);
  if (startId > endId) {
    double temp = startId;
    startId = endId;
    endId = temp;
    backwards = true;
  }

  endId += 1;

  if (rep == nullptr)
    rep = getCardinalRep();

  double trackLen(0);
  StateOnPlane state;

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (! (*pointIt)->hasFitterInfo(rep)) {
      Exception e("Track::getTracklength: trackPoint has no fitterInfo", __LINE__,__FILE__);
      throw e;
    }

    if (pointIt != trackPoints_.begin() + startId) {
      trackLen += rep->extrapolateToPlane(state, (*pointIt)->getFitterInfo(rep)->getPlane());
    }

    state = (*pointIt)->getFitterInfo(rep)->getFittedState();
  }

  if (backwards)
    trackLen *= -1.;

  return trackLen;
}


template<unsigned int dimMeas>
TrackCand* Track<dimMeas>::constructTrackCand() const {
  TrackCand* cand = new TrackCand();

  cand->setTime6DSeedAndPdgCode(timeSeed_, stateSeed_, getCardinalRep()->getPDG());
  cand->setCovSeed(covSeed_);
  cand->setMcTrackId(mcTrackId_);

  for (unsigned int i = 0; i < trackPointsWithMeasurement_.size(); ++i) {
    const TrackPoint<dimMeas>* tp = trackPointsWithMeasurement_[i];
    const std::vector< AbsMeasurement<dimMeas>* >& measurements = tp->getRawMeasurements();

    for (unsigned int j = 0; j < measurements.size(); ++j) {
      const AbsMeasurement<dimMeas>* m = measurements[j];
      TrackCandHit* tch;

      int planeId = -1;
      if (dynamic_cast<const PlanarMeasurement*>(m)) {
        planeId = static_cast<const PlanarMeasurement*>(m)->getPlaneId();
      }

      if (m->isLeftRightMeasurement()) {
        tch = new WireTrackCandHit(m->getDetId(), m->getHitId(), planeId,
            tp->getSortingParameter(), m->getLeftRightResolution());
      }
      else {
        tch = new TrackCandHit(m->getDetId(), m->getHitId(), planeId,
            tp->getSortingParameter());
      }
      cand->addHit(tch);
    }
  }

  return cand;
}


template<unsigned int dimMeas>
double Track<dimMeas>::getTOF(AbsTrackRep* rep, int startId, int endId) const {

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();

  if (startId > endId) {
    std::swap(startId, endId);
  }

  endId += 1;

  if (rep == nullptr)
    rep = getCardinalRep();

  StateOnPlane state;

  const TrackPoint<dimMeas>* startPoint(trackPoints_[startId]);
  const TrackPoint<dimMeas>* endPoint(trackPoints_[endId]);
  
  if (!startPoint->hasFitterInfo(rep)
      || !endPoint->hasFitterInfo(rep)) {
      Exception e("Track::getTOF: trackPoint has no fitterInfo", __LINE__,__FILE__);
      throw e;
    }

  double tof = (rep->getTime(endPoint->getFitterInfo(rep)->getFittedState())
		- rep->getTime(startPoint->getFitterInfo(rep)->getFittedState()));
  return tof;
}

template<unsigned int dimMeas>
void Track<dimMeas>::prune(const Option_t* option) {

  PruneFlags f;
  f.setFlags(option);

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=fitStatuses_.begin(); it!=fitStatuses_.end(); ++it) {
    it->second->getPruneFlags().setFlags(option);
  }

  // prune trackPoints
  if (f.hasFlags("F") || f.hasFlags("L")) {
    TrackPoint<dimMeas>* firstPoint = getPointWithFitterInfo(0);
    TrackPoint<dimMeas>* lastPoint = getPointWithFitterInfo(-1);
    for (unsigned int i = 0; i<trackPoints_.size(); ++i) {
      if (trackPoints_[i] == firstPoint && f.hasFlags("F"))
        continue;

      if (trackPoints_[i] == lastPoint && f.hasFlags("L"))
        continue;

      delete trackPoints_[i];
      trackPoints_.erase(trackPoints_.begin()+i);
      --i;
    }
  }

  // prune TrackReps
  if (f.hasFlags("C")) {
    for (unsigned int i = 0; i < trackReps_.size(); ++i) {
      if (i != cardinalRep_) {
        deleteTrackRep(i);
        --i;
      }
    }
  }


  // from remaining trackPoints: prune measurementsOnPlane, unneeded fitterInfoStuff
  for (unsigned int i = 0; i<trackPoints_.size(); ++i) {
    if (f.hasFlags("W"))
      trackPoints_[i]->deleteRawMeasurements();

    std::vector< AbsFitterInfo* > fis =  trackPoints_[i]->getFitterInfos();
    for (unsigned int j = 0; j<fis.size(); ++j) {

      if (i == 0 && f.hasFlags("FLI"))
        fis[j]->deleteForwardInfo();
      else if (i == trackPoints_.size()-1 && f.hasFlags("FLI"))
        fis[j]->deleteBackwardInfo();
      else if (f.hasFlags("FI"))
        fis[j]->deleteForwardInfo();
      else if (f.hasFlags("LI"))
        fis[j]->deleteBackwardInfo();

      // also delete reference info if points have been removed since it is invalid then!
      if (f.hasFlags("R") or f.hasFlags("F") or f.hasFlags("L"))
        fis[j]->deleteReferenceInfo();
      if (f.hasFlags("M"))
        fis[j]->deleteMeasurementInfo();
    }
  }

  fillPointsWithMeasurement();

  #ifdef DEBUG
  debugOut << "pruned Track: "; Print();
  #endif

}


template<unsigned int dimMeas>
void Track<dimMeas>::Print(const Option_t* option) const {
  TString opt = option;
  opt.ToUpper();
  if (opt.Contains("C")) { // compact

    printOut << "\n    ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {

      int color = 32*(size_t)(trackPoints_[i]) % 15;
      switch (color) {
        case 0:
          printOut<<"\033[1;30m";
          break;
        case 1:
          printOut<<"\033[0;34m";
          break;
        case 2:
          printOut<<"\033[1;34m";
          break;
        case 3:
          printOut<<"\033[0;32m";
          break;
        case 4:
          printOut<<"\033[1;32m";
          break;
        case 5:
          printOut<<"\033[0;36m";
          break;
        case 6:
          printOut<<"\033[1;36m";
          break;
        case 7:
          printOut<<"\033[0;31m";
          break;
        case 8:
          printOut<<"\033[1;31m";
          break;
        case 9:
          printOut<<"\033[0;35m";
          break;
        case 10:
          printOut<<"\033[1;35m";
          break;
        case 11:
          printOut<<"\033[0;33m";
          break;
        case 12:
          printOut<<"\033[1;33m";
          break;
        case 13:
          printOut<<"\033[0;37m";
          break;
        default:
          ;
      }
      printOut << trackPoints_[i] << "\033[00m  ";
    }
    printOut << "\n";

    printOut << "   ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {
      printf("% -9.3g  ", trackPoints_[i]->getSortingParameter());
    }

    for (std::vector<AbsTrackRep*>::const_iterator rep = trackReps_.begin(); rep != trackReps_.end(); ++rep) {
      printOut << "\n" << getIdForRep(*rep) << "   ";
      for (unsigned int i=0; i<trackPoints_.size(); ++i) {
        if (! trackPoints_[i]->hasFitterInfo(*rep)) {
          printOut << "           ";
          continue;
        }
        AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(*rep);
        if (fi->hasMeasurements())
          printOut << "M";
        else
          printOut << " ";

        if (fi->hasReferenceState())
          printOut << "R";
        else
          printOut << " ";

        printOut << "         ";
      }
      printOut << "\n";

      printOut << " -> ";
      for (unsigned int i=0; i<trackPoints_.size(); ++i) {
        if (! trackPoints_[i]->hasFitterInfo(*rep)) {
          printOut << "           ";
          continue;
        }
        AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(*rep);
        if (fi->hasForwardPrediction())
          printOut << "P";
        else
          printOut << " ";

        if (fi->hasForwardUpdate())
          printOut << "U";
        else
          printOut << " ";

        printOut << "         ";
      }
      printOut << "\n";

      printOut << " <- ";
      for (unsigned int i=0; i<trackPoints_.size(); ++i) {
        if (! trackPoints_[i]->hasFitterInfo(*rep)) {
          printOut << "           ";
          continue;
        }
        AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(*rep);
        if (fi->hasBackwardPrediction())
          printOut << "P";
        else
          printOut << " ";

        if (fi->hasBackwardUpdate())
          printOut << "U";
        else
          printOut << " ";

        printOut << "         ";
      }

      printOut << "\n";

    } //end loop over reps

    printOut << "\n";
    return;
  }



  printOut << "=======================================================================================\n";
  printOut << "genfit::Track, containing " << trackPoints_.size() << " TrackPoints and " << trackReps_.size() << " TrackReps.\n";
  printOut << " Seed state: "; stateSeed_.Print(printOut);

  for (unsigned int i=0; i<trackReps_.size(); ++i) {
    printOut << " TrackRep Nr. " << i;
    if (i == cardinalRep_)
      printOut << " (This is the cardinal rep)";
    printOut << "\n";
    trackReps_[i]->Print();
  }

  printOut << "---------------------------------------------------------------------------------------\n";

  for (unsigned int i=0; i<trackPoints_.size(); ++i) {
    printOut << "TrackPoint Nr. " << i << "\n";
    trackPoints_[i]->Print();
    printOut << "..........................................................................\n";
  }

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=fitStatuses_.begin(); it!=fitStatuses_.end(); ++it) {
    it->second->Print();
  }

  printOut << "=======================================================================================\n";

}


template<unsigned int dimMeas>
void Track<dimMeas>::checkConsistency() const {

  // cppcheck-suppress unreadVariable
  bool consistent = true;
  std::stringstream failures;

  if (*(std::max_element(covSeed_.begin(), covSeed_.end())) == 0.) {
    // Nota bene: The consistency is not set to false when this occurs, because it does not break the consistency of
    // the track. However, when something else fails we keep this as additional error information.
    failures << "Track::checkConsistency(): Warning: covSeed_ is zero" << std::endl;
  }

  // check if correct number of fitStatuses
  if (fitStatuses_.size() != trackReps_.size()) {
    failures << "Track::checkConsistency(): Number of fitStatuses is != number of TrackReps " << std::endl;
    // cppcheck-suppress unreadVariable
    consistent = false;
  }

  // check if cardinalRep_ is in range of trackReps_
  if (trackReps_.size() && cardinalRep_ >= trackReps_.size()) {
    failures << "Track::checkConsistency(): cardinalRep id " << cardinalRep_ << " out of bounds" << std::endl;
    // cppcheck-suppress unreadVariable
    consistent = false;
  }

  for (std::vector<AbsTrackRep*>::const_iterator rep = trackReps_.begin(); rep != trackReps_.end(); ++rep) {
    // check for nullptr
    if ((*rep) == nullptr) {
      failures << "Track::checkConsistency(): TrackRep is nullptr" << std::endl;
      // cppcheck-suppress unreadVariable
      consistent = false;
    }

    // check for valid pdg code
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle((*rep)->getPDG());
    if (particle == nullptr) {
      failures << "Track::checkConsistency(): TrackRep pdg ID " << (*rep)->getPDG() << " is not valid" << std::endl;
      // cppcheck-suppress unreadVariable
      consistent = false;
    }

    // check if corresponding FitStatus is there
    if (fitStatuses_.find(*rep) == fitStatuses_.end() and fitStatuses_.find(*rep)->second != nullptr) {
      failures << "Track::checkConsistency(): No FitStatus for Rep or FitStatus is nullptr" << std::endl;
      // cppcheck-suppress unreadVariable
      consistent = false;
    }
  }

  // check TrackPoints
  for (std::vector<TrackPoint<dimMeas>*>::const_iterator tp = trackPoints_.begin(); tp != trackPoints_.end(); ++tp) {
    // check for nullptr
    if ((*tp) == nullptr) {
      failures << "Track::checkConsistency(): TrackPoint is nullptr" << std::endl;
      // cppcheck-suppress unreadVariable
      consistent = false;
    }
    // check if trackPoint points back to this track
    if ((*tp)->getTrack() != this) {
      failures << "Track::checkConsistency(): TrackPoint does not point back to this track" << std::endl;
      // cppcheck-suppress unreadVariable
      consistent = false;
    }

    // check rawMeasurements
    const std::vector<AbsMeasurement<dimMeas>*>& rawMeasurements = (*tp)->getRawMeasurements();
    for (std::vector<AbsMeasurement<dimMeas>*>::const_iterator m = rawMeasurements.begin(); m != rawMeasurements.end(); ++m) {
      // check for nullptr
      if ((*m) == nullptr) {
        failures << "Track::checkConsistency(): Measurement is nullptr" << std::endl;
	// cppcheck-suppress unreadVariable
        consistent = false;
      }
      // check if measurement points back to TrackPoint
      if ((*m)->getTrackPoint() != *tp) {
        failures << "Track::checkConsistency(): Measurement does not point back to correct TrackPoint" << std::endl;
	// cppcheck-suppress unreadVariable
        consistent = false;
      }
    }

    // check fitterInfos
    const std::vector<AbsFitterInfo*>& fitterInfos = (*tp)->getFitterInfos();
    for (std::vector<AbsFitterInfo*>::const_iterator fi = fitterInfos.begin(); fi != fitterInfos.end(); ++fi) {
      // check for nullptr
      if ((*fi) == nullptr) {
        failures << "Track::checkConsistency(): FitterInfo is nullptr. TrackPoint: " << *tp << std::endl;
	// cppcheck-suppress unreadVariable
        consistent = false;
      }

      // check if fitterInfos point to valid TrackReps in trackReps_
      int mycount (0);
      for (std::vector<AbsTrackRep*>::const_iterator rep = trackReps_.begin(); rep != trackReps_.end(); ++rep) {
        if ( (*rep) == (*fi)->getRep() ) {
          ++mycount;
        }
      }
      if (mycount ==  0) {
        failures << "Track::checkConsistency(): fitterInfo points to TrackRep which is not in Track" << std::endl;
	// cppcheck-suppress unreadVariable
        consistent = false;
      }

      if (!( (*fi)->checkConsistency(&(this->getFitStatus((*fi)->getRep())->getPruneFlags())) ) ) {
        failures << "Track::checkConsistency(): FitterInfo not consistent. TrackPoint: " << *tp << std::endl;
	// cppcheck-suppress unreadVariable
        consistent = false;
      }

    } // end loop over FitterInfos

  } // end loop over TrackPoints


  // check trackPointsWithMeasurement_
  std::vector<TrackPoint<dimMeas>*> trackPointsWithMeasurement;
  trackPointsWithMeasurement.reserve(trackPoints_.size());

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
    if ((*it)->hasRawMeasurements()) {
      trackPointsWithMeasurement.push_back(*it);
    }
  }

  if (trackPointsWithMeasurement.size() != trackPointsWithMeasurement_.size()) {
    failures << "Track::checkConsistency(): trackPointsWithMeasurement_ has incorrect size" << std::endl;
    // cppcheck-suppress unreadVariable
    consistent = false;
  }

  for (unsigned int i = 0; i < trackPointsWithMeasurement.size(); ++i) {
    if (trackPointsWithMeasurement[i] != trackPointsWithMeasurement_[i]) {
      failures << "Track::checkConsistency(): trackPointsWithMeasurement_ is not correct" << std::endl;
      failures << "has         id " << i << ", address " << trackPointsWithMeasurement_[i] << std::endl;
      failures << "should have id " << i << ", address " << trackPointsWithMeasurement[i] << std::endl;
      // cppcheck-suppress unreadVariable
      consistent = false;
    }
  }

  if (not consistent) {
    throw genfit::Exception(failures.str(), __LINE__, __FILE__);
  }
}


template<unsigned int dimMeas>
void Track<dimMeas>::trackHasChanged() {

  #ifdef DEBUG
  debugOut << "Track::trackHasChanged \n";
  #endif

  if (fitStatuses_.empty())
    return;

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=fitStatuses_.begin(); it!=fitStatuses_.end(); ++it) {
    it->second->setHasTrackChanged();
  }
}


template<unsigned int dimMeas>
void Track<dimMeas>::fillPointsWithMeasurement() {
  trackPointsWithMeasurement_.clear();
  trackPointsWithMeasurement_.reserve(trackPoints_.size());

  for (std::vector<TrackPoint<dimMeas>*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
    if ((*it)->hasRawMeasurements()) {
      trackPointsWithMeasurement_.push_back(*it);
    }
  }
}

template<unsigned int dimMeas>
void Track<dimMeas>::deleteTrackPointsAndFitStatus() {
  for (size_t i = 0; i < trackPoints_.size(); ++i)
    delete trackPoints_[i];

  trackPoints_.clear();
  trackPointsWithMeasurement_.clear();

  for (std::map< const AbsTrackRep*, FitStatus* >::iterator it = fitStatuses_.begin(); it!= fitStatuses_.end(); ++it)
    delete it->second;
  fitStatuses_.clear();
}

} /* End of namespace genfit */
