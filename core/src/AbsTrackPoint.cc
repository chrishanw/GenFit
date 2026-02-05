/**
 * Author: Christian Wessel
 * Creation date: 16.04.2025
 */

#include "Track.h"

#include "AbsTrackPoint.h"
#include "Exception.h"
#include "KalmanFitterInfo.h"
#include "IO.h"

#include <algorithm>

namespace genfit {

  AbsTrackPoint::AbsTrackPoint() :
  dim_(0), sortingParameter_(0), track_(nullptr), thinScatterer_(nullptr)
{
  ;
}

AbsTrackPoint::AbsTrackPoint(const unsigned int dim) :
  dim_(dim), sortingParameter_(0), track_(nullptr), thinScatterer_(nullptr)
{
  ;
}

AbsTrackPoint::AbsTrackPoint(const unsigned int dim, const Track* track) :
  dim_(dim), sortingParameter_(0), track_(track), thinScatterer_(nullptr)
{
  ;
}
AbsTrackPoint::AbsTrackPoint(const unsigned int dim, const Track* track, const double sortingParameter) :
  dim_(dim), sortingParameter_(sortingParameter), track_(track), thinScatterer_(nullptr)
{
  ;
}

AbsTrackPoint::AbsTrackPoint(const AbsTrackPoint& rhs) :
  dim_(rhs.dim_), sortingParameter_(rhs.sortingParameter_), track_(rhs.track_), thinScatterer_(nullptr)
{
  // copy fitterInfos
  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    AbsFitterInfo* fi = it->second->clone();
    fi->setAbsTrackPoint(this);
    setFitterInfo(fi);
  }

  if (rhs.thinScatterer_ != nullptr)
    thinScatterer_.reset(new ThinScatterer(*(rhs.thinScatterer_)));
}

AbsTrackPoint::AbsTrackPoint(const AbsTrackPoint& rhs,
    const std::map<const AbsTrackRep*, AbsTrackRep*>& map,
    const std::vector<const genfit::AbsTrackRep*> * repsToIgnore) :
    dim_(rhs.dim_), sortingParameter_(rhs.sortingParameter_), track_(rhs.track_), thinScatterer_(nullptr)
{

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    if (repsToIgnore != nullptr) {
      if (std::find(repsToIgnore->begin(), repsToIgnore->end(), it->first) != repsToIgnore->end())
        continue;
    }
    AbsFitterInfo* fi = it->second->clone();
    fi->setRep(map.at(it->first));
    fi->setAbsTrackPoint(this);
    setFitterInfo(fi);
  }

  if (rhs.thinScatterer_ != nullptr)
    thinScatterer_.reset(new ThinScatterer(*(rhs.thinScatterer_)));
}


AbsTrackPoint& AbsTrackPoint::operator=(AbsTrackPoint rhs) {
  swap(rhs);

  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    it->second->setAbsTrackPoint(this);
  }

  return *this;
}


void AbsTrackPoint::swap(AbsTrackPoint& other) {
  std::swap(this->dim_, other.dim_);
  std::swap(this->sortingParameter_, other.sortingParameter_);
  std::swap(this->track_, other.track_);
  std::swap(this->rawMeasurements_, other.rawMeasurements_);
  std::swap(this->fitterInfos_, other.fitterInfos_);
  this->thinScatterer_.swap(other.thinScatterer_);
}


AbsTrackPoint::~AbsTrackPoint() {
  // FIXME: We definitely need some smart containers or smart pointers that
  // take care of this, but so far we haven't found a convincing
  // option (2013-07-05).

  std::map< const AbsTrackRep*, AbsFitterInfo* >::iterator it;
  for (it = fitterInfos_.begin(); it != fitterInfos_.end(); ++it)
    delete it->second;
}


std::vector< AbsFitterInfo* > AbsTrackPoint::getFitterInfos() const {
  std::vector< AbsFitterInfo* > retVal;

  if (fitterInfos_.empty())
    return retVal;

  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    retVal.push_back(it->second);
  }

  return retVal;
}


AbsFitterInfo* AbsTrackPoint::getFitterInfo(const AbsTrackRep* rep) const {
  if (!rep)
    rep = track_->getCardinalRep();
  std::map<const AbsTrackRep*, AbsFitterInfo*>::const_iterator it = fitterInfos_.find(rep);
  if (it == fitterInfos_.end())
    return nullptr;
  return fitterInfos_.at(rep);
}


void AbsTrackPoint::setFitterInfo(genfit::AbsFitterInfo* fitterInfo) {
  assert (fitterInfo != nullptr);
  if (hasFitterInfo(fitterInfo->getRep()))
    delete fitterInfos_[fitterInfo->getRep()];

  fitterInfos_[fitterInfo->getRep()] = fitterInfo;
}


void AbsTrackPoint::Print(const Option_t*) const {
  printOut << "genfit::AbsTrackPoint, belonging to Track " << track_ << "; sorting parameter = " << sortingParameter_ << "\n";
  printOut << "contains " << rawMeasurements_.size() << " rawMeasurements and " << getFitterInfos().size() << " fitterInfos for " << fitterInfos_.size() << " TrackReps.\n";

  for (unsigned int i=0; i<rawMeasurements_.size(); ++i) {
    printOut << "RawMeasurement Nr. " << i << "\n";
    rawMeasurements_[i]->Print();
    printOut << "............\n";
  }

  for (std::map< const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    printOut << "FitterInfo for TrackRep " << it->first << "\n";
    it->second->Print();
    printOut << "............\n";
  }

  if (thinScatterer_)
    thinScatterer_->Print();

}


} /* End of namespace genfit */
