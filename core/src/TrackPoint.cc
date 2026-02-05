/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#include "AbsMeasurement.h"
#include "TrackPoint.h"
#include "Exception.h"
#include "IO.h"

#include <algorithm>

namespace genfit {

template<unsigned int dimMeas>
TrackPoint<dimMeas>::TrackPoint() : AbsTrackPoint(dimMeas)
{
  ;
}

template<unsigned int dimMeas>
TrackPoint<dimMeas>::TrackPoint(Track* track) :
  AbsTrackPoint(dimMeas, track)
{
  ;
}

template<unsigned int dimMeas>
TrackPoint<dimMeas>::TrackPoint(const std::vector< genfit::AbsMeasurement<dimMeas>* >& rawMeasurements, Track* track) :
  AbsTrackPoint(dimMeas, track)
{
  rawMeasurements_.reserve(rawMeasurements.size());

  for (typename std::vector<AbsMeasurement<dimMeas>*>::const_iterator m = rawMeasurements.begin(); m != rawMeasurements.end(); ++m) {
    addRawMeasurement(*m);
  }
}

template<unsigned int dimMeas>
TrackPoint<dimMeas>::TrackPoint(AbsMeasurement<dimMeas>* rawMeasurement, Track* track) :
  AbsTrackPoint(dimMeas, track)
{
  addRawMeasurement(rawMeasurement);
}


template<unsigned int dimMeas>
TrackPoint<dimMeas>::TrackPoint(const TrackPoint<dimMeas>& rhs) : 
  AbsTrackPoint(dimMeas, rhs.track_, rhs.sortingParameter_)
{
  // clone rawMeasurements
  for (typename std::vector<AbsMeasurement<dimMeas>*>::const_iterator it = rhs.rawMeasurements_.begin(); it != rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement<dimMeas>* tp = (*it)->clone();
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (typename std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    AbsFitterInfo* fi = it->second->clone();
    fi->setTrackPoint(this);
    setFitterInfo(fi);
  }

  if (rhs.thinScatterer_ != nullptr)
    Super::thinScatterer_.reset(new ThinScatterer(*(rhs.thinScatterer_)));
}

template<unsigned int dimMeas>
TrackPoint<dimMeas>::TrackPoint(const TrackPoint<dimMeas>& rhs,
    const std::map<const AbsTrackRep*, AbsTrackRep*>& map,
    const std::vector<const genfit::AbsTrackRep*> * repsToIgnore) : AbsTrackPoint(dimMeas, rhs.track_, rhs.sortingParameter_)
{
  // clone rawMeasurements
  for (typename std::vector<AbsMeasurement<dimMeas>*>::const_iterator it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement<dimMeas>* m = (*it)->clone();
    addRawMeasurement(m);
  }

  if (rhs.thinScatterer_ != nullptr)
    Super::thinScatterer_.reset(new ThinScatterer(*(rhs.thinScatterer_)));
}


template<unsigned int dimMeas>
TrackPoint<dimMeas>& TrackPoint<dimMeas>::operator=(TrackPoint<dimMeas> rhs) {
  swap(rhs);
  // TODO: have to swap Super, currently I'm not sure how

  for (typename std::vector<AbsMeasurement<dimMeas>*>::const_iterator it = rawMeasurements_.begin(); it!=rawMeasurements_.end(); ++it) {
    (*it)->setTrackPoint(this);
  }

  for (typename std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    it->second->setTrackPoint(this);
  }

  return *this;
}


template<unsigned int dimMeas>
void TrackPoint<dimMeas>::swap(TrackPoint<dimMeas>& other) {
  std::swap(this->rawMeasurements_, other.rawMeasurements_);
}


template<unsigned int dimMeas>
TrackPoint<dimMeas>::~TrackPoint() {
  // FIXME: We definitely need some smart containers or smart pointers that
  // take care of this, but so far we haven't found a convincing
  // option (2013-07-05).
  
  for (size_t i = 0; i < rawMeasurements_.size(); ++i)
    delete rawMeasurements_[i];

  typename std::map< const AbsTrackRep*, AbsFitterInfo* >::iterator it;
  for (it = Super::fitterInfos_.begin(); it != Super::fitterInfos_.end(); ++it)
    delete it->second;
}


template<unsigned int dimMeas>
AbsMeasurement<dimMeas>* TrackPoint<dimMeas>::getRawMeasurement(int i) const {
  if (i < 0)
    i += rawMeasurements_.size();

  return rawMeasurements_.at(i);
}


// template<unsigned int dimMeas>
// std::vector< AbsFitterInfo* > TrackPoint<dimMeas>::getFitterInfos() const {
//   return Super::getFitterInfos();
// }


// template<unsigned int dimMeas>
// AbsFitterInfo* TrackPoint<dimMeas>::getFitterInfo(const AbsTrackRep* rep) const {
//   // if (!rep)
//   //   rep = track_->getCardinalRep();
//   // typename std::map<const AbsTrackRep*, AbsFitterInfo*>::const_iterator it = fitterInfos_.find(rep);
//   // if (it == fitterInfos_.end())
//   //   return nullptr;
//   // return fitterInfos_.at(rep);
//   return Super::getFitterInfo(rep);
// }



template<unsigned int dimMeas>
void TrackPoint<dimMeas>::deleteRawMeasurements() {
  for (size_t i = 0; i < rawMeasurements_.size(); ++i)
    delete rawMeasurements_[i];

  rawMeasurements_.clear();
}


// template<unsigned int dimMeas>
// void TrackPoint<dimMeas>::setFitterInfo(genfit::AbsFitterInfo* fitterInfo) {
//   assert (fitterInfo != nullptr);
//   // if (hasFitterInfo(fitterInfo->getRep()))
//   //   delete fitterInfos_[fitterInfo->getRep()];

//   // fitterInfos_[fitterInfo->getRep()] = fitterInfo;

//   Super::setFitterInfo(fitterInfo);
// }


template<unsigned int dimMeas>
void TrackPoint<dimMeas>::Print(const Option_t*) const {
  printOut << "genfit::TrackPoint, belonging to Track " << track_ << "; sorting parameter = " << sortingParameter_ << "\n";
  printOut << "contains " << rawMeasurements_.size() << " rawMeasurements and " << getFitterInfos().size() << " fitterInfos for " << fitterInfos_.size() << " TrackReps.\n";

  for (unsigned int i=0; i<rawMeasurements_.size(); ++i) {
    printOut << "RawMeasurement Nr. " << i << "\n";
    rawMeasurements_[i]->Print();
    printOut << "............\n";
  }

  for (typename std::map< const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    printOut << "FitterInfo for TrackRep " << it->first << "\n";
    it->second->Print();
    printOut << "............\n";
  }

  if (thinScatterer_)
    thinScatterer_->Print();

}


} /* End of namespace genfit */
