/**
 * Author: Christian Wessel
 * Creation date: 25.03.2025
 */

#pragma once

#include <Rtypes.h>
#include "AbsFitterInfo.h"
#include "ThinScatterer.h"

#include <map>
#include <vector>
#include <memory>


namespace genfit {

  class Track;
  class AbsTrackRep;

  /**
   * @brief Abstract TrackPoint class.
   *
   */
  class AbsTrackPoint {

  public:

    AbsTrackPoint();
    AbsTrackPoint(const unsigned int dim = 1);
    explicit AbsTrackPoint(const unsigned int dim = 1, const Track* track = nullptr);
    explicit AbsTrackPoint(const unsigned int dim = 1, const Track* track = nullptr, const double sortingParameter = 0);

    AbsTrackPoint(const AbsTrackPoint&); // copy constructor
    AbsTrackPoint& operator=(AbsTrackPoint); // assignment operator
    void swap(AbsTrackPoint& other);
  
    /**
     * custom copy constructor where all TrackRep pointers are exchanged according to the map.
     * FitterInfos with a rep in repsToIgnore will NOT be copied.
     */
    AbsTrackPoint(const AbsTrackPoint& rhs,
        const std::map<const genfit::AbsTrackRep*, genfit::AbsTrackRep*>& map,
        const std::vector<const genfit::AbsTrackRep*> * repsToIgnore = nullptr);

    virtual ~AbsTrackPoint() {};

    constexpr unsigned int getDim() const { return dim_;}
    // void setDim(unsigned int dim) { dim_ = dim;}

    Track* getTrack() const {return track_;}
    void setTrack(Track* track) {track_ = track;}

    //! Get list of all fitterInfos
    std::vector< genfit::AbsFitterInfo* > getFitterInfos() const;
    //! Get fitterInfo for rep. Per default, use cardinal rep
    AbsFitterInfo* getFitterInfo(const AbsTrackRep* rep = nullptr) const;
    bool hasFitterInfo(const AbsTrackRep* rep) const {
      return (fitterInfos_.find(rep) != fitterInfos_.end());
    }
    //! Takes Ownership
    void setFitterInfo(genfit::AbsFitterInfo* fitterInfo);
    void deleteFitterInfo(const AbsTrackRep* rep) {delete fitterInfos_[rep]; fitterInfos_.erase(rep);}
  
    ThinScatterer* getMaterialInfo() const {return thinScatterer_.get();}
    bool hasThinScatterer() const {return thinScatterer_.get() != nullptr;}
    void setScatterer(ThinScatterer* scatterer) {thinScatterer_.reset(scatterer);}

    double getSortingParameter() const {return sortingParameter_;}
    void setSortingParameter(double sortingParameter) {sortingParameter_ = sortingParameter;}

  private:
    //! Dimension of the measurement of the point
    unsigned int dim_;

    //! Sorting parameter, essentially the number of the point in the track
    double sortingParameter_;

    //! Pointer to Track where TrackPoint belongs to
    Track* track_ = nullptr;

    std::map< const AbsTrackRep*, AbsFitterInfo* > fitterInfos_; //! Ownership over FitterInfos

    std::unique_ptr<ThinScatterer> thinScatterer_; // Ownership

    ClassDef(AbsTrackPoint,1)

  };

} /* End of namespace genfit */
/** @} */
