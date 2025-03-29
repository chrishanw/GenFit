#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace genfit;

// These need no special treatment.
#pragma link C++ class genfit::AbsFinitePlane+;
#pragma link C++ class genfit::AbsHMatrix+;
#pragma link C++ class genfit::RectangularFinitePlane+;
#pragma link C++ class genfit::FitStatus+;
#pragma link C++ class genfit::Material+;
#pragma link C++ class genfit::PruneFlags+;
#pragma link C++ class genfit::TrackCand+;
#pragma link C++ class genfit::TrackCandHit+;
#pragma link C++ class genfit::SharedPlanePtrCreator-;
#pragma link C++ class genfit::AbsTrackRep+;
#pragma link C++ class genfit::MeasuredStateOnPlane+;
#pragma link C++ class genfit::AbsMeasurement+;
#pragma link C++ class genfit::AbsFitterInfo+;
#pragma link C++ class genfit::DetPlane+;
#pragma link C++ class genfit::MeasurementOnPlane+;
#pragma link C++ class genfit::StateOnPlane+;
#pragma link C++ class genfit::ThinScatterer+;
#pragma link C++ class genfit::Track+;
#pragma link C++ class genfit::TrackPoint+;
#pragma link C++ class vector<genfit::TrackPoint*>-;

// Schema Evolution rules.  The official documentation appears to be
// 2010 J. Phys.: Conf. Ser. 219 032004
// http://iopscience.iop.org/1742-6596/219/3/032004
//
// Old versions couldn't actually prune the track, so we ignore the old incarnation
#pragma read sourceClass="genfit::FitStatus" version="[1]" \
  targetClass="genfit::FitStatus"                          \
  source="bool trackIsPruned_;" target="pruneFlags_"       \
  code="{ pruneFlags_.setFlags(); }"
// Prune flag wasn't actually written as no streamer was available.
#pragma read sourceClass="genfit::FitStatus" version="[2]" \
  targetClass="genfit::FitStatus"                          \
  source="" target="pruneFlags_"                           \
  code="{ pruneFlags_.setFlags(); }"

// Time for the TrackCand was only introduced in version 2.  Default to zero.
#pragma read sourceClass="genfit::TrackCand" version="[1]" \
  targetClass="genfit::TrackCand"                          \
  source="" target="time_"                                 \
  code="{ time_ = 0; }"
