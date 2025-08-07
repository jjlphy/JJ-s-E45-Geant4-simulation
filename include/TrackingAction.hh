//
//    **********************************
//    *                                *
//    *    TrackingAction.hh    *
//    *                                *
//    **********************************
//

#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"


class TPCAnaManager;

class TrackingAction : public G4UserTrackingAction {
public:  
  TrackingAction(TPCAnaManager* ana);
  virtual ~TrackingAction();
   
    virtual void BeginOfTrackingAction(const G4Track* aTrack);
    virtual void EndOfTrackingAction(const G4Track* aTrack);
  //  virtual void BeginOfTrackingAction(const G4Event* anEvent);
  //  virtual void EndOfTrackingAction(const G4Event* anEvent);
  private:
  TPCAnaManager* AnaManager;
  //    RunAction* Run;  
};

#endif
