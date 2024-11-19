#ifndef TrackHistoryRecorder_h
#define TrackHistoryRecorder_h

#include "G4Step.hh"
#include "SimTrack.hh"
#include "SimVertex.hh"

class TrackHistoryRecorder {
  public:
    static TrackHistoryRecorder* GetInstance() {
      if (TrackHistoryRecorder::Instance == nullptr) {
        TrackHistoryRecorder::Instance = new TrackHistoryRecorder();
      }
      return TrackHistoryRecorder::Instance;
    }
    
    ~TrackHistoryRecorder() {}

    void RecordStep(std::vector<SimVertex>&, const G4Step*);
 
   private:
    TrackHistoryRecorder() {}
    static TrackHistoryRecorder* Instance;

    unsigned int trackIdx = 0;

};
#endif
