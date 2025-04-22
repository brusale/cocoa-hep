#ifndef __HISTORYTOOL_H__
#define __HISTORYTOOL_H__

#include <vector>
#include <algorithm>

#include "SimTrack.hh"
#include "SimVertex.hh"
#include "G4Track.hh"
#include "G4Step.hh"  

class HistoryTool {
  public:
    HistoryTool() {
      simTracks_.clear();
      simVertices_.clear();
      trackID_.clear();
    };

    static HistoryTool* GetInstance() {
      if (HistoryTool::Instance == nullptr) {
        HistoryTool::Instance = new HistoryTool();
      }
      return HistoryTool::Instance;
    }
 
    void clear() {
      simTracks_.clear();
      simVertices_.clear();
      vertexID_ = 0;
    }

    void AddStep(const G4Track* track) {
      if (std::find(trackID_.begin(), trackID_.end(), track->GetTrackID()) != trackID_.end()) return;
      const G4Step& step = *track->GetStep();
      SimVertex simVertex(step.GetPreStepPoint()->GetPosition(), 
                          step.GetPreStepPoint()->GetGlobalTime(), 
                          track->GetParentID(), 
                          vertexID_);
      simVertices_.push_back(simVertex);
      SimTrack simTrack(track->GetTrackID(), 
                        track->GetParentID(), 
                        track->GetDynamicParticle()->GetPDGcode(), 
                        vertexID_, 
                        track->GetMomentum(), 
                        track->GetPosition(), 
                        track->GetDynamicParticle()->GetCharge());
      vertexID_++;
      simTracks_.push_back(simTrack);
      trackID_.push_back(track->GetTrackID());
    }

    std::vector<SimTrack>& simTracks() { return simTracks_; }
    std::vector<SimVertex>& simVertices() { return simVertices_; }
  private:
    static HistoryTool* Instance;
    std::vector<SimTrack> simTracks_;
    std::vector<SimVertex> simVertices_;    
    std::vector<int> trackID_;
    int vertexID_ = 0;
};
#endif
