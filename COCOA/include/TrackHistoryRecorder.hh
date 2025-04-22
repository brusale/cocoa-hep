#ifndef TrackHistoryRecorder_h
#define TrackHistoryRecorder_h

#include "G4Step.hh"
#include "SimTrack.hh"
#include "SimVertex.hh"
#include <unordered_map>

struct G4TrackInfo {
  std::vector<int> TrackID;
	std::vector<int> ParentID;
	std::vector<int> ParentAtBoundaryID;
  std::vector<int> TrackPdgId;

  // for sim tracks and vertices
  std::vector<float> trackX;
  std::vector<float> trackY;
  std::vector<float> trackZ;
  std::vector<float> trackPx;
  std::vector<float> trackPy;
  std::vector<float> trackPz;
  std::vector<float> trackE;
  std::vector<int> trackGenParticleIdx;
  std::vector<int> trackCharge;

  std::vector<float> vertexX;
  std::vector<float> vertexY;
  std::vector<float> vertexZ;
  std::vector<float> vertexT;
  std::vector<int> vertexGenParticleIdx;
  std::vector<int> vertexIdx;

  std::unordered_map<int, int> trackId2Parent;

  void Add(int trackID, int parentID, bool isCrossing, bool inCalo, int& n_particles, int trackPdgId, G4Track* track, const G4Step* step, bool isBackScattering) {
    auto parent_it = trackId2Parent.find(parentID);
    TrackPdgId.push_back(trackPdgId);
    TrackID.push_back(trackID);
    int parentAtBoundaryID = -1;
    if ((isBackScattering && isCrossing) || isCrossing || parentID == 0) {
      trackId2Parent[trackID] = parentID;
      parentAtBoundaryID = trackID;
      trackX.push_back(step->GetPostStepPoint()->GetPosition().x());
      trackY.push_back(step->GetPostStepPoint()->GetPosition().y());
      trackZ.push_back(step->GetPostStepPoint()->GetPosition().z());
      trackPx.push_back(track->GetMomentum().x());
      trackPy.push_back(track->GetMomentum().y());
      trackPz.push_back(track->GetMomentum().z());
      trackE.push_back(track->GetTotalEnergy());
      trackGenParticleIdx.push_back(n_particles);
      trackCharge.push_back(track->GetDynamicParticle()->GetCharge());
      vertexX.push_back(step->GetPostStepPoint()->GetPosition().x());
      vertexY.push_back(step->GetPostStepPoint()->GetPosition().y());
      vertexZ.push_back(step->GetPostStepPoint()->GetPosition().z());
      vertexT.push_back(step->GetPostStepPoint()->GetGlobalTime());
      vertexGenParticleIdx.push_back(n_particles);
      vertexIdx.push_back(n_particles);
    } else if (isBackScattering && !isCrossing) {
      parentAtBoundaryID = trackID;
    } else {
      if (parent_it != trackId2Parent.end()) {
        ParentAtBoundaryID.push_back(parent_it->second);
        trackId2Parent[trackID] = parent_it->second;
      }
    }
    ParentAtBoundaryID.push_back(parentAtBoundaryID);
  }

	void add(int trackID, int parentID, bool isCrossing, bool inCalo, int& n_particles, int trackPdgId, G4Track* track, const G4Step* step, bool isBackScattering) {

    TrackID.push_back(trackID);
    if (isCrossing) {
        ParentAtBoundaryID.push_back(parentID);
        TrackPdgId.push_back(trackPdgId);
    } else if (inCalo) {
      auto it = std::find(ParentAtBoundaryID.begin(), ParentAtBoundaryID.end(), parentID);
      if (it == ParentAtBoundaryID.end()) return; 
      int index = std::distance(ParentAtBoundaryID.begin(), it);
      ParentAtBoundaryID.push_back(ParentAtBoundaryID[index]);
      TrackPdgId.push_back(TrackPdgId[index]);
    } else {
      ParentAtBoundaryID.push_back(-1);
      TrackPdgId.push_back(trackPdgId);
    }

    if (isCrossing) {
      trackX.push_back(step->GetPostStepPoint()->GetPosition().x());
      trackY.push_back(step->GetPostStepPoint()->GetPosition().y());
      trackZ.push_back(step->GetPostStepPoint()->GetPosition().z());
      trackPx.push_back(track->GetMomentum().x());
      trackPy.push_back(track->GetMomentum().y());
      trackPz.push_back(track->GetMomentum().z());
      trackE.push_back(track->GetTotalEnergy());
      trackGenParticleIdx.push_back(trackID);
      trackCharge.push_back(track->GetDynamicParticle()->GetCharge());
      vertexX.push_back(step->GetPostStepPoint()->GetPosition().x());
      vertexY.push_back(step->GetPostStepPoint()->GetPosition().y());
      vertexZ.push_back(step->GetPostStepPoint()->GetPosition().z());
      vertexT.push_back(step->GetPostStepPoint()->GetGlobalTime());
      vertexGenParticleIdx.push_back(trackID);
      vertexIdx.push_back(trackID);
    }
	}

  void clear() {
    TrackID.clear();
    ParentID.clear();
    ParentAtBoundaryID.clear();
    trackX.clear();
    trackY.clear();
    trackZ.clear();
    trackPx.clear();
    trackPy.clear();
    trackPz.clear();
    trackE.clear();
    trackGenParticleIdx.clear();
    vertexX.clear();
    vertexY.clear();
    vertexZ.clear();
    vertexT.clear();
    vertexGenParticleIdx.clear();
    vertexIdx.clear();
    trackCharge.clear();
    TrackPdgId.clear();
    trackId2Parent.clear();
  }
};

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

    G4TrackInfo& GetTrackInfo() { return trackInfo; }

   private:
    TrackHistoryRecorder() {}
    static TrackHistoryRecorder* Instance;

    unsigned int trackIdx = 0;

    G4TrackInfo trackInfo;

};
#endif
