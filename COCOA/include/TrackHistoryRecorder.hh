#ifndef TrackHistoryRecorder_h
#define TrackHistoryRecorder_h

#include "G4Step.hh"
#include "SimTrack.hh"
#include "SimVertex.hh"

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

	void add(int trackID, int parentID, bool isCrossing, bool inCalo, int& n_particles, int trackPdgId, G4Track* track, const G4Step* step) {

    auto parent_it = std::find(ParentID.begin(), ParentID.end(), parentID);
    if ((isCrossing && parent_it == ParentID.end()) || parentID == 0) {
      n_particles++;
      TrackPdgId.push_back(trackPdgId);
      //if (isCrossing && parent_it == ParentID.end()) {
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
      //}
    } else {
      TrackID.push_back(trackID);
		  ParentID.push_back(parentID);
      ParentAtBoundaryID.push_back(n_particles);
    }
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
