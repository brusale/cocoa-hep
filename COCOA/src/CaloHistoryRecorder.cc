#include "CaloHistoryRecorder.hh"

#include <algorithm>
#include "G4Track.hh"

CaloHistoryRecorder* CaloHistoryRecorder::Instance = nullptr;

void CaloHistoryRecorder::RecordStep(std::vector<std::pair<int, float>>& parentIDEnergyPerTrack,
                                     const G4Step* aStep) {
  // NB G4Track ID starts from 1
  G4Track* aTrack = aStep->GetTrack();
  int trackID = aTrack->GetTrackID();
  int parentID = aTrack->GetParentID();
  float depEne = aStep->GetTotalEnergyDeposit();

  std::pair<int, float> parentIDEnergy = std::make_pair(parentID, depEne);
  AddToRecord(parentIDEnergyPerTrack, trackID, parentIDEnergy);

  // changing the track ID to the parent ID should collapse the graph
  // this way then the energy computation should be O(n)
  // maybe there is a way of making it O(1)
  aTrack->SetTrackID(parentID);
}

bool CaloHistoryRecorder::crossedCaloBoundary(G4Step* aStep) {
  bool crossed = false;
  auto preStepStatus = aStep->GetPreStepPoint()->GetStepStatus();
  if (preStepStatus == fGeomBoundary) { // need to be calo boundary
    crossed = true;
  }
  return crossed;
}

void CaloHistoryRecorder::AddToRecord(std::vector<std::pair<int, float>>& vec,
                                      unsigned int index, 
                                      std::pair<int, float> value) {
  // indeces may not be sorted
  // check on size before resizing
  if (vec.size() <= index) {
    vec.resize(index, std::make_pair(0, 0));
    vec[index - 1] = value;
  } else {
    vec[index] = value;
  }
}
