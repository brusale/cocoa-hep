#include "TrackHistoryRecorder.hh"
#include "TrackEventAction.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include <numeric>

// to be used during SteppingAction
// creates the simvertices to create the graph

TrackHistoryRecorder* TrackHistoryRecorder::Instance = nullptr;

void TrackHistoryRecorder::RecordStep(std::vector<SimVertex>& simVertices,
                                      const G4Step* step) {

  G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
  G4LogicalVolume* lvol = touch->GetVolume()->GetLogicalVolume();

  TrackEventAction* trackEventAction = TrackEventAction::GetInstance();

  std::vector<std::pair<int, G4Track>>& tracksStack = trackEventAction->tracksStack();
  VerticesStack& verticesStack = trackEventAction->verticesStack();


  G4Track* track = step->GetTrack();
  tracksStack.push_back(std::make_pair(track->GetTrackID(),*track));
  //tracksStack.push_back(track);
  std::unordered_map<int, std::vector<int>>& trackid_to_idx = trackEventAction->getTrackIdToIdx();
  //trackid_to_idx.emplace(track->GetTrackID(), tracksStack.size() - 1);
  //trackid_to_idx.emplace(track->GetTrackID(), trackIdx);
  auto track_it = trackid_to_idx.find(track->GetTrackID());
  if (track_it == trackid_to_idx.end()) {
    std::vector<int> tmp{tracksStack.size() - 1};
    trackid_to_idx.emplace(track->GetTrackID(), tmp);
  } else {
    trackid_to_idx[track->GetTrackID()].push_back(tracksStack.size() - 1);
  }

  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::unordered_map<int, int>& trackIdToPrimary = trackEventAction->getTrackIdToPrimaryParticle();
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  int primaryID = track->GetDynamicParticle()->GetPrimaryParticle()->GetTrackID();
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  G4PrimaryParticleInfos& primaryParticleInfos = trackEventAction->getPrimaryParticleInfos();
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  /*if (trackIdToPrimary.find(track->GetTrackID()) == trackIdToPrimary.end()) {
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::vector<int> primaryIDs = primaryParticleInfos.primaryParticleTrackID;
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::cout << "dynamicParticle: " << track->GetDynamicParticle() << std::endl;
    const G4DynamicParticle* dynamicParticle = track->GetDynamicParticle();
    std::cout << "primaryParticle: " << dynamicParticle->GetPrimaryParticle() << std::endl;
    int primaryIDFromTrack = track->GetDynamicParticle()->GetPrimaryParticle()->GetTrackID();
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    auto it = std::find(primaryIDs.begin(), primaryIDs.end(), primaryIDFromTrack);
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    if (it != primaryIDs.end()) {
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      int index = std::distance(primaryIDs.begin(), it);
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      trackIdToPrimary.emplace(track->GetTrackID(), index);
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
    } else {
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      G4PrimaryParticle* particle = track->GetDynamicParticle()->GetPrimaryParticle();
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      primaryParticleInfos.push_back(particle);
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      trackIdToPrimary.emplace(track->GetTrackID(), primaryParticleInfos.primaryParticleTrackID.size() - 1);
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
    }
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
  }*/
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  

  //trackEventAction->fillTrackIdToIdx(track->GetTrackID(), trackIdx);
  //trackIdx++;
  int trackID = track->GetTrackID();
  int trackParentID = track->GetParentID();
  G4StepPoint* trackPreStepPoint = step->GetPreStepPoint();
  SimVertex simVertex(*trackPreStepPoint, trackParentID);
  simVertex.setParticleID(track->GetDefinition()->GetPDGEncoding());

  auto it = std::find_if(simVertices.begin(), simVertices.end(),
                         [&](SimVertex& sv) {
                          return sv == simVertex;
                         });
  if (it == simVertices.end()) {
    simVertex.addDaughter(trackID);
    simVertices.push_back(simVertex);
  } else {
    it->addDaughter(trackID); 
  }

  trackEventAction->stashVertex(simVertex);

  //verticesStack.push_back(simVertex);
  return;
}
