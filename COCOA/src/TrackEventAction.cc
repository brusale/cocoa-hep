#include "TrackEventAction.hh"
#include "Full_trajectory_info_data.hh"
#include "DetectorConstruction.hh"
#include "G4PrimaryParticle.hh"
#include <algorithm>
#include <numeric>
#include <unordered_map>

TrackEventAction* TrackEventAction::Instance = nullptr;

void TrackEventAction::BeginOfTrackEventAction() {
  simVertices_.clear();
  simTracks_.clear();
}

void TrackEventAction::EndOfTrackEventAction() {
  //CreateGraph(simVertices_);
  std::vector<std::vector<int>> connectedVertices = BuildTracks();

  TrackEventAction* trackEventAction = TrackEventAction::GetInstance();
  //TracksStack& tracksStack_ = trackEventAction->tracksStack();
  std::vector<std::pair<int,G4Track>>& tracksStack_ = trackEventAction->tracksStack();
  VerticesStack& verticesStack_ = trackEventAction->verticesStack();
  // Traverse the graph starting from the first node
  // it returns a list of the vertices. the first and last vertices are
  // the start and end point of a track
  //std::vector<std::vector<int>> connectedVertices = TraverseGraph(simVertices_);
  uint32_t nTracks = connectedVertices.size();
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  for (uint32_t i = 0; i < nTracks; i++) {
    std::vector<int>& connectedTracks = connectedVertices[i];
    Full_trajectory_info_data& trajectories = Full_trajectory_info_data::GetInstance();
    for (uint32_t j = 0; j < connectedTracks.size(); j++) {
      G4Track& track = tracksStack_[connectedTracks[j]].second;
      if (track.GetParentID() == 0) {
        FullTrajectoryInfo trjInfo;
        trjInfo.is_conversion_track = false;
        trjInfo.fParentID = track.GetParentID();
        trjInfo.fTrackID = track.GetTrackID();
        trjInfo.fPDGCharge = track.GetDynamicParticle()->GetCharge();
        trjInfo.fPDGCode = track.GetDefinition()->GetPDGEncoding();
        trjInfo.fMomentum = track.GetDynamicParticle()->GetMomentum();
        trjInfo.fMomentumDir = track.GetDynamicParticle()->GetMomentumDirection();
        trjInfo.fEnergy = track.GetDynamicParticle()->GetTotalEnergy();
        trjInfo.fMass = track.GetDynamicParticle()->GetMass();
        trjInfo.fVertexPosition = track.GetVertexPosition();
        trjInfo.fGlobalTime = track.GetGlobalTime();
        trjInfo.caloExtrapolMaxEkin = 0.0;
        trjInfo.caloExtrapolEta = trjInfo.fMomentum.getEta();
        trjInfo.caloExtrapolPhi = GetPhi(trjInfo.fMomentum.x(), trjInfo.fMomentum.y());
        trjInfo.idExtrapolMaxEkin = trjInfo.caloExtrapolMaxEkin;
        trjInfo.idExtrapolEta = trjInfo.caloExtrapolEta;
        trjInfo.idExtrapolPhi = trjInfo.caloExtrapolPhi;
        trjInfo.vTrackMomentumDir.push_back(track.GetMomentum());
        trjInfo.vTrackID.push_back(track.GetTrackID());
        trjInfo.vParentID.push_back(track.GetParentID());
        trjInfo.vTrackPos.push_back(track.GetPosition());
        trjInfo.vTrackTime.push_back(track.GetGlobalTime());
        trjInfo.vTrackPDGID.push_back(track.GetDefinition()->GetPDGEncoding());
        trajectories.fAllTrajectoryInfo.push_back(trjInfo);
      } else {
	      for (std::vector < FullTrajectoryInfo>* _trajectories:{&trajectories.fAllTrajectoryInfo, &trajectories.fAllConvElectrons }){
		      bool foundTraj(false);
		      int mTraj(-1); //, mParent(-1);
		      for (int iTraj = (int)_trajectories->size() - 1; iTraj >= 0; iTraj--) {
			      for (int iParent = (int)_trajectories->at(iTraj).vTrackID.size() - 1; iParent >= 0; iParent--) {
				      if (track.GetParentID() == _trajectories->at(iTraj).vTrackID.at(iParent)) {
					      foundTraj = true;
					      mTraj = iTraj;
					      break;
					    } // if( _trajectories->at(iTraj).vParentID.at(iParent) == ParentID )
				    } // for(int iParent = 0; iParent < (int)_trajectories->at(iTraj).vParentID.size(); iParent++  )
			      if (foundTraj)
				      break;
			    } // for(int iTraj = 0; iTraj < (int)_trajectories->size(); iTraj++)
		      if (foundTraj) {
			      _trajectories->at(mTraj).vTrackMomentumDir.push_back(track.GetMomentum());
			      _trajectories->at(mTraj).vParentID.push_back(track.GetParentID());
			      _trajectories->at(mTraj).vTrackID.push_back(track.GetTrackID());
			      _trajectories->at(mTraj).vTrackPos.push_back(track.GetPosition());
			      _trajectories->at(mTraj).vTrackTime.push_back(track.GetGlobalTime());
			      _trajectories->at(mTraj).vTrackPDGID.push_back(track.GetDefinition()->GetPDGEncoding());
			    } //  if(foundTraj)
		    }
	    }	
    }
  }

  /*for (uint32_t i = 0; i < nTracks; i++) {

    std::vector<int>& connectedVertices_v = connectedVertices[i];
    Full_trajectory_info_data& trajectories = Full_trajectory_info_data::GetInstance();
    for (uint32_t j = 0; j < connectedVertices_v.size(); j++) { // -1 because we do not count the last vertex
      //std::vector<int> daughters = simVertices_[connectedVertices_v[j]].getDaughters();
      //assert(daughters.size() == 1); // there should only be one daughter by construction
      int daughter = daughters[0];
      auto it = std::find_if(tracksStack_.trackID.begin(), tracksStack_.trackID.end(),
                             [&](int trackID) {
                              return trackID == daughter;
                             });
      int idx = std::distance(tracksStack_.trackID.begin(), it);
      FullTrajectoryInfo trjInfo;
      trjInfo.is_conversion_track = false;
      trjInfo.fParentID = tracksStack_.parentID[connectedVertices_v[i]];
      trjInfo.fTrackID = tracksStack_.trackID[idx];
      trjInfo.fPDGCharge = tracksStack_.charge[idx];
      trjInfo.fPDGCode = tracksStack_.pdgCode[idx];
      trjInfo.fMomentum = G4ThreeVector(tracksStack_.momentum_x[idx],
                                        tracksStack_.momentum_y[idx],
                                        tracksStack_.momentum_z[idx]);
      trjInfo.fMomentumDir = G4ThreeVector(tracksStack_.momentumDirection_x[idx],
                                            tracksStack_.momentumDirection_y[idx],
                                            tracksStack_.momentumDirection_z[idx]);
      trjInfo.fEnergy = tracksStack_.energy[idx];
      trjInfo.fMass = tracksStack_.mass[idx];
      trjInfo.fVertexPosition = G4ThreeVector(verticesStack_.position_x[connectedVertices_v[j]],
                                              verticesStack_.position_y[connectedVertices_v[j]],
                                              verticesStack_.position_z[connectedVertices_v[j]]);
      trjInfo.fGlobalTime = tracksStack_.globalTime[idx];

      trjInfo.caloExtrapolMaxEkin = 0.0;
      trjInfo.caloExtrapolEta = trjInfo.fMomentum.getEta();
      trjInfo.caloExtrapolPhi = GetPhi(trjInfo.fMomentum.x(), trjInfo.fMomentum.y());

      trjInfo.idExtrapolMaxEkin = trjInfo.caloExtrapolMaxEkin;
      trjInfo.idExtrapolEta = trjInfo.caloExtrapolEta;
      trjInfo.idExtrapolPhi = trjInfo.caloExtrapolPhi;

      trjInfo.vTrackMomentumDir.push_back(trjInfo.fMomentum);
      trjInfo.vTrackID.push_back(trjInfo.fTrackID);
      trjInfo.vParentID.push_back(trjInfo.fParentID);
      trjInfo.vTrackPos.push_back(G4ThreeVector(tracksStack_.position_x[idx],
                                                tracksStack_.position_y[idx],
                                                tracksStack_.position_z[idx]));
      trjInfo.vTrackTime.push_back(tracksStack_.globalTime[idx]);
      trjInfo.vTrackPDGID.push_back(tracksStack_.pdgCode[idx]);

      trajectories.fAllTrajectoryInfo.push_back(trjInfo);
    }
  }

  // build tracks from the connected vertices
  for (uint32_t i = 0; i < nTracks; i++) {
    SimTrack track;

    // reference track to the start vertex 
    int simVertexID = connectedVertices[i][0];
    SimVertex& startVertex = simVertices_[simVertexID];
    // use end vertex to set the track ID
    SimVertex& endVertex = simVertices_[connectedVertices[i][connectedVertices[i].size() - 1]];
    
    track.setTrackID(endVertex.getParentID());
    track.setSimVertexID(simVertexID);
    track.setParentID(startVertex.getParentID());
    track.setParticleID(endVertex.getParticleID()); 

    simTracks_.push_back(track); 
  }*/
}

std::vector<std::vector<int>> TrackEventAction::BuildTracks() {
  std::vector<std::vector<int>> tracks;
  TrackEventAction* trackEventAction = TrackEventAction::GetInstance();
  std::vector<std::pair<int, G4Track>>& tracksStack = trackEventAction->tracksStack();
  std::unordered_map<int, std::vector<int>>& trackIdToIdx = trackEventAction->getTrackIdToIdx();

  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::cout << "tracksStack.size(): " << tracksStack.size() << std::endl;
  std::cout << "trackIdToIdx.size(): " << trackIdToIdx.size() << std::endl;
  for (auto it = trackIdToIdx.begin(); it != trackIdToIdx.end(); ++it) {
    std::cout << "trackIdToIdx[" << it->first <<"].size(): " << it->second.size() << std::endl;
  }
  for (unsigned int i = 0; i < tracksStack.size(); ++i) {
    std::cout << "track ID: " << tracksStack[i].first << std::endl;
    G4Track& track = tracksStack[i].second;
    std::cout << "track charge: " << track.GetDynamicParticle()->GetCharge() << std::endl;
    std::cout << "track energy: " << track.GetTotalEnergy() << std::endl;
    //std::vector<int> track;
    //if (track.GetDynamicParticle()->GetCharge() == 0) continue;
    int origTrackID = tracksStack[i].first;
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::cout << "origTrackID: " << origTrackID << std::endl;
    std::cout << "trackIdToIdx[origTrackID].size(): " << trackIdToIdx[origTrackID].size() << std::endl;
    tracks.push_back(trackIdToIdx[origTrackID]); 
    for (auto t : trackIdToIdx[origTrackID]) {
      std::cout << t << " ";
    } 
    std::cout << std::endl;
  }
  return tracks;
}

/*std::vector<std::vector<int>> TrackEventAction::BuildTracks() {
  std::vector<std::vector<int>> tracks;
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  TrackEventAction* trackEventAction = TrackEventAction::GetInstance();
  TracksStack& tracks_stack = trackEventAction->tracksStack();
  std::unordered_map<int, std::vector<int>>& track_id_to_idx = trackEventAction->getTrackIdToIdx();
  for (size_t i = 0; i < tracks_stack.trackID.size(); i++) {
    std::vector<int> track;
    if (!tracks_stack.crossedBoundary[i] || tracks_stack.charge[i] == 0) continue;
    track.push_back(tracks_stack.trackID[i]);
    int parentID = tracks_stack.parentID[i];
    while (parentID != 0) {
      int parentIdx = -1;
      auto parentIdxIt = track_id_to_idx.find(parentID);
      if (parentIdxIt == track_id_to_idx.end()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::cout << "Parent ID not found" << std::endl;
      } else {
        for (auto& parent : track_id_to_idx[parentID]) {
          std::cout << parent << " ";
        }
        parentIdx = parentIdxIt->second[0];
      }
      track.push_back(parentID);
      parentID = tracks_stack.parentID[parentIdx];
    }
    tracks.push_back(track); 
  }

  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  for (auto& track : tracks) {
    for (auto t : track) {
      std::cout << t << " ";
    }
    std::cout << std::endl;
  }
  return tracks;
}*/



/*void TrackEventAction::CreateGraph(std::vector<SimVertex>& simVertices) {
  TrackEventAction* trackEventAction = TrackEventAction::GetInstance();
  auto trackIDs = trackEventAction->tracksStack().trackID;
  std::sort(trackIDs.begin(), trackIDs.end());
  auto last = std::unique(trackIDs.begin(), trackIDs.end());
  trackIDs.erase(last, trackIDs.end());

  std::vector<std::vector<int>> connectedVertices;
  for (auto trackID : trackIDs) {
    std::vector<int> connectedVertices_v;
    auto it = std::find_if(simVertices.begin(), simVertices.end(),
                           [&](SimVertex& sv) {
                            return sv.getParentID() == trackID;
                           });
    if (it != simVertices.end()) {
      int v_index = std::distance(simVertices.begin(), it);
      connectedVertices_v.push_back(v_index);
    }
    connectedVertices.push_back(connectedVertices_v);
  }
  for (auto connectedVertices_v : connectedVertices) {
    for (auto v : connectedVertices_v) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }
}

// builds nodes and edges from the simvertices
void TrackEventAction::CreateGraph(std::vector<SimVertex>& simVertices) {
  //std::vector<Node> nodes(simVertices.size());
  nodes_.resize(simVertices.size());
  auto nTracks = std::accumulate(simVertices.begin(), simVertices.end(),
                                 0, [&](int n, SimVertex& sv) {
                                  return n + sv.getDaughters().size();
                                 });
  std::unordered_map<unsigned int, int> trackIndexToTrackID;
  //std::vector<Edge> edges(nTracks);
  //std::vector<Edge> edges;
  edges_.resize(nTracks);
  unsigned int trackIdx = 0;
  for (size_t i = 0; i < simVertices.size(); i++) {
    //node->id_ = i;
    Node node(i);
    auto daughters = simVertices[i].getDaughters();
    nodes_[i] = node;
  }

  for (size_t i = 0; i < simVertices.size(); i++) {
    auto daughters = simVertices[i].getDaughters();
    for (auto daughter : daughters) {
      Edge edge;
      trackIndexToTrackID.emplace(trackIdx, daughter);
      edges_[trackIdx].from = &nodes_[i];
      auto it = std::find_if(simVertices.begin(), simVertices.end(),
                             [&](SimVertex& sv) {
                              return sv.getParentID() == daughter;
                             });
      if (it != simVertices.end()) {
        int v_index = std::distance(simVertices.begin(), it);
        edges_[trackIdx].to = &nodes_[v_index];
      } else {
        edges_[trackIdx].to = nullptr;
      }
      trackIdx++;
      //edges.push_back(edge);
    }
  }
  return;
}

// wrapper function for the graph traversal exposed to the user
// returns a list of the connected vertices
std::vector<std::vector<int>> TrackEventAction::TraverseGraph(std::vector<SimVertex>& simVertices) {
  std::vector<Node> nodes;
  std::vector<Edge> edges;
  std::vector<int> connectedComponents;
  std::vector<std::vector<std::vector<int>>> connectedComponentsList(nodes_.size());

  // Loop over all root vertices (vertices with no parent) and traverse each 
  // corresponding graph. This function returns a list of connected vertices
  for (size_t i = 0; i < simVertices.size(); i++) {
    if (!simVertices[i].noParent()) continue;
    TraverseGraphImpl(nodes_, 
                      edges_, 
                      simVertices, 
                      connectedComponents, 
                      connectedComponentsList, 
                      i, i);
  }
  std::cout << __FILE__ << " " << __LINE__ << std::endl;

  std::vector<std::vector<int>> connectedVertices;
  for (size_t i = 0; i < connectedComponentsList.size(); i++) {
    auto& connectedComponents_v = connectedComponentsList[i];
    if (connectedComponents_v.size() == 0) continue;
    for (auto connectedComponent : connectedComponents_v) {
      std::vector<int> tmp;
      for (auto component : connectedComponent) {
        tmp.push_back(component);
      }
      connectedVertices.push_back(tmp);
    }
  }

  return connectedVertices;
}

void TrackEventAction::TraverseGraphImpl(std::vector<Node>& nodes,
                                             std::vector<Edge>& edges,
                                             std::vector<SimVertex>& simVertices,
                                             std::vector<int>& connectedComponents,
                                             std::vector<std::vector<std::vector<int>>>& connectedComponentsList,
                                             int index, 
                                             size_t root) {
  if (!nodes[index].visited and index != -1) {
    nodes[index].visited = true;
    if(simVertices[index].getDaughters().size() <= 1) {
      auto edge = std::find_if(edges.begin(), edges.end(),
                               [&](Edge& e) {
                                return e.from == &nodes_[index];
                               });
      if (edge == edges.end()) {
        connectedComponentsList[root].push_back(connectedComponents);
        connectedComponents.clear();
        return;
      }
      if (edge->to == nullptr) {
        connectedComponentsList[root].push_back(connectedComponents);
        connecte/dComponents.clear();
        return;
      }
      Node* nextNode = edge->to;
      if (!nextNode->visited) {
        connectedComponents.push_back(nextNode->id_);
        if (simVertices[nextNode->id_].getDaughters().size() > 1) {
          connectedComponentsList[root].push_back(connectedComponents);
          connectedComponents.clear();
        }
        TraverseGraphImpl(nodes, edges, simVertices, connectedComponents, connectedComponentsList, nextNode->id_, root);
      }
    } else {
      connectedComponents.clear();
      root = index;
      for (auto daughter : simVertices[index].getDaughters()) {
        auto node = edges[daughter].to->id_;
        connectedComponents.push_back(node);
        TraverseGraphImpl(nodes, edges, simVertices, connectedComponents, connectedComponentsList, node, root);
      }
    }
  }
}*/