#include "../include/HistoryData.hh"

void HistoryData::set_tree_branches(TTree* outTree) {
  outTree->Branch("history_vertex_x", &vertex_x);
  outTree->Branch("history_vertex_y", &vertex_y);
  outTree->Branch("history_vertex_z", &vertex_z);
  outTree->Branch("history_vertex_t", &vertex_t);
  outTree->Branch("history_vertex_parent", &vertex_parent);
  outTree->Branch("history_vertex_idx", &vertex_idx);
  outTree->Branch("history_track_id", &track_id);
  outTree->Branch("history_track_parent", &track_parent);
  outTree->Branch("history_track_pdg_id", &track_pdg_id);
  outTree->Branch("history_track_sim_vertex_id", &track_sim_vertex_id);
  outTree->Branch("history_track_momentum_x", &track_momentum_x);
  outTree->Branch("history_track_momentum_y", &track_momentum_y);
  outTree->Branch("history_track_momentum_z", &track_momentum_z);
  outTree->Branch("history_track_position_x", &track_position_x);
  outTree->Branch("history_track_position_y", &track_position_y);
  outTree->Branch("history_track_position_z", &track_position_z);
}

void HistoryData::fill_history_var() {
  for (auto& vertex : simVertices_) {
    vertex_x.push_back(vertex.position().x());
    vertex_y.push_back(vertex.position().y());
    vertex_z.push_back(vertex.position().z());
    vertex_t.push_back(vertex.time());
    vertex_idx.push_back(vertex.getVertexID());
    vertex_parent.push_back(vertex.getParentID());
  }
  for (auto& track : simTracks_) {
    track_id.push_back(track.getTrackID());
    track_parent.push_back(track.getParentID());
    track_pdg_id.push_back(track.getParticleID());
    track_sim_vertex_id.push_back(track.getSimVertexID());
    track_momentum_x.push_back(track.getMomentum().x());
    track_momentum_y.push_back(track.getMomentum().y());
    track_momentum_z.push_back(track.getMomentum().z());
    track_position_x.push_back(track.getTrackXYZ().x());
    track_position_y.push_back(track.getTrackXYZ().y());
    track_position_z.push_back(track.getTrackXYZ().z());
  }
}
