#include "../include/SimTracksData.hh"
#include <numeric>

void SimTracksData::set_tree_branches(TTree* outTree) {
  outTree->Branch("sim_track_tkx", &track_x);
  outTree->Branch("sim_track_tky", &track_y);
  outTree->Branch("sim_track_tkz", &track_z);
  outTree->Branch("sim_track_momentum_x", &track_px);
  outTree->Branch("sim_track_momentum_y", &track_py);
  outTree->Branch("sim_track_momentum_z", &track_pz);
  outTree->Branch("sim_track_momentum_e", &track_e);
  outTree->Branch("sim_track_genp", &track_gen_particle_idx);
  outTree->Branch("sim_track_charge", &track_charge);
  outTree->Branch("sim_track_type", &track_type);
}

void SimTracksData::fill_track_var() {
  for (auto& track : simTracks_) {
    track_x.push_back(track.getTrackXYZ().x());
    track_y.push_back(track.getTrackXYZ().y());
    track_z.push_back(track.getTrackXYZ().z());
    track_px.push_back(track.getMomentum().x());
    track_py.push_back(track.getMomentum().y());
    track_pz.push_back(track.getMomentum().z());
    track_e.push_back(track.getMomentum().mag());
    track_gen_particle_idx.push_back(track.getParticleID());
    track_charge.push_back(track.getCharge());
  }
  track_type.resize(simTracks_.size(), 0);
  track_vertex_idx.resize(simTracks_.size(), 0);
  std::iota(track_type.begin(), track_type.end(), 0);
  std::iota(track_vertex_idx.begin(), track_vertex_idx.end(), 0); // by construction
}

void SimTracksData::clear() {
  track_x.clear();
  track_y.clear();
  track_z.clear();
  track_px.clear();
  track_py.clear();
  track_pz.clear();
  track_e.clear();
  track_gen_particle_idx.clear();
  track_charge.clear();
  track_type.clear();
  track_vertex_idx.clear();
  simTracks_.clear();
}