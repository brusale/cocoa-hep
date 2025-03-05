#include "SimShowersData.hh"
#include <numeric>
#include "TInterpreter.h"

void SimShowersData::set_tree_branches(TTree* outTree) {
  outTree->Branch("sim_shower_energy", &shower_energy_);
  outTree->Branch("sim_shower_barycenter_x", &shower_barycenter_x_);
  outTree->Branch("sim_shower_barycenter_y", &shower_barycenter_y_);
  outTree->Branch("sim_shower_barycenter_z", &shower_barycenter_z_);
  outTree->Branch("sim_shower_posBoundary_x", &shower_posBoundary_x_);
  outTree->Branch("sim_shower_posBoundary_y", &shower_posBoundary_y_);
  outTree->Branch("sim_shower_posBoundary_z", &shower_posBoundary_z_);
  outTree->Branch("sim_shower_pdg_id", &shower_pdg_id_);
  outTree->Branch("sim_shower_hits", "vector<vector<unsigned int>>", &shower_hits_);
  outTree->Branch("sim_shower_fractions", "vector<vector<float>>", &shower_fractions_);
  outTree->Branch("sim_shower_hit_idxs", "vector<vector<unsigned int>>", &shower_hit_idxs_);
  outTree->Branch("sim_shower_genp", &shower_genp_);
  outTree->Branch("sim_shower_track_id", &shower_track_id_);
}

void SimShowersData::fill_shower_var() {
  for (auto& shower : simShowers_) {
    if (shower.energy() <= 0) continue;
    shower_energy_.push_back(shower.energy());
    shower_barycenter_x_.push_back(shower.barycenter().x());
    shower_barycenter_y_.push_back(shower.barycenter().y());
    shower_barycenter_z_.push_back(shower.barycenter().z());
    shower_posBoundary_x_.push_back(shower.positionAtBoundary().x());
    shower_posBoundary_y_.push_back(shower.positionAtBoundary().y());
    shower_posBoundary_z_.push_back(shower.positionAtBoundary().z());
    shower_pdg_id_.push_back(shower.GetPDGID());
    shower_track_id_.push_back(shower.GetTrackID());

    std::vector<uint32_t> hits, hits_idxs;
    std::vector<float> fractions;
    for (unsigned int i = 0; i < shower.hitsAndFractions().size(); i++) {
      hits.push_back(shower.hitsAndFractions()[i].first->get_cell_id());
      fractions.push_back(shower.hitsAndFractions()[i].second);
      hits_idxs.push_back(shower.hit_idx()[i]); 
    }
    /*for (auto& hit : shower.hitsIdxAndFractions()) {
      hits.push_back(hit.first);
      fractions.push_back(hit.second);
    }*/
    shower_hits_.push_back(hits);
    shower_fractions_.push_back(fractions);
    shower_hit_idxs_.push_back(hits_idxs);
  }
  shower_genp_.resize(simShowers_.size(), 0);
  std::iota(shower_genp_.begin(), shower_genp_.end(), 0); // by construction
}

void SimShowersData::clear() {
  shower_energy_.clear();
  shower_barycenter_x_.clear();
  shower_barycenter_y_.clear();
  shower_barycenter_z_.clear();
  shower_posBoundary_x_.clear();
  shower_posBoundary_y_.clear();
  shower_posBoundary_z_.clear();
  shower_pdg_id_.clear();
  shower_hits_.clear();
  shower_fractions_.clear();
  simShowers_.clear();
  shower_genp_.clear();
  shower_hit_idxs_.clear();
  shower_track_id_.clear();
}
