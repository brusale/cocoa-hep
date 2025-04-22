#include "CaloEventAction.hh"

#include <algorithm>
#include <numeric>
#include <iostream>
void CaloEventAction::BeginOfCaloEventAction() {
  parentIDs_.clear();
  totDepEne_.clear();
  parentIDEnergyPerTrack_.clear();

  total_deposit_energy = 0.f;
  particles_soa.clear();
  n_particles = -1;
  n_particles_in_event = 0;
  tracksFromBackscattering.clear();
}

void CaloEventAction::AccumulateStepData(const G4Step *step) {
  G4Track *track = step->GetTrack();
  auto pos1 = track->GetVertexPosition();
  auto pos2 = track->GetPosition();

  G4StepPoint *postStepPoint = step->GetPostStepPoint();
  G4StepPoint *preStepPoint = step->GetPreStepPoint();
  //auto calo_start_z = detector_->calo_start_z();  //TODO
  auto calo_start_z = 100.f;
  float limit = 1.f * calo_start_z;

  int particle_index = -1;
  int parent_idx = -1;

  if (particles_soa.trackid_to_idx.find(track->GetParentID()) != particles_soa.trackid_to_idx.end()) {
    parent_idx = particles_soa.trackid_to_idx[track->GetParentID()];
  }

  G4TouchableHandle touch = preStepPoint->GetTouchableHandle();
  G4TouchableHandle touch2 = postStepPoint->GetTouchableHandle();
  G4LogicalVolume *volume = touch->GetVolume()->GetLogicalVolume();
  G4LogicalVolume *volume2 = touch2->GetVolume()->GetLogicalVolume();
  std::string volume_name = touch->GetVolume()->GetName();
  std::string volume_name2 = touch2->GetVolume()->GetName();
  bool isOutside = (volume_name.substr(1,3) != "CAL") and (volume_name2.substr(1,3) != "CAL");
  if (isOutside) return;
  bool isCrossing = (volume_name.substr(1, 3) != "CAL") and (volume_name2.substr(1,3) == "CAL") ;
  bool catch_position = false;
  if (particles_soa.trackid_to_idx.find(track->GetTrackID()) == particles_soa.trackid_to_idx.end()) {
    if (isCrossing) {
      particle_index = (int)particles_soa.particles_vertex_position_x.size();
      tracks_stack[track->GetTrackID()] = particle_index;
      parent_idx = track->GetParentID(); 
      particles_soa.add(
        pos1.x(),
        pos1.y(),
        pos1.z(),
        track->GetMomentumDirection().x(),
        track->GetMomentumDirection().y(),
        track->GetMomentumDirection().z(),
        track->GetKineticEnergy() / 1000.f,
        track->GetParticleDefinition()->GetPDGEncoding(),
        track->GetTrackID(),
        true,
        track->GetTrackID() // we care about the particle that crossed the boundary (i.e. this one)
      );
      catch_position = true;
      n_particles++;
    }
  } else {
    particle_index = particles_soa.trackid_to_idx[track->GetTrackID()];
    catch_position = true;
    if (step->GetTotalEnergyDeposit() != 0) {
      if (particles_soa.particles_first_impact_kinetic_energy[particle_index] == -1) {
        particles_soa.particles_first_impact_kinetic_energy[particle_index] = track->GetKineticEnergy() / 1000.f;
      }
    }
  } 

  if (particle_index == -1) {
    if (tracks_stack.find(track->GetParentID()) != tracks_stack.end()) {
      particle_index = (int)tracks_stack[track->GetParentID()];
      tracks_stack[track->GetTrackID()] = particle_index;
    }
  }

  auto energy = step->GetTotalEnergyDeposit() / 1000.f;
  total_deposit_energy += energy;

  if (particle_index != -1 and energy != 0) {
    if (particles_soa.particles_first_active_impact_position_x[particle_index] == -1) {
      particles_soa.particles_first_active_impact_position_x[particle_index] = pos2.x();
      particles_soa.particles_first_active_impact_position_y[particle_index] = pos2.y();
      particles_soa.particles_first_active_impact_position_z[particle_index] = pos2.z();
      particles_soa.particles_first_active_impact_momentum_direction_x[particle_index] = track->GetMomentumDirection().x();
      particles_soa.particles_first_active_impact_momentum_direction_y[particle_index] = track->GetMomentumDirection().y();
      particles_soa.particles_first_active_impact_momentum_direction_z[particle_index] = track->GetMomentumDirection().z();
    }
    particles_soa.particles_total_energy_deposited_active[particle_index] += energy;

  }

  trackId2Particle.emplace(std::make_pair(track->GetTrackID(), n_particles));
}

void CaloEventAction::EndOfCaloEventAction() {
  unsigned int nTracks = parentIDEnergyPerTrack_.size();
  for (unsigned int i = 0; i < nTracks; ++i) {
    unsigned int currentTrackID = nTracks - i - 1;
    auto parentID = parentIDEnergyPerTrack_[currentTrackID].first;
    float depositedEnergy = 0.f;
    int originalParentID = parentID;
    while (parentID != 0) {
      depositedEnergy += parentIDEnergyPerTrack_[currentTrackID].second;
      parentIDEnergyPerTrack_[currentTrackID].second = 0.f;
      originalParentID = parentID;
      parentID = parentIDEnergyPerTrack_[parentID - 1].first;
    }
    depositedEnergyPerParticle_.push_back(std::make_pair(originalParentID, depositedEnergy));
  }

  for (auto& depositedEnergy : depositedEnergyPerParticle_) {
    int parentID = depositedEnergy.first;
    auto energy = depositedEnergy.second;
    auto it = std::find_if(simClusters_.begin(), simClusters_.end(),
                           [&](SimCluster& simCluster) {
                             return simCluster.parentID() == parentID;
                           });
    if (it == simClusters_.end()) {
      simClusters_.push_back(SimCluster(parentID, energy));
    } else {
      std::cerr << "Error: SimCluster already exists for parent ID " << parentID << std::endl;
    }
  }
  tracksFromBackscattering.clear();
}
