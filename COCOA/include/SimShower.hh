#ifndef SimShower_h
#define SimShower_h

#include <vector>
#include "G4ThreeVector.hh" 
#include "Cell_var.hh"

class SimShower {
  public:
    SimShower() {}
    ~SimShower() {}

    unsigned int GetParentID();

    void AddHit(Cell* cell, float fraction, int hit_idx) {
      auto it = std::find(hit_idx_.begin(), hit_idx_.end(), hit_idx);
      if (fraction <= 0 || cell->get_total_energy() <= 0 || it != hit_idx_.end()) return;
      if (hits_.size() == 0) {
        positionAtBoundary_ = G4ThreeVector(cell->get_x(), cell->get_y(), cell->get_z()); 
      }
      hits_.push_back(cell);
      hit_idx_.push_back(hit_idx);
      energy_ += cell->get_total_energy() * fraction;
      hits_fractions_.push_back(fraction); 
    }

    void CalculateBarycenter() {
      float x = 0;
      float y = 0;
      float z = 0;
      for (auto& hit : hits_) {
        float energy = hit->get_total_energy();
        x += hit->get_x() * energy;
        y += hit->get_y() * energy;
        z += hit->get_z() * energy;
      }
      barycenter_ = G4ThreeVector(x, y, z);
      barycenter_ /= energy_;
    }

    float energy() { return energy_; };
    float et() { return et_; };
    std::vector<int> hit_idx() { return hit_idx_; }

    std::vector<std::pair<Cell*, float>> hitsAndFractions() {
      std::vector<std::pair<Cell*, float>> hitsAndFractions;
      for (unsigned int i = 0; i < hits_.size(); ++i) {
        hitsAndFractions.push_back(std::make_pair(hits_[i], hits_fractions_[i]));
      }
      return hitsAndFractions;
    }

    std::vector<std::pair<unsigned int, float>> hitsIdxAndFractions() {
      std::vector<std::pair<unsigned int, float>> hitsIdxAndFractions;
      for (unsigned int i = 0; i < hit_idx_.size(); ++i) {
        hitsIdxAndFractions.push_back(std::make_pair(hit_idx_[i], hits_fractions_[i]));
      }
      return hitsIdxAndFractions;
    }

    std::vector<std::pair<Cell*, float>> hitsAndEnergies() {
      std::vector<std::pair<Cell*, float>> hitsAndEnergies;
      for (unsigned int i = 0; i < hits_.size(); ++i) {
        hitsAndEnergies.push_back(std::make_pair(hits_[i], hits_fractions_[i] * hits_[i]->get_total_energy()));
      }
      return hitsAndEnergies;
    }

    G4ThreeVector barycenter() { return barycenter_; }
    G4ThreeVector positionAtBoundary() { return positionAtBoundary_; }
    
    void SetPDGID(int pdg_id) { pdg_id_ = pdg_id; }
    int GetPDGID() { return pdg_id_; }

    void SetTrackID(unsigned int track_id) { track_id_ = track_id; }
    unsigned int GetTrackID() { return track_id_; }

  private:
    unsigned int parentID_;
    std::vector<Cell*> hits_;
    std::vector<int> hit_idx_;
    std::vector<float> hits_fractions_;
    float energy_ = 0;
    float et_;
    G4ThreeVector barycenter_;
    G4ThreeVector positionAtBoundary_;
    int pdg_id_;
    unsigned int track_id_;

};
#endif