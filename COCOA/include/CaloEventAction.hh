#ifndef CaloEventAction_h
#define CaloEventAction_h

#include <vector>
#include <unordered_map>
#include "SimCluster.hh"
//#include "SteppingAction.hh"
#include "CaloHistoryRecorder.hh"

class CaloEventAction {
  public:
    CaloEventAction() {}
    ~CaloEventAction() {}

    void BeginOfCaloEventAction();
    void EndOfCaloEventAction();

    std::vector<int> GetParentIDs() { return parentIDs_; }
    std::vector<float> GetTotDepEne() { return totDepEne_; }

    std::vector<SimCluster> simClusters() { return simClusters_; } 

    std::vector<std::pair<int, float>>& parentIDEnergyPerTrack() { return parentIDEnergyPerTrack_; }

    ParticlesSoA& GetParticlesSoA() { return particles_soa; }
    std::unordered_map<int, int>& GetTracksStack() { return tracks_stack; }
    std::unordered_map<int, int>& GetTrackIDToParticle() { return trackId2Particle; }

    void AccumulateStepData(const G4Step*);

  private:
    std::vector<int> parentIDs_; // size: number of tracks; entry i: parent ID of track i
    std::vector<float> totDepEne_; // size: number of tracks; entry i: total deposited energy of track i
    std::vector<std::pair<int, float>> parentIDEnergyPerTrack_;
    std::vector<std::pair<int, float>> depositedEnergyPerParticle_;
    float total_deposit_energy;
    int n_particles;

    std::vector<SimCluster> simClusters_;

    ParticlesSoA particles_soa;
    std::unordered_map<int, int> tracks_stack;

    std::unordered_map<int, int> trackId2Particle;
};

#endif
