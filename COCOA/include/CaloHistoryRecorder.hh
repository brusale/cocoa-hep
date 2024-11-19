#ifndef CaloHistoryRecorder_h
#define CaloHistoryRecorder_h

#include <vector>
#include <unordered_map>
#include "G4Step.hh"

class CaloHistoryRecorder {
  public:
    static CaloHistoryRecorder* GetInstance() {
      if (CaloHistoryRecorder::Instance == nullptr) {
        CaloHistoryRecorder::Instance = new CaloHistoryRecorder();
      }
      return CaloHistoryRecorder::Instance;
    }

    ~CaloHistoryRecorder() {}

    void RecordStep(std::vector<std::pair<int, float>>&,
                    const G4Step*);

    bool crossedCaloBoundary(G4Step*);

    void AddToRecord(std::vector<std::pair<int, float>>&, unsigned int, std::pair<int, float>);
  private:
    CaloHistoryRecorder() {}
    static CaloHistoryRecorder* Instance;
};


struct  ParticlesSoA {
  std::vector<float> particles_vertex_position_x;
  std::vector<float> particles_vertex_position_y;
  std::vector<float> particles_vertex_position_z;
  
  std::vector<float> particles_momentum_direction_x;
  std::vector<float> particles_momentum_direction_y;
  std::vector<float> particles_momentum_direction_z;

  std::vector<float> particles_kinetic_energy;
  std::vector<int> particles_pdgid;

  std::vector<float> particles_total_energy_deposited_active;
  std::vector<float> particles_total_energy_deposited_all;

  std::vector<float> particles_first_impact_kinetic_energy;
  std::vector<float> particles_first_active_impact_position_x;
  std::vector<float> particles_first_active_impact_position_y;
  std::vector<float> particles_first_active_impact_position_z;

  std::vector<int> particles_first_active_impact_sensor_idx;

  std::vector<float> particles_first_active_impact_momentum_direction_x;
  std::vector<float> particles_first_active_impact_momentum_direction_y;
  std::vector<float> particles_first_active_impact_momentum_direction_z;

  std::unordered_map<int, int> trackid_to_idx;
  //std::vector<std::pair<int, int>> trackid_to_idx;

  std::vector<bool> particles_tagged;

  std::vector<int> particles_parent_idx;

  void add(float vertex_position_x,
           float vertex_position_y,
           float vertex_position_z,
           float momentum_direction_x,
           float momentum_direction_y,
           float momentum_direction_z,
           float kinetic_energy,
           int pdgid,
           int trackid,
           bool tagged, 
           int parent_idx) {
    
    trackid_to_idx.insert(std::make_pair(trackid, (int)particles_vertex_position_x.size()));
    //trackid_to_idx[trackid] = (int)particles_vertex_position_x.size();
    //trackid_to_idx.push_back(std::make_pair(trackid, (int)particles_vertex_position_x.size()));
    particles_vertex_position_x.push_back(vertex_position_x);
    particles_vertex_position_y.push_back(vertex_position_y);
    particles_vertex_position_z.push_back(vertex_position_z);
    particles_momentum_direction_x.push_back(momentum_direction_x);
    particles_momentum_direction_y.push_back(momentum_direction_y);
    particles_momentum_direction_z.push_back(momentum_direction_z);
    particles_kinetic_energy.push_back(kinetic_energy);
    particles_pdgid.push_back(pdgid);
    particles_total_energy_deposited_active.push_back(0.);
    particles_total_energy_deposited_all.push_back(0.);
    particles_tagged.push_back(tagged);
    particles_parent_idx.push_back(parent_idx);
    particles_first_impact_kinetic_energy.push_back(-1.);
    particles_first_active_impact_position_x.push_back(-1.);
    particles_first_active_impact_position_y.push_back(-1.);
    particles_first_active_impact_position_z.push_back(-1.);
    particles_first_active_impact_sensor_idx.push_back(-1);
    particles_first_active_impact_momentum_direction_x.push_back(-1.);
    particles_first_active_impact_momentum_direction_y.push_back(-1.);
    particles_first_active_impact_momentum_direction_z.push_back(-1.);
  }

  void clear() {
    particles_vertex_position_x.clear();
    particles_vertex_position_y.clear();
    particles_vertex_position_z.clear();
    particles_momentum_direction_x.clear();
    particles_momentum_direction_y.clear();
    particles_momentum_direction_z.clear();
    particles_kinetic_energy.clear();
    particles_pdgid.clear();
    particles_total_energy_deposited_active.clear();
    particles_total_energy_deposited_all.clear();
    particles_first_impact_kinetic_energy.clear();
    particles_first_active_impact_position_x.clear();
    particles_first_active_impact_position_y.clear();
    particles_first_active_impact_position_z.clear();
    particles_first_active_impact_sensor_idx.clear();
    particles_first_active_impact_momentum_direction_x.clear();
    particles_first_active_impact_momentum_direction_y.clear();
    particles_first_active_impact_momentum_direction_z.clear();
    trackid_to_idx.clear();
    particles_tagged.clear();
    particles_parent_idx.clear();
  }   
};

#endif
