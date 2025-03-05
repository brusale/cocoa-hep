#ifndef __SIMSHOWERSDATA_H__
#define __SIMSHOWERSDATA_H__

#include "TTree.h"
#include "SimShower.hh"

#include <vector>
#include <cstdint>
class SimShowersData {
  public:
    SimShowersData() {}
    SimShowersData(std::vector<SimShower> simShowers) : simShowers_(simShowers) {}
    ~SimShowersData() {}

    static SimShowersData& GetInstance() {
      static SimShowersData instance;
      return instance;
    }
    void set_tree_branches(TTree* outTree);
    void fill_shower_var();
    void clear();

    void getSimShowers(std::vector<SimShower>& simShowers) { simShowers_ = simShowers; }

    private:
      std::vector<float> shower_barycenter_x_;
      std::vector<float> shower_barycenter_y_;
      std::vector<float> shower_barycenter_z_;
      std::vector<float> shower_posBoundary_x_;
      std::vector<float> shower_posBoundary_y_;
      std::vector<float> shower_posBoundary_z_;
      std::vector<float> shower_energy_;
      std::vector<int> shower_pdg_id_;
      std::vector<unsigned int> shower_track_id_;
      std::vector<std::vector<uint32_t>> shower_hits_;
      std::vector<std::vector<float>> shower_fractions_;
      std::vector<int> shower_genp_;
      std::vector<std::vector<unsigned int>> shower_hit_idxs_;

      std::vector<SimShower> simShowers_;
};
#endif