#ifndef __SIMTRACKSDATA_H__
#define __SIMTRACKSDATA_H__

#include <vector>
#include "SimTrack.hh"
#include "TTree.h"

class SimTracksData {
  public:
    SimTracksData() {}
    SimTracksData(std::vector<SimTrack> simTracks) : simTracks_(simTracks) {}
    ~SimTracksData() {}

    static SimTracksData& GetInstance() {
      static SimTracksData instance;
      return instance;
    }

    void set_tree_branches(TTree* outTree);
    void fill_track_var();
    void clear();

    void getSimTracks(std::vector<SimTrack>& simTracks) { simTracks_ = simTracks; }

  private:
    std::vector<float> track_x;
    std::vector<float> track_y;
    std::vector<float> track_z;
    std::vector<float> track_px;
    std::vector<float> track_py;
    std::vector<float> track_pz;
    std::vector<float> track_e;
    std::vector<int> track_gen_particle_idx;
    std::vector<int> track_charge;
    std::vector<int> track_type;
    std::vector<int> track_vertex_idx;

    std::vector<SimTrack> simTracks_;
};

#endif