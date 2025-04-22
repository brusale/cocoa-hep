#ifndef __HISTORYDATA_H__
#define __HISTORYDATA_H__

#include <vector>
#include "TTree.h"
#include "SimTrack.hh"
#include "SimVertex.hh"

class HistoryData {
  public:
    HistoryData() {}
    HistoryData(std::vector<SimTrack> simTracks, std::vector<SimVertex> simVertices) : simTracks_(simTracks), simVertices_(simVertices) {}
    ~HistoryData() {}

    static HistoryData& GetInstance() {
      static HistoryData instance;
      return instance;
    }

    void clear() {
      simTracks_.clear();
      simVertices_.clear();
      vertex_x.clear();
      vertex_y.clear();
      vertex_z.clear();
      vertex_t.clear();
      vertex_parent.clear();
      vertex_idx.clear();
      track_id.clear();
      track_parent.clear();
      track_pdg_id.clear();
      track_sim_vertex_id.clear();
      track_momentum_x.clear();
      track_momentum_y.clear();
      track_momentum_z.clear();
      track_position_x.clear();
      track_position_y.clear();
      track_position_z.clear();

    }
    void set_tree_branches(TTree* outTree);
    void fill_history_var();

    void getSimTracks(std::vector<SimTrack>& simTracks) { simTracks_ = simTracks; }
    void getSimVertices(std::vector<SimVertex>& simVertices) { simVertices_ = simVertices; }
  private:
    std::vector<SimTrack> simTracks_;
    std::vector<SimVertex> simVertices_;

    std::vector<float> vertex_x;
    std::vector<float> vertex_y;
    std::vector<float> vertex_z;
    std::vector<float> vertex_t;
    std::vector<int> vertex_parent;
    std::vector<int> vertex_idx;
    std::vector<int> track_id;
    std::vector<int> track_parent;
    std::vector<int> track_pdg_id;
    std::vector<int> track_sim_vertex_id;
    std::vector<float> track_momentum_x;
    std::vector<float> track_momentum_y;
    std::vector<float> track_momentum_z;
    std::vector<float> track_position_x;
    std::vector<float> track_position_y;
    std::vector<float> track_position_z;

};
#endif
