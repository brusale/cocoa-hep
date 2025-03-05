#ifndef __SIMVERTICESDATA_H__
#define __SIMVERTICESDATA_H__

#include "SimVertex.hh"
#include <vector>
#include "TTree.h"
class SimVerticesData {
  public:
    SimVerticesData() {}
    SimVerticesData(std::vector<SimVertex> simVertices) : simVertices_(simVertices) {}
    ~SimVerticesData() {}

    static SimVerticesData& GetInstance() {
      static SimVerticesData instance;
      return instance;
    }

    void set_tree_branches(TTree* outTree);
    void fill_vertex_var();
    void clear();

    void getSimVertices(std::vector<SimVertex>& simVertices) { simVertices_ = simVertices; }
  private:
    std::vector<float> vertex_x;
    std::vector<float> vertex_y;
    std::vector<float> vertex_z;
    std::vector<float> vertex_t;
    std::vector<int> vertex_idx;
    std::vector<int> vertex_parent;
    std::vector<SimVertex> simVertices_;
};
#endif