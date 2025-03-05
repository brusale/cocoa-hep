#include "../include/SimVerticesData.hh"
#include <numeric>

void SimVerticesData::set_tree_branches(TTree* outTree) {
  outTree->Branch("sim_vertex_x", &vertex_x);
  outTree->Branch("sim_vertex_y", &vertex_y);
  outTree->Branch("sim_vertex_z", &vertex_z);
  outTree->Branch("sim_vertex_t", &vertex_t);
  outTree->Branch("sim_vertex_parent", &vertex_parent);
  outTree->Branch("sim_vertex_idx", &vertex_idx);
}

void SimVerticesData::fill_vertex_var() {
  for (auto& vertex : simVertices_) {
    vertex_x.push_back(vertex.position().x());
    vertex_y.push_back(vertex.position().y());
    vertex_z.push_back(vertex.position().z());
    vertex_t.push_back(vertex.time());
    vertex_idx.push_back(vertex.getVertexID());
  }
  vertex_parent.resize(simVertices_.size(), 0);
  std::iota(vertex_parent.begin(), vertex_parent.end(), 0); // by construction
}

void SimVerticesData::clear() {
  vertex_x.clear();
  vertex_y.clear();
  vertex_z.clear();
  vertex_t.clear();
  vertex_parent.clear();
  vertex_idx.clear();
  simVertices_.clear();
}