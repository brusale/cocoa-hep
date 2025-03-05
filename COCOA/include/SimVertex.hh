#ifndef SimVertex_h
#define SimVertex_h 

#include <vector>
#include "G4ThreeVector.hh"
#include "G4StepPoint.hh"

class SimVertex {
  public:
    SimVertex() = default;
    SimVertex(G4ThreeVector position, float time, int parentID, int vertexId) : position_(position), time_(time), parentID_(parentID), vertexId_(vertexId) {}

    int getParentID() { return parentID_; } 
    void setParentID(int parentID) { this->parentID_ = parentID; }

    int getParticleID() { return particleID_; }
    void setParticleID(int particleID) { this->particleID_ = particleID; }

    bool noParent() { return parentID_ == 0; }

    bool isEndVertex() { return daughters_.size() == 0; }

    G4ThreeVector position() { return position_; }

    float time() { return time_; }

    int getVertexID() { return vertexId_; }

    bool operator==(SimVertex& other) { return ((this->position().x() == other.position().x()) && 
                                                (this->position().y() == other.position().y()) &&
                                                (this->position().z() == other.position().z())) ;}// and (parentID_ == other.parentID_)); }

    void addDaughter(int daughterID) { daughters_.push_back(daughterID); }
    std::vector<int> getDaughters() { return daughters_; }
  private:
    G4ThreeVector position_;
    float time_;
    int parentID_;
    int particleID_; // PDG code of the parent particle
    int vertexId_; // ID of the vertex
    std::vector<int> daughters_; // contains the IDs of the daughter tracks of this vertex
};
#endif
