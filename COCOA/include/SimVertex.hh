#ifndef SimVertex_h
#define SimVertex_h 

#include <vector>
#include "G4ThreeVector.hh"
#include "G4StepPoint.hh"

class SimVertex {
  public:
    SimVertex() = default;
    SimVertex(G4StepPoint position, int parentID) : position_(position), parentID_(parentID) {}

    int getParentID() { return parentID_; } 
    void setParentID(int parentID) { this->parentID_ = parentID; }

    int getParticleID() { return particleID_; }
    void setParticleID(int particleID) { this->particleID_ = particleID; }

    bool noParent() { return parentID_ == 0; }

    bool isEndVertex() { return daughters_.size() == 0; }

    G4ThreeVector position() { return position_.GetPosition(); }

    bool operator==(SimVertex& other) { return ((this->position().x() == other.position().x()) && 
                                                (this->position().y() == other.position().y()) &&
                                                (this->position().z() == other.position().z())) ;}// and (parentID_ == other.parentID_)); }

    void addDaughter(int daughterID) { daughters_.push_back(daughterID); }
    std::vector<int> getDaughters() { return daughters_; }
  private:
    G4StepPoint position_;
    int parentID_;
    int particleID_; // PDG code of the parent particle
    std::vector<int> daughters_; // contains the IDs of the daughter tracks of this vertex
};
#endif
