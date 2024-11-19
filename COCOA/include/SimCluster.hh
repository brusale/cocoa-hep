#ifndef SimCluster_h
#define SimCluster_h

class SimCluster {
  public:
    SimCluster() = default;
    SimCluster(int parentID, float energy) :
      parentID_(parentID), energy_(energy) {}

    int parentID() { return parentID_; }
    void setParentID(int parentID) { parentID_ = parentID; }

    float energy() { return energy_; }
    void setEnergy(float energy) { energy_ = energy; }

  private:
    int parentID_;
    float energy_;
};
#endif