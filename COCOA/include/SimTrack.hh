#ifndef SimTrack_h
#define SimTrack_h

#include <cstdint>
#include "G4ThreeVector.hh"

class SimTrack {
  public:
    SimTrack() = default;
    SimTrack(uint32_t trackID, int parentID, int particleID, int simVertexID, G4ThreeVector momentum) :
      trackID_(trackID), parentID_(parentID), particleID_(particleID), simVertexID_(simVertexID), momentum_(momentum) {}

    uint32_t getTrackID();
    int getParentID();
    int getParticleID();
    int getSimVertexID();
    G4ThreeVector getMomentum();

    void setTrackID(uint32_t trackID) { this->trackID_ = trackID; }
    void setParentID(int parentID) { this->parentID_ = parentID; }
    void setParticleID(int particleID) { this->particleID_ = particleID; }
    void setSimVertexID(int simVertexID) { this->simVertexID_ = simVertexID; }
    void setMomentum(G4ThreeVector momentum) { this->momentum_ = momentum; }

    bool operator==(const SimTrack& other) const { return trackID_ == other.trackID_; }


  private:
    uint32_t trackID_; // from G4
    int parentID_; // ID of the parent track
    int particleID_; // PDG code of the particle
    int simVertexID_; // ID of the vertex created by this track
    G4ThreeVector momentum_; // momentum of the particle

};
#endif
