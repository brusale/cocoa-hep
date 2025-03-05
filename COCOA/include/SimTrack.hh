#ifndef SimTrack_h
#define SimTrack_h

#include <cstdint>
#include "G4ThreeVector.hh"

class SimTrack {
  public:
    SimTrack() = default;
    SimTrack(uint32_t trackID, int parentID, int particleID, int simVertexID, G4ThreeVector momentum, G4ThreeVector trackXYZ, int charge) :
      trackID_(trackID), parentID_(parentID), particleID_(particleID), simVertexID_(simVertexID), momentum_(momentum), trackXYZ_(trackXYZ), charge_(charge) {}

    uint32_t getTrackID() { return trackID_; };
    int getParentID() { return parentID_; };
    int getParticleID() { return particleID_; };
    int getSimVertexID() { return simVertexID_; };
    G4ThreeVector getMomentum() { return momentum_; }
    G4ThreeVector getTrackXYZ() { return trackXYZ_; }
    int getCharge() { return charge_; }

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
    G4ThreeVector trackXYZ_; // position of the particle
    int charge_; // charge of the particle
};
#endif
