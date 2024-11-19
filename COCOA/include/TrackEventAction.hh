#ifndef TrackEventAction_h
#define TrackEventAction_h

#include "SimVertex.hh"
#include "SimTrack.hh" 
#include "TrackHistoryRecorder.hh"
#include "G4Track.hh"
#include "SimVertex.hh"
#include "G4PrimaryParticle.hh"
//#include "SteppingAction.hh"
#include <vector>
#include <unordered_map>
struct Node {
  Node() : id_(-9999999), visited(false) {}
  Node(int id) : id_(id), visited(false) {}
  int id_;
  bool visited;
  bool is_root;
};

struct Edge {
  Edge() : from(nullptr), to(nullptr) {}

  Node* from;
  Node* to;
};

struct TracksStack {
  std::vector<int> trackID;
  std::vector<int> parentID;
  std::vector<float> position_x;
  std::vector<float> position_y;
  std::vector<float> position_z;
  std::vector<float> momentum_x;
  std::vector<float> momentum_y;
  std::vector<float> momentum_z;
  std::vector<float> momentumDirection_x;
  std::vector<float> momentumDirection_y;
  std::vector<float> momentumDirection_z;
  std::vector<float> energy;
  std::vector<float> mass;
  std::vector<float> globalTime;
  std::vector<int> charge;
  std::vector<int> pdgCode;
  std::vector<bool> crossedBoundary;

  void push_back(G4Track* track) {
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    std::cout << "trackID: " << track->GetTrackID() << std::endl;
    std::cout << "parentID: " << track->GetParentID() << std::endl;
    std::cout << "postiion_x: " << track->GetPosition().x() << std::endl;
    std::cout << "postiion_y: " << track->GetPosition().y() << std::endl;
    std::cout << "postiion_z: " << track->GetPosition().z() << std::endl;
    std::cout << "momentum_x: " << track->GetMomentum().x() << std::endl;
    std::cout << "momentum_y: " << track->GetMomentum().y() << std::endl;
    std::cout << "momentum_z: " << track->GetMomentum().z() << std::endl;
    std::cout << "momentum_dir_x: " << track->GetMomentumDirection().x() << std::endl;
    std::cout << "momentum_dir_y: " << track->GetMomentumDirection().y() << std::endl;
    std::cout << "momentum_dir_z: " << track->GetMomentumDirection().z() << std::endl;
    std::cout << "energy: " << track->GetTotalEnergy() << std::endl;
    std::cout << "mass: " << track->GetDynamicParticle()->GetMass() << std::endl;
    std::cout << "charge: " << track->GetDynamicParticle()->GetCharge() << std::endl;
    std::cout << "pdgCode: " << track->GetDefinition()->GetPDGEncoding() << std::endl;
    std::cout << "globalTime: " << track->GetGlobalTime() << std::endl;
    trackID.push_back(track->GetTrackID());
    parentID.push_back(track->GetParentID());
    position_x.push_back(track->GetPosition().x());
    position_y.push_back(track->GetPosition().y());
    position_z.push_back(track->GetPosition().z());
    momentum_x.push_back(track->GetMomentum().x());
    momentum_y.push_back(track->GetMomentum().y());
    momentum_z.push_back(track->GetMomentum().z());
    momentumDirection_x.push_back(track->GetMomentumDirection().x());
    momentumDirection_y.push_back(track->GetMomentumDirection().y());
    momentumDirection_z.push_back(track->GetMomentumDirection().z());
    energy.push_back(track->GetTotalEnergy());
    mass.push_back(track->GetDynamicParticle()->GetMass());
    charge.push_back(track->GetDynamicParticle()->GetCharge());   
    pdgCode.push_back(track->GetDefinition()->GetPDGEncoding()); 
    crossedBoundary.push_back(track->GetNextVolume() != track->GetVolume());
    globalTime.push_back(track->GetGlobalTime());
  }
};

struct VerticesStack {
  std::vector<int> vertexID;
  std::vector<int> parentID;
  std::vector<float> position_x;
  std::vector<float> position_y;
  std::vector<float> position_z;
  std::vector<int> particleID;
  std::vector<std::vector<int>> daughters;

  void push_back(SimVertex& simVertex) {
    vertexID.push_back(parentID.size());
    parentID.push_back(simVertex.getParentID());
    position_x.push_back(simVertex.position().x());
    position_y.push_back(simVertex.position().y());
    position_z.push_back(simVertex.position().z());
    particleID.push_back(simVertex.getParticleID());
    daughters.push_back(simVertex.getDaughters());
  }
};

struct G4PrimaryParticleInfos {
  std::vector<int> primaryParticleTrackID;
  std::vector<double> primaryParticleCharge;
  std::vector<float> primaryParticleMomentum_x;
  std::vector<float> primaryParticleMomentum_y;
  std::vector<float> primaryParticleMomentum_z;
  std::vector<float> primaryParticleMomentumDirection_x;
  std::vector<float> primaryParticleMomentumDirection_y;
  std::vector<float> primaryParticleMomentumDirection_z;
  std::vector<float> primaryParticleTotalEnergy;
  std::vector<float> primaryParticleMass;

  void push_back(G4PrimaryParticle* particle) {
    primaryParticleTrackID.push_back(particle->GetTrackID());
    primaryParticleCharge.push_back(particle->GetCharge());
    primaryParticleMomentum_x.push_back(particle->GetMomentum().x());
    primaryParticleMomentum_y.push_back(particle->GetMomentum().y());
    primaryParticleMomentum_z.push_back(particle->GetMomentum().z());
    primaryParticleMomentumDirection_x.push_back(particle->GetMomentumDirection().x());
    primaryParticleMomentumDirection_y.push_back(particle->GetMomentumDirection().y());
    primaryParticleMomentumDirection_z.push_back(particle->GetMomentumDirection().z());
    primaryParticleTotalEnergy.push_back(particle->GetTotalEnergy());
    primaryParticleMass.push_back(particle->GetMass());
  }
};

class TrackEventAction {
  public:
    TrackEventAction() {}
    ~TrackEventAction() {}

    void BeginOfTrackEventAction();
    void EndOfTrackEventAction();

    std::vector<SimVertex>& simVertices() { return simVertices_; }
    std::vector<SimTrack> simTracks() { return simTracks_; }

    //void CreateGraph(std::vector<SimVertex>& simVertices);
    //std::vector<std::vector<int>> TraverseGraph(std::vector<SimVertex>& simVertices);

    void stashTrack(G4Track* track) { tracksStack_.push_back(std::make_pair(track->GetTrackID(), *track)); }
    void stashVertex(SimVertex& simVertex) { verticesStack_.push_back(simVertex); }

    //TracksStack& tracksStack() { return tracksStack_; }
    std::vector<std::pair<int, G4Track>>& tracksStack() { return tracksStack_; }
    VerticesStack& verticesStack() { return verticesStack_; }
    G4PrimaryParticleInfos& getPrimaryParticleInfos() { return primaryParticleInfos; }

    void fillTrackIdToIdx(int trackID, int idx) {
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      if (trackIdToIdx.find(trackID) == trackIdToIdx.end()) {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        std::vector<int> tmp = {idx};
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        trackIdToIdx.emplace(std::make_pair(trackID, tmp));
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      } else {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        trackIdToIdx[trackID].push_back(idx);
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
      }
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
     }
    std::unordered_map<int, std::vector<int>>& getTrackIdToIdx() { return trackIdToIdx; }
    std::unordered_map<int, int>& getTrackIdToPrimaryParticle() { return trackIdToPrimaryParticle; }

    static TrackEventAction* Instance;
    static TrackEventAction* GetInstance() {
      if (TrackEventAction::Instance == nullptr) {
        TrackEventAction::Instance = new TrackEventAction();
      }
      return TrackEventAction::Instance;
    }

  private:
    std::vector<SimVertex> simVertices_;
    std::vector<SimTrack> simTracks_;

    std::vector<Node> nodes_;
    std::vector<Edge> edges_;

    VerticesStack verticesStack_;
    //TracksStack tracksStack_;
    G4PrimaryParticleInfos primaryParticleInfos;

    std::vector<std::pair<int, G4Track>> tracksStack_;

    std::unordered_map<int, std::vector<int>> trackIdToIdx;
    std::unordered_map<int, int> trackIdToPrimaryParticle;

    void TraverseGraphImpl(std::vector<Node>& nodes, 
                           std::vector<Edge>& edges, 
                           std::vector<SimVertex>& simVertices, 
                           std::vector<int>& connectedComponents, 
                           std::vector<std::vector<std::vector<int>>>& connectedComponentsList, 
                           int index, 
                           size_t root);

    void addEdge(Node* from, Node* to);

    std::vector<std::vector<int>> BuildTracks();

};
#endif
