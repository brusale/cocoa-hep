#ifndef __CALORIMETERCONSTRUCTION_H__
#define __CALORIMETERCONSTRUCTION_H__


#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Types.hh"
#include "G4Material.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "Config_reader_var.hh"

#include <vector>
#include <unordered_map>
#include <utility>

struct CellCoordinates {
    CellCoordinates(int ilayer, int ieta, int iphi) {
        layer = ilayer;
        eta = ieta;
        phi = iphi;
    }

    CellCoordinates() : layer(-1), eta(-1), phi(-1) {}

    int layer;
    int eta;
    int phi;
};

class CalorimeterConstruction {

	public:
        CalorimeterConstruction(G4LogicalVolume* expHallLV, bool fCheck_Overlaps, Geometry_definition Geometry);
        ~CalorimeterConstruction();
        void EndCap_Calorimeter();
        void Barrel_Calorimeter();
        void Build_Barrel_CAL ( int iLayer, int NumberOfPixel, int cellMergeFactor,
				long double d_eta, long double d_phi, 
                                long double r_inn, long double r_out, long double previous_layers_delta_r,
                                G4Material *Material_CAL, G4VisAttributes* VisAtt, 
                                const char *LV,  const char *PL, int direction);
        void Build_EndCap_CAL( int iLayer, int NumberOfPixel, int cellMergeFactor, long double d_phi, 
			       long double r_inn_barrel, long double depth, long double previous_layers_depths,
			       G4Material *Material_CAL, G4VisAttributes* VisAtt, 
			       const char *LV,  const char *PL, int direction );

        std::unordered_map<G4LogicalVolume*, std::vector<std::pair<int, CellCoordinates>>> GetCellsMap() const {
            return cells_map_;
        }
        //std::unordered_map<G4LogicalVolume*, std::vector<std::pair<int, int*>>>& GetCellsMap() {
        //    return cells_map_;
        //}

        //int *CellIndex(const char* cellName , double XPos, double YPos, double ZPos);
        CellCoordinates CellIndex(const char* cellName, double XPos, double YPos, double ZPos);

    private:
        char* Name_creation(char *name, int low_layer, int high_layer);
        template <typename T>
	T           GetMinOrMax(const std::vector<std::vector<T > >& array_2d, bool chooseMin) const;
        int         GetNPixelsMax() const;
        long double GetMinDPhi() const;
        bool        CheckMergeFactor( int NumberOfPixel, int cellMergeFactor ) const;
        std::vector<G4VSolid*>* MergeCells( std::vector<G4VSolid*>* cells_to_merge,
					    int NumberOfPixel,
					    int cellMergeFactor );
        G4LogicalVolume* GlobalLV;
        long double theta_min;
        bool fCheckOverlaps;
        long double length;
        Geometry_definition geometry;
	G4Material *m_iron;

        Config_reader_var& config_json_var = Config_reader_var::GetInstance();
            //*Preporation for visualization
        G4VisAttributes* ECAL1_VisAtt=
        new G4VisAttributes(true, G4Colour(100, 0, 255));
        G4VisAttributes* ECAL2_VisAtt=
        new G4VisAttributes(true, G4Colour(30, 143, 255));
        G4VisAttributes* ECAL3_VisAtt=
        new G4VisAttributes(true, G4Colour(0, 0, 255));
        G4VisAttributes* HCAL1_VisAtt=
        new G4VisAttributes(true, G4Colour(0, 255, 0));
        G4VisAttributes* HCAL2_VisAtt=
        new G4VisAttributes(true, G4Colour(255, 255, 0));
        G4VisAttributes* HCAL3_VisAtt=
        new G4VisAttributes(true, G4Colour(255, 0, 0));
        G4VisAttributes* IronGap_VisAtt=
        new G4VisAttributes(true, G4Colour(255, 255, 255));

        //std::unordered_map<G4LogicalVolume*, std::vector<std::pair<int, int*>>> cells_map_;
        std::unordered_map<G4LogicalVolume*, std::vector<std::pair<int, CellCoordinates>>> cells_map_;
};



#endif // __CALORIMETERCONSTRUCTION_H__
