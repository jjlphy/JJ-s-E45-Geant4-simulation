// -*- C++ -*-

#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <map>

#include <G4VUserDetectorConstruction.hh>
#include <G4RotationMatrix.hh>
#include <G4String.hh>

class G4Element;
class G4Material;
class G4LogicalVolume;
class G4PVPlacement;
class MagneticField;

//_____________________________________________________________________________
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  static G4String ClassName();
  DetectorConstruction();
  ~DetectorConstruction();
  static const std::vector<G4String>& GetSDList() { return s_detector_list; }

private:
  G4int                           m_experiment;
  std::map<G4String, G4Element*>  m_element_map;
  std::map<G4String, G4Material*> m_material_map;
  G4LogicalVolume*                m_world_lv;
  G4double                        m_rotation_angle;
  G4RotationMatrix*               m_rotation_matrix;
  G4bool                          m_check_overlaps;
  MagneticField*                  m_field;
  static std::vector<G4String>    s_detector_list;

private:
  virtual G4VPhysicalVolume* Construct();
  // Materials
  void ConstructElements();
  void ConstructMaterials();
  // Detectors
  void AddNewDetector(G4VSensitiveDetector* sd);
  void ConstructBAC();
  void ConstructBH2();
  void ConstructFTOF();
  void ConstructHTOF();
  void ConstructHypTPC();
  void ConstructKVC();
  void ConstructBVH_U();
  void ConstructBVH_D();
  void ConstructT0(); 
  void ConstructShsMagnet();
  void ConstructTarget();
  void ConstructFieldOutline();
  void ConstructVP();
  // For K1.8
  void ConstructAreaTent();
  void ConstructBC3();
  void ConstructBC4();
  void ConstructK18BeamlineSpectrometer();
  void ConstructKuramaMagnet();
  void ConstructLAC();
  void ConstructSCH();
  void ConstructSDC1();
  void ConstructSDC2();
  void ConstructSDC3();
  void ConstructSDC4();
  void ConstructWC();

  void CheckOverlaps(G4bool flag) { m_check_overlaps = flag; }
};

//_____________________________________________________________________________
inline G4String
DetectorConstruction::ClassName()
{
  static G4String s_name("DetectorConstruction");
  return s_name;
}

#endif
