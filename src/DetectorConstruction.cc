// -*- C++ -*-

#include "DetectorConstruction.hh"

#include <G4Box.hh>
#include <G4ChordFinder.hh>
#include <G4Element.hh>
#include <G4FieldManager.hh>
#include <G4IntersectionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Polyhedra.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4SDManager.hh>
#include <G4SubtractionSolid.hh>
#include <G4Transform3D.hh>
#include <G4TransportationManager.hh>
#include <G4Trd.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4UserLimits.hh>
#include <G4VisAttributes.hh>
#include <G4TwoVector.hh>
#include <G4GenericTrap.hh>

#include "BeamMan.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "BACSD.hh"
#include "BH2SD.hh"
#include "FTOFSD.hh"
#include "HTOFSD.hh"
#include "KVCSD.hh"
#include "MagneticField.hh"
#include "TPCSD.hh"
#include "TPCPadSD.hh"
#include "TPCEdepSD.hh"
#include "TargetSD.hh"
#include "VPSD.hh"
#include "padHelper.hh"

#include "BVH_U_SD.hh"
#include "BVH_D_SD.hh"
#include "T0SD.hh"
#include "SCHSD.hh"


namespace
{
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gSize = DetSizeMan::GetInstance();
// color
const G4Colour AQUA(0.247, 0.8, 1.0);
const G4Colour ORANGE(1.0, 0.55, 0.0);
const G4Colour LAVENDER(0.901, 0.901, 0.98);
const G4Colour MAROON(0.5, 0.0, 0.0);
const G4Colour PINK(1.0, 0.753, 0.796);
}

std::vector<G4String> DetectorConstruction::s_detector_list;

//_____________________________________________________________________________
DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    m_experiment(gConf.Get<Int_t>("Experiment")),
    m_element_map(),
    m_material_map(),
    m_world_lv(),
    m_rotation_angle(gConf.Get<Double_t>("SpectrometerAngle")*CLHEP::deg),
    m_rotation_matrix(new G4RotationMatrix),
    m_check_overlaps(false),
    m_field()
{
  m_rotation_matrix->rotateY(- m_rotation_angle);
}

//_____________________________________________________________________________
DetectorConstruction::~DetectorConstruction()
{
}

//_____________________________________________________________________________
G4VPhysicalVolume*
DetectorConstruction::Construct()
{
  using CLHEP::m;

  ConstructElements();
  ConstructMaterials();

  auto world_solid = new G4Box("WorldSolid", 10.*m/2, 6.*m/2, 16.*m/2);
  m_world_lv = new G4LogicalVolume(world_solid, m_material_map["Air"],
                                   "World");
  m_world_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto world_pv = new G4PVPlacement(nullptr, G4ThreeVector(), m_world_lv,
                                    "World", nullptr, false, 0, m_check_overlaps);

  m_field = new MagneticField;
  auto transMan = G4TransportationManager::GetTransportationManager();
  auto fieldMan = transMan->GetFieldManager();
  fieldMan->SetDetectorField(m_field);
  fieldMan->CreateChordFinder(m_field);
  // fieldMan->GetChordFinder()->SetDeltaChord(1.e-3*CLHEP::mm);

  if(gConf.Get<G4int>("Generator") == 10)
    ConstructK18BeamlineSpectrometer();

#if 1
  ConstructBH2();
  ConstructBAC();
  ConstructKVC();
  ConstructBVH_U();
  ConstructBVH_D();
  ConstructT0();
  //ConstructSCH();

#endif

#if 1
  ConstructShsMagnet();
  m_field->Initialize();
#endif

#if 0
  ConstructFieldOutline();
#endif

#if 1
  ConstructTarget();
  ConstructHypTPC();
  ConstructHTOF();
#endif

#if 0
  ConstructVP();
#endif

  G4cout << FUNC_NAME << " SD List Tree" << G4endl;
  G4SDManager::GetSDMpointer()->ListTree();
  return world_pv;
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructElements()
{
  using CLHEP::g;
  using CLHEP::mole;
  /* G4Element(name, symbol, Z, A) */
  G4String name, symbol;
  G4double Z, A;
  name = "Hydrogen";
  m_element_map[name] = new G4Element(name, symbol="H",  Z=1.,
                                      A=1.00794 *g/mole);
  name = "Carbon";
  m_element_map[name] = new G4Element(name, symbol="C",  Z=6.,
                                      A=12.011 *g/mole);
  name = "Nitrogen";
  m_element_map[name] = new G4Element(name, symbol="N",  Z=7.,
                                      A=14.00674 *g/mole);
  name = "Oxygen";
  m_element_map[name] = new G4Element(name, symbol="O",  Z=8.,
                                      A=15.9994 *g/mole);
  name = "Sodium";
  m_element_map[name] = new G4Element(name, symbol="Na", Z=11.,
                                      A=22.989768 *g/mole);
  name = "Silicon";
  m_element_map[name] = new G4Element(name, symbol="Si", Z=14.,
                                      A=28.0855 *g/mole);
  name = "Phoshorus";
  m_element_map[name] = new G4Element(name, symbol="P", Z=15.,
                                      A=30.973762 *g/mole);
  name = "Sulfur";
  m_element_map[name] = new G4Element(name, symbol="S", Z=16.,
                                      A=32.066 *g/mole);
  name = "Argon";
  m_element_map[name] = new G4Element(name, symbol="Ar", Z=18.,
                                      A=39.948 *g/mole);
  name = "Chrominum";
  m_element_map[name] = new G4Element(name, symbol="Cr", Z=24.,
                                      A=51.9961 *g/mole);
  name = "Manganese";
  m_element_map[name] = new G4Element(name, symbol="Mn", Z=25.,
                                      A=54.93805 *g/mole);
  name = "Iron";
  m_element_map[name] = new G4Element(name, symbol="Fe", Z=26.,
                                      A=55.847 *g/mole);
  name = "Nickel";
  m_element_map[name] = new G4Element(name, symbol="Ni", Z=28.,
                                      A=58.69 *g/mole);
  name = "Iodine";
  m_element_map[name] = new G4Element(name, symbol="I",  Z=53.,
                                      A=126.90447 *g/mole);
  name = "Cesium";
  m_element_map[name] = new G4Element(name, symbol="Cs", Z=55.,
                                      A=132.90543 *g/mole);
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructMaterials()
{
  using CLHEP::g;
  using CLHEP::mg;
  using CLHEP::cm3;
  using CLHEP::mole;
  using CLHEP::STP_Temperature;
  G4cout << FUNC_NAME << G4endl;
  /*
    G4Material(name, density, nelement, state, temperature, pressure);
    G4Material(name, z, a, density, state, temperature, pressure);
  */
  G4String name;
  G4double Z, A, density, massfraction;
  G4int natoms, nel, ncomponents;
  const G4double room_temp = STP_Temperature + 20.*CLHEP::kelvin;
  // Vacuum
  name = "Vacuum";
  m_material_map[name] =
    new G4Material(name, density=CLHEP::universe_mean_density, nel=2);
  m_material_map[name]->AddElement(m_element_map["Nitrogen"], 0.7);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], 0.3);
  // Air
  name = "Air";
  m_material_map[name] = new G4Material(name, density=1.2929e-03*g/cm3,
                                        nel=3, kStateGas, room_temp);
  G4double fracN  = 75.47;
  G4double fracO  = 23.20;
  G4double fracAr =  1.28;
  G4double denominator = fracN + fracO + fracAr;
  m_material_map[name]->AddElement(m_element_map["Nitrogen"],
                                   massfraction=fracN/denominator);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],
                                   massfraction=fracO/denominator);
  m_material_map[name]->AddElement(m_element_map["Argon"],
                                   massfraction=fracAr/denominator);
  // Water
  name = "Water";
  m_material_map[name] = new G4Material(name, density=1.*g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], 2);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], 1);
  // Aluminum
  name = "Aluminum";
  m_material_map[name] = new G4Material(name, Z=13., A=26.9815*g/mole,
                                        density=2.699*g/cm3);
  // Iron
  name = "Iron";
  m_material_map[name] = new G4Material(name, Z=26., A=55.85*g/mole,
                                        density=7.87*g/cm3);
  // Copper
  name = "Copper";
  m_material_map[name] = new G4Material(name, Z=29., A=63.546*g/mole,
                                        density=8.96*g/cm3);
  // Carbon
  name = "Carbon";
  m_material_map[name] = new G4Material(name, Z=6., A = 12.0107*g/mole,
                                        density=2.265*g/cm3);
  // Diamond
  name = "Diamond";
  m_material_map[name] = new G4Material(name, Z=6., A = 12.0107*g/mole,
                                        density=3.34*g/cm3);
  // LH2
  name = "LH2";
  m_material_map[name] = new G4Material(name, Z=1., A=1.008*g/mole,
                                        density=70.99*mg/cm3);
  // LD2
  name ="LD2";
  m_material_map[name] = new G4Material(name, Z=1., A=2.01410*g/mole,
                                        density=166.0*mg/cm3);
  // Ar gas
  name = "Argon";
  G4double densityAr = 1.782e-03 * g/cm3 * STP_Temperature / room_temp;
  density = densityAr;
  m_material_map[name] = new G4Material(name, Z=18., A=39.948*g/mole,
                                        density, kStateGas, room_temp);
  // Ethane (C2H6)
  name = "Ethane";
  G4double densityEthane = 1.356e-3 *g/cm3 * STP_Temperature / room_temp;
  density = densityEthane;
  m_material_map[name] = new G4Material(name, density, nel=2,
                                        kStateGas, room_temp);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=2);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms=6);
  // Methane (CH4)
  name = "Methane";
  G4double densityMethane = 0.717e-3 *g/cm3 * STP_Temperature / room_temp;
  density = densityMethane;
  m_material_map[name] = new G4Material(name, density, nel=2,
                                        kStateGas, room_temp);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms=4);
  // Ar(50%) + Ethane(50%) mixture
  name = "ArEthane";
  density = 0.5*(densityAr + densityEthane);
  m_material_map[name] = new G4Material(name, density, ncomponents=2,
                                        kStateGas, room_temp);
  m_material_map[name]->AddMaterial(m_material_map["Argon"],
                                    massfraction=0.5*densityAr/density);
  m_material_map[name]->AddMaterial(m_material_map["Ethane"],
                                    massfraction=0.5*densityEthane/density);
  // P10 gas Ar(90%) + Methane(10%) mixture
  name = "P10";
  density = 0.9*densityAr + 0.1*densityMethane;
  m_material_map[name] = new G4Material(name, density, nel=2,
                                        kStateGas, room_temp);
  m_material_map[name]->AddMaterial(m_material_map["Argon"],
                                    massfraction=0.9*densityAr/density);
  m_material_map[name]->AddMaterial(m_material_map["Methane"],
                                    massfraction=0.1*densityMethane/density);
  // G10 epoxy glass
  name = "G10";
  m_material_map[name] = new G4Material(name, density=1.700*g/cm3,
                                        ncomponents=4);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"] , natoms=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"] , natoms=3);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"] , natoms=3);
  // Kapton
  name = "Kapton";
  m_material_map[name] = new G4Material(name, density=1.42*g/cm3,
                                        ncomponents=4);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"],
                                   massfraction=0.0273);
  m_material_map[name]->AddElement(m_element_map["Carbon"],
                                   massfraction=0.7213);
  m_material_map[name]->AddElement(m_element_map["Nitrogen"],
                                   massfraction=0.0765);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],
                                   massfraction=0.1749);

  // Scintillator (Polystyene(C6H5CH=CH2))
  name = "Scintillator";
  m_material_map[name] = new G4Material(name, density=1.032*g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=8);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms=8);
  // EJ-232 (Plastic Scintillator, Polyvinyltoluene)
  name = "EJ232";
  m_material_map[name] = new G4Material(name, density=1.023*g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=9);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms=12);
  // CH2 Polyethelene
  name = "CH2";
  m_material_map[name] = new G4Material(name, density=0.95*g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms=2);
  // Silica Aerogel for LAC
  name = "SilicaAerogelLAC";
  m_material_map[name] = new G4Material(name, density=0.18 *g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],  natoms=2);
  // Silica Aerogel for BAC
  name = "AerogelBAC";
  m_material_map[name] = new G4Material(name, density=0.377 *g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],  natoms=2);
  // Quartz for KVC (SiO2, crystalline)
  name = "QuartzKVC";
  m_material_map[name] = new G4Material(name, density=2.64 *g/cm3, nel=2);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms=1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],  natoms=2);
  // Acrylic for WC
  name = "Acrylic";
  m_material_map[name] = new G4Material(name, density=1.18 *g/cm3, nel=3);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms=5);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"],  natoms=8);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],  natoms=2);
  

  G4String target_material = gConf.Get<G4String>("TargetMaterial");
  G4cout << "   Target material : " << target_material << G4endl;
  auto itr_target = m_material_map.find(target_material);
  if (itr_target == m_material_map.end()) {
    G4String e(FUNC_NAME + " No such material : " + target_material);
    throw std::invalid_argument(e);
  } else {
    m_material_map["Target"] = itr_target->second;
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::AddNewDetector(G4VSensitiveDetector* sd)
{
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  s_detector_list.push_back(sd->GetName());
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructBAC()
{
  using CLHEP::mm;
  using CLHEP::deg;
  const auto& ra2 = gGeom.GetRotAngle2("BAC")*deg;
  const auto& half_size = gSize.GetSize("BacRadiator")*mm/2.;
  auto pos = gGeom.GetGlobalPosition("BAC");
  auto bacSD = new BACSD("BAC");
  bacSD->SetRefractiveIndex(1.10);
  AddNewDetector(bacSD);
  auto mother_solid = new G4Box("BacMotherSolid",
                                half_size.x() + 10*mm,
                                half_size.y() + 10*mm,
                                half_size.z() + 10*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "BacMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(- ra2 - m_rotation_angle);
  pos.rotateY(m_rotation_angle);
  new G4PVPlacement(rot, pos, mother_lv,
                    "BacMotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Radiator
  auto bac_solid = new G4Box("BacSolid",
                             half_size.x(), half_size.y(), half_size.z());
  auto bac_lv = new G4LogicalVolume(bac_solid, m_material_map["AerogelBAC"], "BacLV");
  new G4PVPlacement(nullptr, G4ThreeVector(), bac_lv, "BacPV",
                    mother_lv, false, 0, m_check_overlaps);
  bac_lv->SetVisAttributes(G4Colour::Yellow());
  bac_lv->SetSensitiveDetector(bacSD);
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructBH2()
{
  using CLHEP::mm;
  using CLHEP::deg;
  auto pos = gGeom.GetGlobalPosition("BH2");
  const auto& ra2 = gGeom.GetRotAngle2("BH2")*deg;
  G4double effective_x = 0.;
  for(G4int i=0; i<NumOfSegBH2; ++i){
    effective_x += gSize.Get("Bh2SegWidth", i)*mm;
  }
  const auto& half_size = gSize.GetSize("Bh2Seg")*mm/2.;
  auto bh2SD = new BH2SD("BH2");
  AddNewDetector(bh2SD);
  // Mother
  auto mother_solid = new G4Box("Bh2MotherSolid",
                                effective_x/2 + 10.*mm,
                                half_size.y() + 10.*mm,
                                half_size.z() + 10.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "Bh2MotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(- ra2 - m_rotation_angle);
  pos.rotateY(m_rotation_angle);
  new G4PVPlacement(rot, pos, mother_lv,
                    "Bh2MotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto inner_segment_solid = new G4Box("Bh2InnerSegmentSolid", half_size.x(),
                                       half_size.y(), half_size.z());
  auto inner_segment_lv = new G4LogicalVolume(inner_segment_solid,
                                              m_material_map["EJ232"],
                                              "Bh2InnerSegmentLV");
  auto outer_segment_solid = new G4Box("Bh2OuterSegmentSolid",
                                       gSize.Get("Bh2SegWidth")*mm/2,
                                       half_size.y(), half_size.z());
  auto outer_segment_lv = new G4LogicalVolume(outer_segment_solid,
                                              m_material_map["EJ232"],
                                              "Bh2OuterSegmentLV");
  pos.set(-effective_x/2 + gSize.Get("Bh2SegWidth")*mm/2, 0., 0.);
  for(G4int i=0; i<NumOfSegBH2; ++i){
    if(i>0){
      pos.setX(pos.x() +
               gSize.Get("Bh2SegWidth", i-1)*mm/2. +
               gSize.Get("Bh2SegWidth", i)*mm/2.);
    }
    if(i==0 || i==NumOfSegBH2-1){ // outer segments
      new G4PVPlacement(nullptr, pos, outer_segment_lv,
                        "Bh2SegmentPV", mother_lv, false, i, m_check_overlaps);
    }
    else{ // inner segments
      new G4PVPlacement(nullptr, pos, inner_segment_lv,
                        "Bh2SegmentPV", mother_lv, false, i, m_check_overlaps);
    }
  }
  inner_segment_lv->SetVisAttributes(G4Colour::Cyan());
  inner_segment_lv->SetSensitiveDetector(bh2SD);
  outer_segment_lv->SetVisAttributes(G4Colour::Cyan());
  outer_segment_lv->SetSensitiveDetector(bh2SD);
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructFTOF()
{
  using CLHEP::mm;
  const auto& ra2 = gGeom.GetRotAngle2("TOF") * CLHEP::deg;
  const auto& half_size = gSize.GetSize("FtofSeg") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("TOF") * mm;
  auto ftofSD = new FTOFSD("FTOF");
  AddNewDetector(ftofSD);
  // Mother
  auto mother_solid = new G4Box("FtofMotherSolid",
                                half_size.x()*NumOfSegFTOF + 50.*mm,
                                half_size.y() + 50.*mm,
                                half_size.z()*2 + 50.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "FtofMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(- ra2 - m_rotation_angle);
  auto pos = (gGeom.GetGlobalPosition("KURAMA") +
              gGeom.GetGlobalPosition("TOF"));
  pos.rotateY(m_rotation_angle);
  new G4PVPlacement(rot, pos, mother_lv,
                    "FtofMotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("FtofSegmentSolid", half_size.x(),
                                 half_size.y(), half_size.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        m_material_map["Scintillator"],
                                        "FtofSegmentLV");
  for(G4int i=0; i<NumOfSegFTOF; ++i){
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(ftofSD);
    pos = G4ThreeVector((-NumOfSegFTOF/2 + i)*pitch,
                        0.0,
                        2.*(- i%2 + 0.5)*half_size.z());
    new G4PVPlacement(nullptr, pos, segment_lv,
                      "FtofSegmentPV", mother_lv, false, i, m_check_overlaps);
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructHTOF()
{
  using CLHEP::mm;
  using CLHEP::deg;
  auto htof_sd = new HTOFSD("HTOF");
  AddNewDetector(htof_sd);
  const auto& htof_pos = gGeom.GetGlobalPosition("HTOF");
  const auto& half_size = gSize.GetSize("HtofSeg") * 0.5 * mm;
  const G4double L = gGeom.GetLocalZ("HTOF");
  const G4double dXdW = gGeom.GetWirePitch("HTOF");
  // Segment
  // Sintillator
  auto htof_scintilltor = new G4Box("HtofScint", half_size.x(),
                                    half_size.y(), half_size.z());
  // Light-guides
  const G4ThreeVector lg_size(17.5,50.0,2.0);
  // Upper one
  std::vector<G4TwoVector> upper_lg_vertices;
  upper_lg_vertices.push_back(G4TwoVector(-half_size.x() , -half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(-half_size.x() , half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(half_size.x() , half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(half_size.x() , -half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(-lg_size.x() , -half_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(-lg_size.x() , -half_size.z() + lg_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(lg_size.x() , -half_size.z() + lg_size.z()));
  upper_lg_vertices.push_back(G4TwoVector(lg_size.x() , -half_size.z()));
  auto htof_upper_lg = new G4GenericTrap("HtofLG_upper", lg_size.y(), upper_lg_vertices);
  auto rotM_upper_lg = new G4RotationMatrix;
  rotM_upper_lg->rotateX(90.0 *deg);
  rotM_upper_lg->rotateZ(- 180.0 *deg);
  G4ThreeVector trans_upper_lg(0.*mm , half_size.y() + lg_size.y(),  0.*mm);

  // Lower one
  std::vector<G4TwoVector> lower_lg_vertices;
  lower_lg_vertices.push_back(G4TwoVector(-half_size.x() , -half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(-half_size.x() , half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(half_size.x() , half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(half_size.x() , -half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(-lg_size.x() , half_size.z() - lg_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(-lg_size.x() , half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(lg_size.x() , half_size.z()));
  lower_lg_vertices.push_back(G4TwoVector(lg_size.x() , half_size.z() - lg_size.z()));
  auto htof_lower_lg = new G4GenericTrap("HtofLG_lower", lg_size.y(), lower_lg_vertices);
  auto rotM_lower_lg = new G4RotationMatrix;
  rotM_lower_lg->rotateX(- 90.0 *deg);
  rotM_lower_lg->rotateZ(- 180.0 *deg);
  G4ThreeVector trans_lower_lg(0.*mm , -half_size.y() - lg_size.y() , 0.*mm);
  auto solid_dummy = new G4UnionSolid("dummysolid", htof_scintilltor, htof_upper_lg,
                                      rotM_upper_lg, trans_upper_lg);

  auto htof_solid = new G4UnionSolid("HtofSolid", solid_dummy,
                                     htof_lower_lg, rotM_lower_lg,
                                     trans_lower_lg);
  //Common slats
  auto htof_lv = new G4LogicalVolume(htof_solid, m_material_map["Scintillator"],
                                     "HtofLV");

  //HTOF beam-through part
  G4double HTOF_window=112.*mm;
  G4double HTOF_upper_lower_y_diff = 12.0*mm;
  
  auto window_dummy_upper = new G4Box("window_dummy", half_size.x(), half_size.y()/2. - HTOF_window/4. + HTOF_upper_lower_y_diff/4., half_size.z());
  auto window_dummy_lower = new G4Box("window_dummy", half_size.x(), half_size.y()/2. - HTOF_window/4. - HTOF_upper_lower_y_diff/4., half_size.z());
  

  //Upper slats(Beam-through)
  G4ThreeVector trans_upper_window(0.*mm , half_size.y()/2. - HTOF_window/4. + lg_size.y() + HTOF_upper_lower_y_diff/4.,  0.*mm);
  auto htof_solid_upper = new G4UnionSolid("HtofSolid_upper", window_dummy_upper, htof_upper_lg,
                                           rotM_upper_lg, trans_upper_window);
  auto htof_upper_lv = new G4LogicalVolume(htof_solid_upper, m_material_map["Scintillator"],
                                           "HtofUpperLV");
  //Lower slats(Beam-through)
  G4ThreeVector trans_lower_window(0.*mm , -half_size.y()/2. + HTOF_window/4. - lg_size.y() + HTOF_upper_lower_y_diff/4.,  0.*mm);
  auto htof_solid_lower = new G4UnionSolid("HtofSolid_lower", window_dummy_lower, htof_lower_lg,
                                           rotM_lower_lg, trans_lower_window);
  auto htof_lower_lv = new G4LogicalVolume(htof_solid_lower, m_material_map["Scintillator"],
                                           "HtofLowerLV");

  for(G4int i=0; i<NumOfPlaneHTOF; ++i){
    for(G4int j=0; j<NumOfSegHTOFOnePlane; ++j){
      G4int seg = i*NumOfSegHTOFOnePlane + j;
      // (lateral, height, radial)
      G4ThreeVector seg_pos(-dXdW * (j - (NumOfSegHTOFOnePlane - 1)/2.),
                            0.*mm,
                            -L);
      auto rotMOutP = new G4RotationMatrix;
      rotMOutP->rotateY(- i * 360./NumOfPlaneHTOF*deg);
      seg_pos.rotateY(i * 360./NumOfPlaneHTOF*deg);
      seg_pos += htof_pos;

      G4int copy_no = seg+2;

      G4ThreeVector window_pos(0.*mm, half_size.y()/2. + HTOF_window/4., 0.*mm);
      //common slats
      if(i!=0)	new G4PVPlacement(rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", copy_no), m_world_lv, false, copy_no, m_check_overlaps);
      else if(j==0) new G4PVPlacement(rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", 0), m_world_lv, false, 0, m_check_overlaps);
      else if(j==3) new G4PVPlacement(rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", 5), m_world_lv, false, 5, m_check_overlaps);
      //Beam-through slats
      else if(j==1){
	new G4PVPlacement(rotMOutP, seg_pos + window_pos, htof_upper_lv, Form("HtofPV%d", 1), m_world_lv, false, 1, m_check_overlaps);
	new G4PVPlacement(rotMOutP, seg_pos - window_pos, htof_lower_lv, Form("HtofPV%d", 2), m_world_lv, false, 2, m_check_overlaps);
      }
      else if(j==2){
	new G4PVPlacement(rotMOutP, seg_pos + window_pos, htof_upper_lv, Form("HtofPV%d", seg), m_world_lv, false, 3, m_check_overlaps);
	new G4PVPlacement(rotMOutP, seg_pos - window_pos, htof_lower_lv, Form("HtofPV%d", seg+31), m_world_lv, false, 4, m_check_overlaps);
      }
    }
  }

  htof_lv->SetSensitiveDetector(htof_sd);
  htof_upper_lv->SetSensitiveDetector(htof_sd);
  htof_lower_lv->SetSensitiveDetector(htof_sd);
  htof_lv->SetVisAttributes(G4Colour::Cyan());
  htof_upper_lv->SetVisAttributes(G4Colour::Cyan());
  htof_lower_lv->SetVisAttributes(G4Colour::Cyan());

  // Supporting frame parts
  // Dummy for subtraction
  auto RingHoleOut = new G4Box("RingHoleOut", 127.5*mm, 6.0*mm, 5.*mm);
  auto RingHoleIn = new G4Box("RingHoleIn", 10.*mm, 6.*mm, 5.*mm);
  auto RingHole = new G4SubtractionSolid("RingHole" , RingHoleOut, RingHoleIn);

  // Top Ring
  const G4double RingIn_zPlane[2] = { -5.*mm, 5.*mm };
  const G4double RingIn_rInner[2] = { 0.*mm, 0.*mm };
  const G4double RingIn_rOuter[2] =  { 334.13*mm, 334.13*mm } ;
  auto TopRingOutSolid = new G4Tubs("TopRingOutSolid", 0.*mm, 391.97*mm,
                                    5.*mm, 0.*deg, 360.*deg);
  auto TopRingInSolid = new G4Polyhedra("TopRingInSolid", 22.5*deg, (360. + 22.5)*deg,
                                        8, 2,
                                        RingIn_zPlane, RingIn_rInner, RingIn_rOuter);

  G4SubtractionSolid* TopRingSubSolid[8];
  for(G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ(i * 45.0 * deg);
    G4ThreeVector trans_hole(0*mm, 348.13*mm, 0.*mm);
    trans_hole.rotateZ(- i * 45.0 * deg);
    if(i==0) TopRingSubSolid[0] = new G4SubtractionSolid("TopRingSubSolid_0", TopRingOutSolid, RingHole, rotMHole, trans_hole);
    else TopRingSubSolid[i] = new G4SubtractionSolid(Form("TopRingSubSolid_%d", i), TopRingSubSolid[i-1], RingHole, rotMHole, trans_hole);
  }

  auto TopRingSolid = new G4SubtractionSolid("TopRingSoild" , TopRingSubSolid[7], TopRingInSolid);
  auto TopRing_lv = new G4LogicalVolume(TopRingSolid, m_material_map["Iron"], "TopRingLV");
  TopRing_lv->SetVisAttributes(ORANGE);

  // Bottom Ring
  const G4double RingOut_zPlane[2] = { -5.*mm, 5.*mm };
  const G4double RingOut_rInner[2] = { 334.13*mm, 334.13*mm };
  const G4double RingOut_rOuter[2] = { 362.13*mm, 362.13*mm };
  auto BotRingOutSolid = new G4Polyhedra("BotRingOutSolid", 22.5*deg, (360. + 22.5)*deg,
                                         8, 2,
                                         RingOut_zPlane, RingOut_rInner, RingOut_rOuter);

  G4SubtractionSolid* BotRingSubSolid[8];
  for(G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ(i * 45.0 * deg);
    G4ThreeVector trans_hole(0.*mm, 348.13*mm, 0.*mm);
    trans_hole.rotateZ(- i * 45.0 * deg);
    if(i==0) BotRingSubSolid[i] = new G4SubtractionSolid(Form("BotRingSubSolid_%d", i) , BotRingOutSolid, RingHole, rotMHole, trans_hole);
    else BotRingSubSolid[i] = new G4SubtractionSolid(Form("BotRingSubSolid_%d", i) , BotRingSubSolid[i-1], RingHole, rotMHole, trans_hole);
  }
  auto BotRing_lv = new G4LogicalVolume(BotRingSubSolid[7], m_material_map["Iron"],"BotRingLV");
  BotRing_lv->SetVisAttributes(G4Colour::Blue());

  // Bracket
  auto BraInSolid = new G4Box("BraInSolid", 137.5*mm, 75.*mm, 2.5*mm);
  auto BraInHole = new G4Box("BraInHole", 127.5*mm, 30.*mm, 2.5*mm);
  G4ThreeVector brain_trans_hole(0.*mm, 25.*mm, 0.*mm);
  auto BraInSubSolid = new G4SubtractionSolid("BraInSubSolid", BraInSolid, BraInHole,
                                              0, brain_trans_hole);

  auto BraOutSolid = new G4Box("BraOutSolid", 137.5*mm, 2.5*mm, 16.5*mm);
  G4ThreeVector braout_trans_hole(0.*mm, 0.*mm, 3.5*mm);
  auto BraOutSubSolid = new G4SubtractionSolid("BraOutSubSolid", BraOutSolid, RingHole,
                                               0, braout_trans_hole);

  G4ThreeVector trans_braket(0.*mm, 72.5*mm, 14.*mm);
  auto BraSolid = new G4UnionSolid("BraSolid", BraInSubSolid, BraOutSubSolid,
                                   0, trans_braket);
  auto Bra_lv = new G4LogicalVolume(BraSolid, m_material_map["Iron"],
                                    "BraLV");
  Bra_lv->SetVisAttributes(G4Colour::Green());

  // Preamp Support Frame
  const G4double PreFrameIn_zPlane[2] = { -37.*mm, 37.*mm };
  const G4double PreFrameIn_rInner[2] = { 334.13*mm, 334.13*mm };
  const G4double PreFrameIn_rOuter[2] = { 339.13*mm, 339.13*mm } ;
  auto PreFrameInSolid = new G4Polyhedra("PreFrameInSolid", 22.5*deg, (360. + 22.5)*deg,
                                         8, 2,
                                         PreFrameIn_zPlane,  PreFrameIn_rInner,  PreFrameIn_rOuter);

  const G4double PreFrameOut_zPlane[2] = { -32.5*mm, 32.5*mm };
  const G4double PreFrameOut_rInner[2] = { 354.13*mm, 354.13*mm };
  const G4double PreFrameOut_rOuter[2] = { 362.13*mm, 362.13*mm } ;
  auto PreFrameOutSolid = new G4Polyhedra("PreFrameOutSolid", 22.5*deg, (360. + 22.5)*deg,
                                          8, 2,
                                          PreFrameOut_zPlane,  PreFrameOut_rInner,  PreFrameOut_rOuter);
  auto PreFrameHole = new G4Box("PreFrameHole", 3.5*mm, 4.0*mm, 32.5*mm);

  G4SubtractionSolid* PreFrameOutSubSolid[8];
  for(G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ(i * 45.0 * deg);
    G4ThreeVector trans_hole(0.*mm, 348.13*mm, 0.*mm);
    trans_hole.rotateZ(- i * 45.0 * deg);
    if(i==0) PreFrameOutSubSolid[0] = new G4SubtractionSolid("PreFrameOutSubSolid_0",
                                                             PreFrameOutSolid, PreFrameHole);
    else PreFrameOutSubSolid[i] = new G4SubtractionSolid(Form("PreFrameOutSubSolid_%d", i),
                                                         PreFrameOutSubSolid[i-1], PreFrameHole,
                                                         rotMHole, trans_hole);
  }
  G4ThreeVector trans_preframe(0.*mm, 0.*mm, -9.0*mm);
  auto PreFrameSolid = new G4UnionSolid("PreFrameSolid", PreFrameInSolid, PreFrameOutSubSolid[7],
                                        0, trans_preframe);
  auto PreFrame_lv = new G4LogicalVolume(PreFrameSolid, m_material_map["Aluminum"],
					 "PreFrameLV");
  PreFrame_lv->SetVisAttributes(G4Colour::Red());

  // Bar
  auto BarMainSolid = new G4Box("BarMainSolid", 45.*mm, 618.25*mm, 5.*mm);
  auto BarSideSolid = new G4Box("BarSideSolid", 5.*mm, 618.25*mm, 5.*mm);
  G4ThreeVector trans_bar_side(-40.*mm, 0.*mm, -10.*mm);
  auto BarUniSolid = new G4UnionSolid("BarUniSolid", BarMainSolid, BarSideSolid,
                                      0, trans_bar_side);

  auto BarBotSolid = new G4Box("BarSideSolid", 45.*mm, 5.*mm, 5.*mm);
  G4ThreeVector trans_bar_bot(0.*mm, -613.25*mm, -10.*mm);

  auto BarSolid = new G4UnionSolid("BarSolid", BarUniSolid, BarBotSolid,
                                   0, trans_bar_bot);
  auto Bar_lv = new G4LogicalVolume(BarSolid, m_material_map["Iron"],
                                    "BarLV");
  Bar_lv->SetVisAttributes(G4Colour::Blue());

  //if(true){   // Htof Frame Placement
  if(false){   // Htof Frame Placement

    auto rotMOutRing = new G4RotationMatrix;
    rotMOutRing->rotateX(- 90. *deg);

    // Top Ring
    G4ThreeVector TopRing_pos(0.*mm, 586.72*mm, 0.*mm);
    TopRing_pos += htof_pos;
    new G4PVPlacement(rotMOutRing, TopRing_pos, TopRing_lv, "TopRingPV", m_world_lv, false, 0, m_check_overlaps);

    // Bottom Ring
    G4ThreeVector BotRing_pos(0.*mm, -586.72*mm, 0.*mm);
    BotRing_pos += htof_pos;
    new G4PVPlacement(rotMOutRing, BotRing_pos, BotRing_lv, "BotRingPV",
                      m_world_lv, false, 0, m_check_overlaps);

    // Preamp Support Frame
    for(G4int i=0; i<2; ++i){
      auto rotMOutP_PreFrame = new G4RotationMatrix;
      rotMOutP_PreFrame->rotateX((1 - 2 * i) * 90.*deg);
      rotMOutP_PreFrame->rotateZ(i * 180.*deg);

      G4ThreeVector PreFrame_pos(0.*mm, 453.22*mm, 0.*mm);
      PreFrame_pos.rotateZ(- i * 180.*deg);
      PreFrame_pos += htof_pos;
      new G4PVPlacement(rotMOutP_PreFrame, PreFrame_pos, PreFrame_lv,
                        Form("PreFramePV%d", i),
                        m_world_lv, false, i, m_check_overlaps);
    }

    G4int Bar_seg=0;
    for(G4int i=0; i<NumOfPlaneHTOF; ++i){
      //Bar
      auto rotMOutP_bar = new G4RotationMatrix;
      rotMOutP_bar->rotateY(- i * 45.*deg);
      G4ThreeVector Bar_pos(0.*mm, - 36.53*mm, - 372.13*mm);
      Bar_pos.rotateY(i * 45.*deg);
      Bar_pos += htof_pos;
      if(i!=0 && i!=4){ //Beam through
	new G4PVPlacement(rotMOutP_bar, Bar_pos, Bar_lv,
                          Form("BarPV%d", Bar_seg),
                          m_world_lv, false, Bar_seg, m_check_overlaps);
	Bar_seg++;
      }
      //Bracket
      for(G4int j=0; j<2; ++j){
	auto rotMOutP_bra = new G4RotationMatrix;
	rotMOutP_bra->rotateY( - (1 - 2 * j) * i * 45.*deg);
	rotMOutP_bra->rotateZ( - j * 180.*deg);
	G4ThreeVector Bra_pos(0.*mm, 506.72*mm, - 364.63*mm);
	Bra_pos.rotateY(i * 45.*deg);
	Bra_pos.rotateZ(j * 180.*deg);
	Bra_pos += htof_pos;
	G4int Bra_seg = 2 * i + j;
	new G4PVPlacement(rotMOutP_bra, Bra_pos, Bra_lv,
                          Form("BraPV%d", Bra_seg),
                          m_world_lv, false, Bra_seg, m_check_overlaps);
      }
    }
  }

}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructHypTPC()
{
  using CLHEP::mm;
  using CLHEP::deg;
  
  const auto tpc_pos = gGeom.GetGlobalPosition("HypTPC")*mm;
  const auto target_pos = gGeom.GetGlobalPosition("SHSTarget")*mm;
  const auto& gas_vessel_size = gSize.GetSize("TpcGasVessel")*mm*0.5;
  const auto& gas_vessel_window = gSize.GetSize("TpcGasVesselWindow")*mm*0.5;
  const auto& target_holder_size = gSize.GetSize("TargetHolder")*mm*0.5;
  const auto& target_holder_yofs = gas_vessel_size[2]-target_holder_size[2];
  const auto& target_holder_window = gSize.GetSize("TargetHolderWindowPhiDz")*mm;
  const auto target_holder_thickness =
    (target_holder_size[1] - target_holder_size[0])/2./2.; // 3*mm/2.
  const auto& p10_size = gSize.GetSize("TpcP10Volume")*mm/2.;
  const auto& field_cage_size = gSize.GetSize("TpcFieldCage")*mm/2.;
  const auto& eff_volume_size = gSize.GetSize("TpcEffectiveVolume")*mm/2.;
  const G4double phiStart = 22.5*deg;
  const G4double phiTotal = 360*deg;
  const G4int numSide   = 8;
  const G4int numZPlane = 2;
  // Gas Vessel
  const G4double zPlaneGV[] = { -gas_vessel_size.z(),
                                gas_vessel_size.z() };
  const G4double rInnerGV[] = { 0., 0. };
  const G4double rOuterGV[] = { gas_vessel_size[1],
                                gas_vessel_size[1] };
  G4VSolid* gv_solid;
  gv_solid = new G4Polyhedra("TpcGasVesselSolid",
                             phiStart, phiTotal, numSide, numZPlane,
                             zPlaneGV, rInnerGV, rOuterGV);
  auto space_solid = new G4Tubs("TpcTargetSpaceSolid",
                                0.,
                                target_holder_size[0],
                                target_holder_size[2],
                                phiStart, phiTotal);
  G4ThreeVector pos = target_pos;
  pos.setY(target_holder_yofs + target_holder_thickness*2.);
  pos.rotateX(90.*deg);
  gv_solid = new G4SubtractionSolid("TpcGasVesselSolid",
                                    gv_solid, space_solid,
                                    nullptr, pos);
  auto rot = new G4RotationMatrix;
  rot->rotateX(90.*deg);
  auto gv_lv = new G4LogicalVolume(gv_solid,
                                   m_material_map["Aluminum"],
                                   "TpcGasVesselLV");
  new G4PVPlacement(rot, tpc_pos, gv_lv, "TpcGasVesselPV",
                    m_world_lv, false, 0, m_check_overlaps);
  gv_lv->SetVisAttributes(G4Colour::White());
  // Target holder
  auto th_side_solid = new G4Tubs("TargetHolderSideSolid",
                                  target_holder_size[0],
                                  target_holder_size[1],
                                  target_holder_size[2],
                                  phiStart, phiTotal);
  auto th_bottom_solid = new G4Tubs("TargetHolderBottomSolid",
                                    0.,
                                    target_holder_size[0],
                                    target_holder_thickness,
                                    phiStart, phiTotal);
  pos.set(0, 0, -target_holder_size[2] + target_holder_thickness);
  auto th_solid = new G4UnionSolid("TargetholderSolid",
                                   th_side_solid, th_bottom_solid,
                                   nullptr, pos);
  auto th_lv = new G4LogicalVolume(th_solid,
                                   m_material_map["G10"],
                                   "TargetHolderLV");
  pos = target_pos;
  pos.setY(target_holder_yofs);
  pos.rotateX(90.*deg);
  new G4PVPlacement(nullptr, pos, th_lv, "TargetHolderPV",
                    gv_lv, false, 0, m_check_overlaps);
  th_lv->SetVisAttributes(G4Colour::Green());
  // Target holder window
  auto th_window_solid = new G4Tubs("TargetHolderWindowSolid",
                                    target_holder_size[0],
                                    target_holder_size[1],
                                    target_holder_window[2]/2.,
                                    target_holder_window[0]*deg,
                                    target_holder_window[1]*deg);
  auto th_window_lv = new G4LogicalVolume(th_window_solid,
                                          m_material_map["P10"],
                                          "TargetHolderWindowLV");
  rot = new G4RotationMatrix;
  rot->rotateZ(180.*deg);
  pos.set(0, 0, -target_holder_yofs);
  new G4PVPlacement(nullptr, pos, th_window_lv,
                    "TargetHolderWindowPV",
                    th_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(rot, pos, th_window_lv,
                    "TargetHolderWindowPV",
                    th_lv, false, 1, m_check_overlaps);
  th_window_lv->SetVisAttributes(G4Colour::White());

  auto gv_window_solid = new G4Box("GasVesselWindowSolid",
                                   gas_vessel_window.x(),
                                   gas_vessel_window.y(),
                                   // gas_vessel_window.z()
                                   (gas_vessel_size[1] - p10_size[1])/2.
                                   );
  auto gv_window_lv = new G4LogicalVolume(gv_window_solid,
                                          m_material_map["P10"],
                                          "GasVesselWindowLV");
  gv_window_lv->SetVisAttributes(G4Colour::White());
  pos.set(0., (gas_vessel_size[1]+p10_size[1])/2., 0.);
  for (G4int i=0; i<numSide; ++i) {
    rot = new G4RotationMatrix;
    rot->rotateX(90.*deg);
    rot->rotateY((i+1)*360.*deg/numSide);
    pos.rotateZ(360.*deg/numSide);
    new G4PVPlacement(rot, pos, gv_window_lv,
                      "GasVesselWindowPV"+std::to_string(i),
                      gv_lv, false, i, m_check_overlaps);
  }
  // P10
  const G4double rInnerP10[] = { p10_size[0], p10_size[0] };
  const G4double rOuterP10[] = { p10_size[1], p10_size[1] };
  const G4double zPlaneP10[] = { -p10_size[2], p10_size[2] };
  G4VSolid* p10_solid;
  p10_solid = new G4Polyhedra("TpcP10Solid",
                              phiStart, phiTotal, numSide, numZPlane,
                              zPlaneP10, rInnerP10, rOuterP10);
  pos = target_pos;
  pos.setY(target_holder_yofs);
  pos.rotateX(90.*deg);
  p10_solid = new G4SubtractionSolid("TpcP10Solid",
                                     p10_solid, space_solid,
                                     nullptr, pos);
  p10_solid = new G4SubtractionSolid("TpcP10Solid",
                                     p10_solid, th_side_solid,
                                     nullptr, pos);
  auto p10_lv = new G4LogicalVolume(p10_solid,
                                    m_material_map["P10"],
                                    "TpcP10LV");
  new G4PVPlacement(nullptr, tpc_pos, p10_lv, "TpcP10PV",
                    gv_lv, false, 1000, m_check_overlaps); //haein
  p10_lv->SetVisAttributes(G4Colour::Yellow());
  // p10_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto tpc_sd = new TPCSD("TPC");
  auto tpcpad_sd = new TPCPadSD("TPCPad");
  auto tpcedep_sd = new TPCEdepSD("TPCEdep");
  if(!gConf.Get<G4bool>("TPCPadOn")){
    AddNewDetector(tpc_sd);
    p10_lv->SetSensitiveDetector(tpc_sd);
  }
  else if(gConf.Get<G4bool>("TPCPadOn")){
    AddNewDetector(tpc_sd);
    p10_lv->SetSensitiveDetector(tpc_sd);
    AddNewDetector(tpcpad_sd);
    AddNewDetector(tpcedep_sd);
  }

  // Field Cage
  const G4double rInnerFC[] = { field_cage_size[0], field_cage_size[0] };
  const G4double rOuterFC[] = { field_cage_size[1], field_cage_size[1] };
  const G4double zPlaneFC[] = { -field_cage_size[2], field_cage_size[2] };
  auto fc_solid = new G4Polyhedra("FieldCageSolid",
                                  phiStart, phiTotal, numSide, numZPlane,
                                  zPlaneFC, rInnerFC, rOuterFC);
  auto fc_lv = new G4LogicalVolume(fc_solid, m_material_map["G10"],
                                   "FieldCageLV");
  new G4PVPlacement(nullptr, G4ThreeVector(), fc_lv,
                    "FieldCagePV", p10_lv, false, 0, m_check_overlaps);
  fc_lv->SetVisAttributes(G4Colour::Green());
  // Virtual pads
  G4LogicalVolume* pad_lv[NumOfPadTPC];
  G4LogicalVolume* pad_edep_lv[NumOfPadTPC];
  [[maybe_unused]] G4double angle[NumOfPadTPC+1] = {};
  const G4double pad_center_z = gSize.Get("TpcPadCenterZ")*mm;
  G4double pad_in[NumOfPadTPC] = {};
  G4double pad_out[NumOfPadTPC] = {};

  G4double pad_in_edep[NumOfPadTPC] = {};
  G4double pad_out_edep[NumOfPadTPC] = {};

  
  G4double tpc_rad = 250;
  // out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  const G4double pad_length_in = gSize.Get("TpcPadLengthIn");
  const G4double pad_length_out = gSize.Get("TpcPadLengthOut");
  const G4double pad_gap = gSize.Get("TpcPadGap");
  const G4int pad_configure = gSize.Get("TpcPadConfigure");
  switch (pad_configure) {
  case 1:
    for (G4int i=0; i<NumOfPadTPC; ++i) {
      if (i<NumOfPadTPCIn) {
	pad_in[i]  = 10.+(pad_length_in+pad_gap)*i;
	pad_out[i] = 10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]   = 360.;
      } else {
	pad_in[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn);
	pad_out[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn) + pad_length_out;
	angle[i] =
	  -std::acos((math::Sq(pad_out[i]) +
                      math::Sq(pad_center_z) -
                      math::Sq(tpc_rad)) /
                     (2*pad_out[i]*pad_center_z))*math::Rad2Deg() + 180.;
      }
    }
    break;
  case 2:
    for (G4int i=0; i<NumOfPadTPC; ++i) {
      if (i<NumOfPadTPCIn) {
	pad_in[i]  = 10.+(pad_length_in+pad_gap)*i;
	pad_out[i] = 10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]   = 360.;
      } else {
	pad_in[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn);
	pad_out[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn) + pad_length_out;
      }
    }
    angle[10] = 180. - 155.35;
    angle[11] = 180. - 144.8;
    angle[12] = 180. - 138.;
    angle[13] = 180. - 116.73;
    angle[14] = 180. - 106.;
    angle[15] = 180. - 98.77;
    angle[16] = 180. - 94.29;
    angle[17] = 180. - 89.8;
    angle[18] = 180. - 87.18;
    angle[19] = 180. - 84.16;
    angle[20] = 180. - 81.48;
    angle[21] = 180. - 73.39;
    angle[22] = 180. - 65.51011;
    angle[23] = 180. - 60.19;
    angle[24] = 180. - 56.35239;
    angle[25] = 180. - 52.85;
    angle[26] = 180. - 50.14;
    angle[27] = 180. - 47.17;
    angle[28] = 180. - 41.24;
    angle[29] = 180. - 29.;
    angle[30] = 180. - 23.23;
    angle[31] = 180. - 18.69;
    break;
  case 3:
    //for tracking analysis
    //If you need the dE/dx information, it should be modified.
    //Thin sensitive detector is introduced.
    for (G4int i=0; i<NumOfPadTPC; ++i) {
      G4double pad_radius = padHelper::getRadius(i);
      pad_in[i] = pad_radius;
      pad_out[i] = pad_radius + 0.1*mm;
      /*
        double pad_halflength = padHelper::getLength(i)/2;
        pad_in[i] = pad_radius-pad_halflength;
        pad_out[i] = pad_radius + pad_halflength;
      */
      if (i<NumOfPadTPCIn) {
	angle[i]   = 360.;
      } else {
	angle[i] = padHelper::getsTheta(i);
      }
    }
    break;
  case 4:
    //for tracking analysis
    // position -> pad center
    // dE/dx -> energy deposition
    for (G4int i=0; i<NumOfPadTPC; ++i) {
      G4double pad_radius = padHelper::getRadius(i);
      G4double pad_halflength = padHelper::getLength(i)/2;
      pad_in[i] = pad_radius;
      pad_out[i] = pad_radius + 0.1*mm;
      pad_in_edep[i] = pad_radius - pad_halflength;
      pad_out_edep[i] = pad_radius + pad_halflength;
      
      if (i<NumOfPadTPCIn) {
	angle[i]   = 360.;
      } else {
	angle[i] = padHelper::getsTheta(i);
      }
    }

    break;
  case 5:
    //for tracking alaysis
    // to get dE/dx information -> short step
    for (G4int i=0;i<NumOfPadTPC; ++i) {
      G4double pad_radius = padHelper::getRadius(i);
      G4double pad_halflength = padHelper::getLength(i)/2;
      pad_in[i] = pad_radius - pad_halflength;
      pad_out[i] = pad_radius + pad_halflength;
      if (i<NumOfPadTPCIn) {
	angle[i]   = 360.;
      } else {
	angle[i] = padHelper::getsTheta(i);
      }
    }
    break;

  default:
    break;
  }
  // Inner Pads
  const G4double rInnerEA[] = { eff_volume_size[0], eff_volume_size[0] };
  const G4double rOuterEA[] = { eff_volume_size[1], eff_volume_size[1] };
  const G4double zPlaneEA[] = { -eff_volume_size[2], eff_volume_size[2] };
  auto eff_volume = new G4Polyhedra("TpcEffectiveVolumeSolid",
                                  phiStart, phiTotal, numSide, numZPlane,
                                  zPlaneEA, rInnerEA, rOuterEA);
  G4VSolid* pad_solid[NumOfPadTPC];
  G4VSolid* pad_edep_solid[NumOfPadTPC];
  pos.set(0, pad_center_z, 0);
  for (G4int i=0; i<NumOfPadTPCIn; ++i) {

    pad_solid[i] = new G4Tubs("TpcPadSolid"+std::to_string(i),
                              pad_in[i]*mm,
                              pad_out[i]*mm,
                              field_cage_size[2],
                              phiStart, phiTotal);
    
    
    pad_solid[i] = new G4IntersectionSolid("TpcPadSolid"+std::to_string(i),
                                           pad_solid[i], p10_solid,
                                           nullptr, pos);
    pad_solid[i] = new G4IntersectionSolid("TpcPadSolid"+std::to_string(i),
                                           pad_solid[i], eff_volume,
                                           nullptr, pos);
    pad_lv[i]  = new G4LogicalVolume(pad_solid[i], m_material_map["P10"],
                                     "TpcPadLV"+std::to_string(i));

    if(pad_configure!=4)
      new G4PVPlacement(nullptr, -pos, pad_lv[i], "TpcPadPV"+std::to_string(i),
			p10_lv, true, i, m_check_overlaps);

    if(pad_configure==4){
      pad_edep_solid[i] = new G4Tubs("TpcPadEdepSolid"+std::to_string(i),
				     pad_in_edep[i]*mm,
				     pad_out_edep[i]*mm,
				     field_cage_size[2],
				     phiStart, phiTotal);

   
      pad_edep_solid[i] = new G4IntersectionSolid("TpcPadEdepSolid"+std::to_string(i),
						  pad_edep_solid[i], p10_solid,
						  nullptr, pos);
      pad_edep_solid[i] = new G4IntersectionSolid("TpcPadEdepSolid"+std::to_string(i),
						  pad_edep_solid[i], eff_volume,
						  nullptr, pos);

      /*
	pad_edep_solid[i] = new G4SubtractionSolid("TpcPadEdepSolid"+std::to_string(i),
	pad_edep_solid[i], pad_solid[i],
	nullptr, G4ThreeVector());
      */
      pad_edep_lv[i]  = new G4LogicalVolume(pad_edep_solid[i], m_material_map["P10"],
					    "TpcPadEdepLV"+std::to_string(i));

    
    
      new G4PVPlacement(nullptr, -pos, pad_edep_lv[i], "TpcPadEdepPV"+std::to_string(i),
			p10_lv, true, 2000 + i, m_check_overlaps);
      
      new G4PVPlacement(nullptr, G4ThreeVector(), pad_lv[i], "TpcPadPV"+std::to_string(i),
			pad_edep_lv[i], true, i, m_check_overlaps);

    }
  }
  // Outer Pads
  
  for(G4int i=NumOfPadTPCIn; i<NumOfPadTPC; ++i){
    pad_solid[i] = new G4Tubs("TpcPadSolid"+std::to_string(i),
			      pad_in[i]*mm,
			      pad_out[i]*mm,
			      field_cage_size[2],
			      (90.+angle[i])*deg,
			      (360.-2.*angle[i])*deg);


    pad_solid[i] = new G4IntersectionSolid("TpcPadSolid"+std::to_string(i),
					   pad_solid[i], p10_solid,
					   nullptr, pos);
    pad_solid[i] = new G4IntersectionSolid("TpcPadSolid"+std::to_string(i),
					   pad_solid[i], eff_volume,
					   nullptr, pos);
    pad_lv[i]  = new G4LogicalVolume(pad_solid[i], m_material_map["P10"],
				     "TpcPadLV"+std::to_string(i));
    if(pad_configure!=4)
      new G4PVPlacement(nullptr, -pos, pad_lv[i], "TpcPadPV"+std::to_string(i),
			p10_lv, true, i, m_check_overlaps);

    if(pad_configure==4){
      pad_edep_solid[i] = new G4Tubs("TpcPadEdepSolid"+std::to_string(i),
                              pad_in_edep[i]*mm,
                              pad_out_edep[i]*mm,
                              field_cage_size[2],
			      (90.+angle[i])*deg,
			      (360.-2.*angle[i])*deg);

    

      pad_edep_solid[i] = new G4IntersectionSolid("TpcPadEdepSolid"+std::to_string(i),
						  pad_edep_solid[i], p10_solid,
						  nullptr, pos);
      pad_edep_solid[i] = new G4IntersectionSolid("TpcPadEdepSolid"+std::to_string(i),
						  pad_edep_solid[i], eff_volume,
						  nullptr, pos);
      /*pad_edep_solid[i] = new G4SubtractionSolid("TpcPadEdepSolid"+std::to_string(i),
	pad_edep_solid[i], pad_solid[i],
	nullptr, G4ThreeVector());
      */
      pad_edep_lv[i]  = new G4LogicalVolume(pad_edep_solid[i], m_material_map["P10"],
					    "TpcPadEdepLV"+std::to_string(i));
      new G4PVPlacement(nullptr, -pos, pad_edep_lv[i], "TpcPadEdepPV"+std::to_string(i),
			p10_lv, true, 2000 + i, m_check_overlaps);

      new G4PVPlacement(nullptr, G4ThreeVector(), pad_lv[i], "TpcPadPV"+std::to_string(i),
			pad_edep_lv[i], true, i, m_check_overlaps);
    }
    
  }
  

  for (G4int i=0; i<NumOfPadTPC; ++i) {
    pad_lv[i]->SetVisAttributes(ORANGE);
    //pad_edep_lv[i]->SetVisAttributes(G4VisAttributes::GetInvisible());
    if(!gConf.Get<G4bool>("TPCPadOn"))
      pad_lv[i]->SetSensitiveDetector(tpc_sd);
    else if(gConf.Get<G4bool>("TPCPadOn")){
      pad_lv[i]->SetSensitiveDetector(tpcpad_sd);
      if(pad_configure==4)
	pad_edep_lv[i]->SetSensitiveDetector(tpcedep_sd);
    }
    
  }
  // Dead area
  auto dead_solid = new G4Box("DeadSolid", 5*mm, 250*mm, 0.001*mm);
  auto dead_lv = new G4LogicalVolume(dead_solid, m_material_map["P10"],
                                     "DeadLV");
  auto rotdead1 = new G4RotationMatrix;
  rotdead1->rotateZ(45.*deg);
  new G4PVPlacement(rotdead1, G4ThreeVector(0., 0.*mm, -300.1*mm),
                    dead_lv, "DeadPV1", p10_lv, true, 0, m_check_overlaps);
  auto rotdead2 = new G4RotationMatrix;
  rotdead2->rotateZ(-45.*deg);
  new G4PVPlacement(rotdead2, G4ThreeVector(0., 0.*mm, -300.1*mm),
                    dead_lv, "DeadPV2", p10_lv, true, 1, m_check_overlaps);
  dead_lv->SetVisAttributes(G4Colour::Gray());
#if 0
  // Virtual pad
  G4Tubs* vpad_solid[NumOfPadTPC];
  G4LogicalVolume* vpad_lv[NumOfPadTPC];
  for(G4int i=0; i<NumOfPadTPCIn; ++i){
    vpad_solid[i] = new G4Tubs(Form("TpcVPadSolid%d", i), pad_in[i]*mm,
                               pad_out[i]*mm, 0.5*mm, 0., 360.*deg);
    vpad_lv[i] = new G4LogicalVolume(vpad_solid[i], m_material_map["P10"],
                                     Form("TpcVPadLV%d", i));
    vpad_lv[i]->SetVisAttributes(ORANGE);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., -pad_center_z, -302.*mm),
                      vpad_lv[i], Form("TpcVPadPV%d", i), p10_lv, true, 0, m_check_overlaps);
  }
  for(G4int i=NumOfPadTPCIn; i<NumOfPadTPC; ++i){
    vpad_solid[i] = new G4Tubs(Form("TpcVPadSolid%d", i),  pad_in[i]*mm,
                               pad_out[i]*mm,
                               0.5*mm, (90.+angle[i])*deg,
                               (360.-2.*angle[i])*deg);
    vpad_lv[i] = new G4LogicalVolume(vpad_solid[i], m_material_map["P10"],
                                     Form("TpcVPadLV%d", i));
    vpad_lv[i]->SetVisAttributes(ORANGE);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., -pad_center_z, -302.*mm),
                      vpad_lv[i], Form("TpcVPadPV%d", i), p10_lv, true, 0, m_check_overlaps);
  }
#endif
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructKVC()
{
  using CLHEP::mm;
  using CLHEP::deg;
  const auto& ra2 = gGeom.GetRotAngle2("KVC")*deg;
  const auto& half_size = gSize.GetSize("KvcRadiator")*mm/2.;
  auto pos = gGeom.GetGlobalPosition("KVC");
  auto kvcSD = new KVCSD("KVC");
  AddNewDetector(kvcSD);
  auto mother_solid = new G4Box("KvcMotherSolid",
                                half_size.x()*NumOfSegKVC + 10*mm,
                                half_size.y() + 10*mm,
                                half_size.z() + 10*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "KvcMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(- ra2 - m_rotation_angle);
  pos.rotateY(m_rotation_angle);
  new G4PVPlacement(rot, pos, mother_lv,
                    "KvcMotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Radiator
  auto kvc_solid = new G4Box("KvcSolid",
                             half_size.x(), half_size.y(), half_size.z());
auto kvc_lv = new G4LogicalVolume(kvc_solid, m_material_map["QuartzKVC"], "KvcLV");
  for(G4int i=0; i<NumOfSegKVC; ++i){
    pos = G4ThreeVector(half_size.x()*(-NumOfSegKVC+1+2*i), 0., 0.);
    new G4PVPlacement(nullptr, pos, kvc_lv, "KvcPV",
                      mother_lv, false, i, m_check_overlaps);
  }
  kvc_lv->SetVisAttributes(G4Colour::Yellow());
  kvc_lv->SetSensitiveDetector(kvcSD);
}
//_____________________________________________________________________________ BVH 추가
void DetectorConstruction::ConstructBVH_U()
{
  using CLHEP::mm;
  using CLHEP::deg;
  const auto& ra2 = gGeom.GetRotAngle2("BVH_U")*deg;
  const auto& half_size = gSize.GetSize("BvhURadiator")*mm/2.;
  auto pos = gGeom.GetGlobalPosition("BVH_U");
  auto bvhUSD = new BVH_U_SD("BVH_U"); // 새로운 SD 클래스 만들거나 KVCSD 재사용
  AddNewDetector(bvhUSD);

  auto mother_solid = new G4Box("BvhUMotherSolid",
                                half_size.x()*NumOfSegBVH_U + 10*mm,
                                half_size.y() + 10*mm,
                                half_size.z() + 10*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "BvhUMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(- ra2 - m_rotation_angle);
  pos.rotateY(m_rotation_angle);
  new G4PVPlacement(rot, pos, mother_lv,
                    "BvhUMotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Radiator
  auto bvh_solid = new G4Box("BvhUSolid",
                             half_size.x(), half_size.y(), half_size.z());
  auto bvh_lv = new G4LogicalVolume(bvh_solid, m_material_map["EJ232"], "BvhULV");

  for(G4int i=0; i<NumOfSegBVH_U; ++i){
    pos = G4ThreeVector(half_size.x()*(-NumOfSegBVH_U+1+2*i), 0., 0.);
    new G4PVPlacement(nullptr, pos, bvh_lv, "BvhUPV",
                      mother_lv, false, i, m_check_overlaps);
  }
  bvh_lv->SetVisAttributes(G4Colour::Green());
  bvh_lv->SetSensitiveDetector(bvhUSD);
}
//_____________________________________________________________________________ BVH 추가
void DetectorConstruction::ConstructBVH_D()
{
  using CLHEP::mm;
  using CLHEP::deg;
  const auto& ra2 = gGeom.GetRotAngle2("BVH_D")*deg;
  const auto& half_size = gSize.GetSize("BvhDRadiator")*mm/2.;
  auto pos = gGeom.GetGlobalPosition("BVH_D");
  auto bvhDSD = new BVH_D_SD("BVH_D");
  AddNewDetector(bvhDSD);
  auto mother_solid = new G4Box("BvhDMotherSolid",
                                half_size.x()*NumOfSegBVH_D + 10*mm,
                                half_size.y() + 10*mm,
                                half_size.z() + 10*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "BvhDMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(- ra2 - m_rotation_angle);
  pos.rotateY(m_rotation_angle);
  new G4PVPlacement(rot, pos, mother_lv,
                    "BvhDMotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto bvh_solid = new G4Box("BvhDSolid",
                             half_size.x(), half_size.y(), half_size.z());
  auto bvh_lv = new G4LogicalVolume(bvh_solid, m_material_map["EJ232"], "BvhDLV");
  for(G4int i=0; i<NumOfSegBVH_D; ++i){
    pos = G4ThreeVector(half_size.x()*(-NumOfSegBVH_D+1+2*i), 0., 0.);
    new G4PVPlacement(nullptr, pos, bvh_lv, "BvhDPV",
                      mother_lv, false, i, m_check_overlaps);
  }
  bvh_lv->SetVisAttributes(G4Colour::Green());
  bvh_lv->SetSensitiveDetector(bvhDSD);
}
//_____________________________________________________________________________
void DetectorConstruction::ConstructT0()
{
  using CLHEP::mm;
  using CLHEP::deg;

  const auto& ra2 = gGeom.GetRotAngle2("T0") * deg;
  const auto& half_size = gSize.GetSize("T0Radiator") * mm / 2.;
  auto pos = gGeom.GetGlobalPosition("T0");

  auto t0SD = new T0SD("T0"); 
  AddNewDetector(t0SD);

  // 모체 볼륨
  auto mother_solid = new G4Box("T0MotherSolid",
                                half_size.x() * NumOfSegT0 + 10 * mm,
                                half_size.y() + 10 * mm,
                                half_size.z() + 10 * mm);

  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "T0MotherLV");

auto rot = new G4RotationMatrix;
rot->rotateZ(45 * deg); // 🌀 z축 기준 반시계 방향 45도 회전
rot->rotateY(- ra2 - m_rotation_angle); // 기존 설정 유지
pos.rotateY(m_rotation_angle);

  new G4PVPlacement(rot, pos, mother_lv,
                    "T0MotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Radiator (Scintillator)
  auto t0_solid = new G4Box("T0Solid",
                            half_size.x(), half_size.y(), half_size.z());

  auto t0_lv = new G4LogicalVolume(t0_solid,
                                   m_material_map["EJ232"], // EJ232로 가정
                                   "T0LV");

  for (G4int i = 0; i < NumOfSegT0; ++i) {
    pos = G4ThreeVector(half_size.x() * (-NumOfSegT0 + 1 + 2*i), 0., 0.);
    new G4PVPlacement(nullptr, pos, t0_lv, "T0PV",
                      mother_lv, false, i, m_check_overlaps);
  }

  t0_lv->SetVisAttributes(G4Colour::Blue());
  t0_lv->SetSensitiveDetector(t0SD);
}
//-----------------------------------------------------------------------------
void DetectorConstruction::ConstructSCH()
{
  using CLHEP::mm;

  const auto& sch_pos  = gGeom.GetGlobalPosition("SCH");
  const auto& sch_size = gSize.GetSize("SchSeg") * 0.5 * mm;
  const G4double dXdW  = gGeom.GetWirePitch("SCH") * mm;

  // 1) SD 이름: "/SCH" -> "SCH"
  auto schSD = new SCHSD("SCH");
  // 2) 내부 등록: 직접 SDManager 호출(X) -> AddNewDetector(O)
  AddNewDetector(schSD);

  G4LogicalVolume* sch_lv[NumOfSegSCH];
  auto sch_solid = new G4Box("SchSolid", sch_size.x(), sch_size.y(), sch_size.z());

  for(G4int i=0; i<NumOfSegSCH; ++i){
    sch_lv[i] = new G4LogicalVolume(sch_solid, m_material_map["Scintillator"],
                                    Form("SchSeg%dLV", i), 0, 0, 0);
    sch_lv[i]->SetVisAttributes(G4Colour::Cyan());

    const G4double ipos_x = dXdW * (i - (NumOfSegSCH - 1)/2.);
    G4ThreeVector pos(sch_pos.x() + ipos_x, sch_pos.y(),
                      sch_pos.z() + 1.*mm*(1 - 2*(i%2)));

    new G4PVPlacement(m_rotation_matrix, pos, sch_lv[i],
                      Form("SchSeg%dPV", i), m_world_lv, false, i, m_check_overlaps);

    sch_lv[i]->SetSensitiveDetector(schSD);
  }
}



//_____________________________________________________________________________
void
DetectorConstruction::ConstructShsMagnet()
{
  using CLHEP::mm;
  using CLHEP::deg;
#if 0 // Old version
  const auto& tpc_pos = gGeom.GetGlobalPosition("HypTPC");
  const G4double DPHI_TPC = 360.*deg;
  auto tube_solid = new G4Tubs("TubeSolid", 45.*cm, 80.*cm, 135./2.*cm, 0.*deg, DPHI_TPC);
  auto hole_solid = new G4Box("HoleSolid", 100./2.*cm, 100./2.*cm, 60./2.*cm);
  G4ThreeVector HTrans(0, -60.*cm, 0);
  G4RotationMatrix* yRot = new G4RotationMatrix;  // Rotates X and Z axes only
  G4Transform3D transform(*yRot, HTrans);
  auto coil_solid = new G4SubtractionSolid("ShsMagnetCoilSolid",
                                           tube_solid, hole_solid, transform);
  auto coil_lv = new G4LogicalVolume(coil_solid, m_material_map["Iron"],
                                     "ShsMagnetCoilLV");
  G4RotationMatrix *rotHelm = new G4RotationMatrix();
  rotHelm->rotateX(90.*deg);
  rotHelm->rotateZ(- m_rotation_angle);
  new G4PVPlacement(rotHelm, tpc_pos, coil_lv, "ShsMagnetCoilPV",
                    m_world_lv, false, 0, m_check_overlaps);
  coil_lv->SetVisAttributes(PINK);
#else // Current version
  const G4ThreeVector yoke_size(1550./2.*mm, 950./2.*mm, 1200./2.*mm);
  // auto shs_solid = new G4Box("ShsMagnetSolid", yoke_size.x(),
  // 			      yoke_size.y(), yoke_size.z());
  // auto shs_lv = new G4LogicalVolume(shs_solid, m_material_map["Air"],
  // 				     "ShsMagnetLV");
  // shs_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // auto shs_pv = new G4PVPlacement(nullptr, tpc_pos, shs_lv,
  // 				   "ShsMagnetPV", m_world_lv, false, 0, m_check_overlaps);
  const G4int NumOfParams = 4;
  G4double yoke_width[NumOfParams] = { 1550*mm, 1530*mm, 1480*mm, 1470*mm };
  G4double yoke_depth[NumOfParams] = { 1200*mm, 1180*mm, 1140*mm, 1130*mm };
  G4double yoke_height[NumOfParams] = { 950*mm, 890*mm, 732*mm, 712*mm };
  G4double coner_size[NumOfParams] = { 1449.569*mm, 1429.569*mm,
				       1357.645*mm, 1347.645*mm };
  G4double chamber_rad[NumOfParams] = { 800./2.*mm, 820./2.*mm,
					850./2.*mm, 860./2.*mm };
  G4double hole_gap[NumOfParams] = { 300*mm, 320*mm, 348*mm, 358*mm };
  G4double side_gap[NumOfParams];
  for(G4int i=0 ; i<NumOfParams; ++i){
    side_gap[i] = yoke_depth[i]/2.0*std::sqrt(2.)*mm;
  }
  G4double dummy_size[3] = { 2000*mm, 2000*mm, 2000*mm };
  G4double hole_pos[NumOfParams];
  for(G4int i=0; i<NumOfParams; ++i){
    hole_pos[i] = yoke_depth[i]/2.;
  }
  G4Box* box_solid[NumOfParams];
  G4Tubs* tube_solid[NumOfParams];
  G4Box* box_outer_solid[NumOfParams];
  G4Box* box_inner_solid[NumOfParams];
  G4Box* box_side_solid[NumOfParams];
  G4SubtractionSolid* octa_solid[NumOfParams];
  G4SubtractionSolid* corner_solid[NumOfParams];
  G4SubtractionSolid* side1_solid[NumOfParams];
  G4SubtractionSolid* side2_solid[NumOfParams];
  G4SubtractionSolid* chamber_solid[NumOfParams];
  G4RotationMatrix* rot_magnet  = new G4RotationMatrix;
  rot_magnet->rotateZ(45.*deg);
  for(G4int i=0; i<NumOfParams; ++i){
    box_solid[i] = new G4Box(Form("Box%dSolid", i), yoke_width[i]/2.,
                             yoke_depth[i]/2., yoke_height[i]/2.);
    box_outer_solid[i] = new G4Box(Form("BoxOuter%dSolid", i),
                                   dummy_size[0], dummy_size[1],
                                   dummy_size[2]);
    box_inner_solid[i] = new G4Box(Form("BoxInner%dSolid", i),
                                   coner_size[i], coner_size[i],
                                   dummy_size[2]);
    box_side_solid[i] = new G4Box(Form("BoxSide%dSolid", i),
                                  side_gap[i]/2., side_gap[i]/2.,
                                  hole_gap[i]/2.);
    tube_solid[i] = new G4Tubs(Form("Tube%dSolid", i), 0., chamber_rad[i],
                               2000*mm, 0.*deg, 360.*deg);
    corner_solid[i] = new G4SubtractionSolid(Form("Corner%dSolid", i),
                                             box_outer_solid[i],
                                             box_inner_solid[i],
                                             nullptr, G4ThreeVector());
    octa_solid[i] = new G4SubtractionSolid(Form("Octa%dSolid", i),
                                           box_solid[i], corner_solid[i],
                                           rot_magnet, G4ThreeVector());
    G4ThreeVector side1_size(0., hole_pos[i], 0.);
    side1_solid[i] = new G4SubtractionSolid(Form("Side1-%dSolid", i),
                                            octa_solid[i], box_side_solid[i],
                                            rot_magnet, side1_size);
    G4ThreeVector side2_size(0., -hole_pos[i], 0.);
    side2_solid[i] = new G4SubtractionSolid(Form("Side2-%dSolid", i),
                                            side1_solid[i], box_side_solid[i],
                                            rot_magnet, side2_size);
    chamber_solid[i] = new G4SubtractionSolid(Form("Chamber%dSolid", i),
                                              side2_solid[i], tube_solid[i],
                                              nullptr, G4ThreeVector());
  }
  auto frame1_solid = new G4SubtractionSolid("Frame1Solid", chamber_solid[0],
                                             chamber_solid[1], nullptr,
                                             G4ThreeVector());
  auto frame2_solid = new G4SubtractionSolid("Frame2Solid", chamber_solid[2],
                                             chamber_solid[3], nullptr,
                                             G4ThreeVector());
  auto magnet_solid = new G4UnionSolid("ShsMagnetSolid", frame1_solid,
                                       frame2_solid, nullptr,
                                       G4ThreeVector());
  auto magnet_lv = new G4LogicalVolume(magnet_solid,
                                       m_material_map["Iron"],
                                       "ShsMagnetLV");
  // magnet_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  magnet_lv->SetVisAttributes(PINK);
  G4RotationMatrix rot_frame;
  rot_frame.rotateX(90.*deg);
  new G4PVPlacement(G4Transform3D(rot_frame, G4ThreeVector()),
                    magnet_lv, "ShsMagnetPV", m_world_lv, false, 0, m_check_overlaps);
  // Coil Support
  G4double CoilSupPos_height = 250*mm;
  G4double RadIn = 445*mm;
  G4double RadOut= 545*mm;
  G4double GapRadIn = 465*mm;
  G4double GapRadOut = 545*mm;
  G4double SupHeight = 112*mm;
  G4double GapHeight = 62*mm;
  G4String fullNameCoilSup = "SCCoilSup";
  auto solidTube_Sup = new G4Tubs("CoilSupMainSolid", RadIn, RadOut, SupHeight/2.,
                                  0*deg, 360*deg);
  auto solidTube_Sub = new G4Tubs("CoilSupSubSolid", GapRadIn, GapRadOut+20*mm,
                                  GapHeight/2., 0*deg, 360*deg);
  auto solidDetectorCS = new G4SubtractionSolid("CoilSupSolid", solidTube_Sup,
                                                solidTube_Sub, nullptr,
                                                G4ThreeVector());
  auto logicDetectorCS = new G4LogicalVolume(solidDetectorCS,
                                             m_material_map["Iron"],
                                             "CoilSupLV");
  // logicDetectorCS->SetVisAttributes(G4Color::Green());
  logicDetectorCS->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4RotationMatrix rot_sup;
  rot_sup.rotateX(90.*deg);
  new G4PVPlacement(G4Transform3D(rot_sup, G4ThreeVector(0, CoilSupPos_height, 0)),
                    logicDetectorCS, "CoilSupUpPV", m_world_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(rot_sup, G4ThreeVector(0, -CoilSupPos_height, 0)),
                    logicDetectorCS, "CoilSupDwPV", m_world_lv, false, 0, m_check_overlaps);
  // Coil
  const G4double coilRad_in = 466*mm;
  const G4double coilRad_out = 535*mm;
  const G4double coilHeight = 60*mm;
  const G4ThreeVector coilu_pos(0*mm, 250*mm, 0*mm);
  const G4ThreeVector coild_pos(0*mm, -250*mm, 0*mm);
  auto solidDetectorCoil = new G4Tubs("CoilSolid", coilRad_in, coilRad_out,
                                      coilHeight/2., 0*deg, 360*deg);
  auto logicDetectorCoil = new G4LogicalVolume(solidDetectorCoil,
                                               m_material_map["Copper"],
                                               "CoilLV");
  // logicDetectorCoil->SetVisAttributes(G4Color::Yellow());
  logicDetectorCoil->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4RotationMatrix rot_coil;
  rot_coil.rotateX(90.*deg);
  new G4PVPlacement(G4Transform3D(rot_coil, coilu_pos),
                    logicDetectorCoil, "CoilUpPV", m_world_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(G4Transform3D(rot_coil, coild_pos),
                    logicDetectorCoil, "CoilDwPV", m_world_lv, false, 0, m_check_overlaps);
#endif
  m_field->SetStatusShsField(true);
  m_field->SetShsFieldMap(gConf.Get<G4String>("SHSFLDMAP"));
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructTarget()
{
  using CLHEP::mm;
  using CLHEP::deg;
  auto target_sd = new TargetSD("TGT");
  AddNewDetector(target_sd);
  const auto target_pos = gGeom.GetGlobalPosition("SHSTarget")*mm;
  const auto target_size = gSize.GetSize("Target")*mm/2.;
  G4VSolid* target_solid;
  auto rot = new G4RotationMatrix;
  switch (m_experiment) {
  case 42: {
    target_solid = new G4Box("Target", target_size.x(),
                             target_size.y(), target_size.z());
  }
    break;
  case 45: case 27: {
    G4double target_r = gSize.Get("Target", ThreeVector::X);
    G4double target_z = gSize.Get("Target", ThreeVector::Z);
    target_solid = new G4Tubs("TargetSolid", 0.*mm,
                              target_r*mm, target_z*mm, 0., 360*deg);
    rot->rotateX(90.*deg);
    break;
  }
  case 72: {
    if (gConf.Get<G4bool>("TargetVP")) {
      target_solid = new G4Box("TargetSolid",
                               target_size[1],
                               0.001*mm,
                               target_size[2]);
    } else {
      target_solid = new G4Tubs("TargetSolid",
                                target_size[0],
                                target_size[1],
                                target_size[2],
                                0.*deg, 360.*deg);
    }
    rot->rotateX(90.*deg);
    const auto kapton_size = gSize.GetSize("TargetKapton")*mm/2.;
    const auto gfrp_size = gSize.GetSize("TargetGFRP")*mm/2.;
    auto kapton = new G4Tubs("TargetKapton",
                             kapton_size[0],
                             kapton_size[1],
                             kapton_size[2],
                             0.*deg, 360.*deg);
    auto kapton_lv = new G4LogicalVolume(kapton,
                                         m_material_map["Kapton"],
                                         "TargetKaptonLV");
    kapton_lv->SetVisAttributes(G4Colour::Red());
    new G4PVPlacement(rot, target_pos, kapton_lv, "TargetKaptonPV",
                      m_world_lv, true, 0, m_check_overlaps);
    auto gfrp = new G4Tubs("TargetGFRP",
                           gfrp_size[0],
                           gfrp_size[1],
                           gfrp_size[2],
                           0.*deg, 360.*deg);
    auto gfrp_lv = new G4LogicalVolume(gfrp,
                                       m_material_map["G10"],
                                       "TargetGFRPLV");
    gfrp_lv->SetVisAttributes(G4Colour::Green());
    new G4PVPlacement(rot, target_pos, gfrp_lv, "TargetGFRPPV",
                      m_world_lv, true, 0, m_check_overlaps);
  }
    break;
  default:
    G4Exception(FUNC_NAME,
                "Invalid experiment", FatalException,
                ("Found invalid experiment "+std::to_string(m_experiment)
                 +" in "+gConf.Get<G4String>("CONF")).c_str());
    return;
  }
  auto target_lv = new G4LogicalVolume(target_solid,
                                       m_material_map["Target"],
                                       "TargetLV");
  target_lv->SetSensitiveDetector(target_sd);
  target_lv->SetVisAttributes(G4Colour::Blue());
  new G4PVPlacement(rot, target_pos,
                    target_lv, "TargetPV",
                    m_world_lv, true, 0, m_check_overlaps);
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructFieldOutline()
{
  using CLHEP::mm;
  using CLHEP::deg;
  if(m_field->GetStatusShsField()){
    const auto ra2 = 0.*deg;
    const auto half_size = m_field->GetSizeShsField()/2.;
    if (half_size.mag() == 0.) return;
    auto pos = G4ThreeVector();
    auto mother_solid = new G4Box("FieldOutlineMotherSolid",
                                  half_size.x() + 1*mm,
                                  half_size.y() + 1*mm,
                                  half_size.z() + 1*mm);
    auto mother_lv = new G4LogicalVolume(mother_solid,
                                         m_material_map["Air"],
                                         "FieldOutlineMotherLV");
    auto rot = new G4RotationMatrix;
    rot->rotateY(- ra2 - m_rotation_angle);
    pos.rotateY(m_rotation_angle);
    new G4PVPlacement(rot, pos, mother_lv,
                      "FieldOutlineMotherPV", m_world_lv, false, 0, m_check_overlaps);
    mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
    auto solid = new G4Box("FieldOutlineSolid",
                           half_size.x(), half_size.y(), half_size.z());
    auto lv = new G4LogicalVolume(solid, m_material_map["Air"],
                                  "FieldOutlineLV");
    new G4PVPlacement(nullptr, G4ThreeVector(), lv, "FieldOutlinePV",
                      mother_lv, false, 0, m_check_overlaps);
    lv->SetVisAttributes(G4Colour::Yellow());
  }
}

//_____________________________________________________________________________
// DetectorConstruction::ConstructVP()

void DetectorConstruction::ConstructVP()
{
  using CLHEP::mm;
  using CLHEP::deg;

  // 1) 여기서 ‘활성화’할 VP 번호만 골라주세요.
  //    지금은 VP4, VP5, VP6만 생성 (VP1~3 비활성화)
  const std::vector<int> enabledVP = {1, 2, 6, 7, 8};

  for (int i : enabledVP) {
    try {
      G4String name = "VP" + std::to_string(i);
      const auto& ra2 = gGeom.GetRotAngle2(name)*deg;
      const auto& half_size = gSize.GetSize(name)*mm/2.;
      auto pos = gGeom.GetGlobalPosition(name);

      // (중요) mother 여유를 줄여 겹침 여지 최소화 (1mm -> 0.1mm)
      auto mother_solid = new G4Box(name+"MotherSolid",
                                    half_size.x() + 0.1*mm,
                                    half_size.y() + 0.1*mm,
                                    half_size.z() + 0.1*mm);
      auto mother_lv = new G4LogicalVolume(mother_solid,
                                           m_material_map["Air"],
                                           name+"MotherLV");
      auto rot = new G4RotationMatrix;
      rot->rotateY(- ra2 - m_rotation_angle);
      pos.rotateY(m_rotation_angle);

      new G4PVPlacement(rot, pos, mother_lv, name+"MotherPV",
                        m_world_lv, false, 0, m_check_overlaps);
      mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

      auto vp_solid = new G4Box(name+"Solid", half_size.x(), half_size.y(), half_size.z());
      auto vp_lv = new G4LogicalVolume(vp_solid, m_material_map["Air"], name+"LV");
      new G4PVPlacement(nullptr, G4ThreeVector(), vp_lv, name+"PV",
                        mother_lv, false, i, m_check_overlaps);
      vp_lv->SetVisAttributes(G4Colour::Yellow());

      static auto vpSD = new VPSD("VP");
      AddNewDetector(vpSD);
      vp_lv->SetSensitiveDetector(vpSD);

    } catch (const std::exception& e) {
      // 해당 VP가 DCGeomParam에 없으면 그냥 건너뜀
      continue;
    }
  }
} 
