// -*- C++ -*-

#include "DetectorConstruction.hh"

#include <G4Box.hh>
#include <G4ChordFinder.hh>
#include <G4Element.hh>
#include <G4FieldManager.hh>
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
#include "LACSD.hh"
#include "MagneticField.hh"
#include "TPCSD.hh"
#include "SCHSD.hh"
#include "SDCSD.hh"
#include "TargetSD.hh"
#include "VPSD.hh"
#include "WCSD.hh"
#include "padHelper.hh"

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

SDCSD* sdc_sd = nullptr;
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructAreaTent()
{
  using CLHEP::m;
  using CLHEP::mm;
  //const G4ThreeVector size(5.0*m/2, 5.0*m/2, 6.0*m/2);
  const G4ThreeVector size(4.8*m/2, 4.8*m/2, 6.3*m/2);
  //const G4ThreeVector pos(0., 0., -143.-1200.-170.*mm+size.z());
  const G4ThreeVector pos(0., 0., -1270.-119.*mm+size.z());
  auto tent_out_solid = new G4Box("AreaTentOutSolid",
                                  size.x(), size.y(), size.z());
  auto tent_in_solid = new G4Box("AreaTentInSolid", size.x()-1.*mm,
                                 size.y()-1.*mm, size.z()-1.*mm);
  auto tent_solid = new G4SubtractionSolid("AreaTentSolid",
                                           tent_out_solid, tent_in_solid);
  auto tent_lv = new G4LogicalVolume(tent_solid,
                                     m_material_map["Air"],
                                     "AreaTentLV");
  tent_lv->SetVisAttributes(G4Color::Blue());
  // tent_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  new G4PVPlacement(nullptr, pos, tent_lv,
                    "AreaTentPV", m_world_lv, false, 0);
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructBC3()
{
  using CLHEP::mm;
  if(!sdc_sd){
    sdc_sd = new SDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector(sdc_sd);
  }

  const auto& bc3_pos = (gGeom.GetGlobalPosition("KURAMA") +
                         (gGeom.GetGlobalPosition("BC3-X1") +
                          gGeom.GetGlobalPosition("BC3-U2")) * 0.5);
  const auto& frame_size = gSize.GetSize("Bc3Frame") * 0.5 * mm;
  const auto& drift_size = gSize.GetSize("Bc3Drift") * 0.5 * mm;
  auto bc3_solid = new G4Box("Bc3Solid", frame_size.x(),
                             frame_size.y(), frame_size.z());
  auto bc3_lv = new G4LogicalVolume(bc3_solid, m_material_map["Argon"],
                                    "Bc3LV", 0, 0, 0);
  bc3_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(m_rotation_matrix, bc3_pos,
                    bc3_lv, "Bc3PV", m_world_lv, false, 0);
  auto bc3pl_solid = new G4Box("Bc3PlSolid", drift_size.x(),
                               drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Bc3X1", "Bc3X2", "Bc3V1",
			    "Bc3V2", "Bc3U1", "Bc3U2" };
  for(G4int i=0; i<NumOfLayersBcOut/2; ++i){
    G4ThreeVector pos;
    switch (i) {
    case 0: pos.setZ(-22.*mm); break;
    case 1: pos.setZ(-18.*mm); break;
    case 2: pos.setZ( -2.*mm); break;
    case 3: pos.setZ(  2.*mm); break;
    case 4: pos.setZ( 18.*mm); break;
    case 5: pos.setZ( 22.*mm); break;
    }
    auto bc3pl_lv = new G4LogicalVolume(bc3pl_solid,
                                        m_material_map["Argon"],
                                        plane_name[i] + "LV", 0, 0, 0);
    bc3pl_lv->SetSensitiveDetector(sdc_sd);
    new G4PVPlacement(nullptr, pos, bc3pl_lv, plane_name[i] + "PV",
                      bc3_lv, false, 1001+i);
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructBC4()
{
  using CLHEP::mm;
  if(!sdc_sd){
    sdc_sd = new SDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector(sdc_sd);
  }

  const auto& bc4_pos = (gGeom.GetGlobalPosition("KURAMA") +
                         (gGeom.GetGlobalPosition("BC4-U1") +
                          gGeom.GetGlobalPosition("BC4-X2")) * 0.5);
  const auto& frame_size = gSize.GetSize("Bc4Frame") * 0.5 * mm;
  const auto& drift_size = gSize.GetSize("Bc4Drift") * 0.5 * mm;
  auto bc4_solid = new G4Box("Bc4Solid", frame_size.x(),
                             frame_size.y(), frame_size.z());
  auto bc4_lv = new G4LogicalVolume(bc4_solid, m_material_map["Argon"],
                                    "Bc4LV", 0, 0, 0);
  bc4_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(m_rotation_matrix, bc4_pos,
                    bc4_lv, "Bc4PV", m_world_lv, false, 0);
  auto bc4pl_solid = new G4Box("Bc4PlSolid", drift_size.x(),
                               drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Bc4U1", "Bc4U2", "Bc4V1",
			    "Bc4V2", "Bc4X1", "Bc4X2" };
  for(G4int i=0; i<NumOfLayersBcOut/2; ++i){
    G4ThreeVector pos;
    switch (i) {
    case 0: pos.setZ(-22.*mm); break;
    case 1: pos.setZ(-18.*mm); break;
    case 2: pos.setZ( -2.*mm); break;
    case 3: pos.setZ(  2.*mm); break;
    case 4: pos.setZ( 18.*mm); break;
    case 5: pos.setZ( 22.*mm); break;
    }
    auto bc4pl_lv = new G4LogicalVolume(bc4pl_solid,
                                        m_material_map["Argon"],
                                        plane_name[i] + "LV", 0, 0, 0);
    bc4pl_lv->SetSensitiveDetector(sdc_sd);
    new G4PVPlacement(nullptr, pos, bc4pl_lv, plane_name[i] + "PV",
                      bc4_lv, false, 2001+i);
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructK18BeamlineSpectrometer()
{
  using CLHEP::m;
  using CLHEP::mm;
  using CLHEP::deg;
  using CLHEP::kilogauss;
  G4cout << FUNC_NAME << G4endl;
  // D4 magnet
  const G4double D4Rho = gSize.Get("D4Rho")*mm;
  const G4double D4BendAngle = gSize.Get("D4BendAngle")*deg; // [deg]
  const auto& D4FieldSize = gSize.GetSize("D4Field")*mm;
  const G4double D4r1 = D4Rho - D4FieldSize.x()/2.;
  const G4double D4r2 = D4Rho + D4FieldSize.x()/2.;
  const G4double D4Width1 = 1700.*mm/2;
  // const auto& D4Coil1Size = gSize.GetSize("D4Coil1")*mm;
  // const G4double D4CoilWidth1 = 500.*mm/2;
  // const G4double D4CoilLength1 = 200.*mm/2;
  // const G4double D4CoilLength2 = 100.*mm/2;
  // const G4double D4Length = 4.468*m;
  const G4double D4B0 = 15.010222 * kilogauss;
  const G4double D4alpha = 23.5; // [deg]
  const G4double D4beta = 0; // [deg]
  // Q10 magnet
  const auto& Q10Size = gSize.GetSize("Q10")*mm;
  const G4double Q10B0 = -9.814096 * kilogauss;
  const G4double Q10a0 = 0.1*m;    // [m]
  // Q11 magnet
  const auto& Q11Size = gSize.GetSize("Q11")*mm;
  const G4double Q11B0 = 7.493605 * kilogauss;
  const G4double Q11a0 = 0.1*m;
  // Q12 magnet
  const auto& Q12Size = gSize.GetSize("Q12")*mm;
  const G4double Q12B0 = -5.87673 * kilogauss;
  const G4double Q12a0 = 0.1*m;
  // Q13 magnet
  const auto& Q13Size = gSize.GetSize("Q13")*mm;
  const G4double Q13B0 = 0. * kilogauss;
  const G4double Q13a0 = 0.1*m;      // [m]
  const G4double driftL0 = gSize.Get("K18L0")*mm; // VI-Q10
  const G4double driftL1 = gSize.Get("K18L1")*mm; // Q10-Q11
  const G4double driftL2 = gSize.Get("K18L2")*mm; // Q11-D4
  const G4double driftL3 = gSize.Get("K18L3")*mm; // D4-Q12
  const G4double driftL4 = gSize.Get("K18L4")*mm; // Q12-Q13
  const G4double driftL5 = gSize.Get("K18L5")*mm; // Q13-VO
  const G4double driftL6 = gSize.Get("K18L6")*mm; // VO-FF:1200, VO-HS
  // const G4double driftL7 = gSize.Get("K18L7")*mm; // HS-KURAMA
  const G4ThreeVector VC1Size(Q10Size.x(), Q10Size.y()/2, driftL1);
  const G4ThreeVector VC2Size(Q10Size.x(), Q10Size.y()/2, driftL2);
  const G4ThreeVector VC3Size(Q12Size.x(), Q12Size.y()/2, driftL3);
  const G4ThreeVector VC4Size(Q13Size.x(), Q13Size.y()/2, driftL4);

  G4double x,y,z;
  G4RotationMatrix rotZero;

  // VI
  x = (D4Rho*(1. - std::cos(D4BendAngle)) +
       (driftL0 + Q10Size.z() + VC1Size.z() + Q11Size.z() + VC2Size.z()) *
       std::sin(D4BendAngle));
  y = 0.*m;
  z = (- D4Rho*std::sin(D4BendAngle)
       - (driftL0 + Q10Size.z() + VC1Size.z() + Q11Size.z() + VC2Size.z()) *
       std::cos(D4BendAngle)
       - driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6);
  // BeamMan::GetInstance().SetVIPosition(G4ThreeVector(x, y, z));

  // Q10 magnet
  auto Q10Solid = new G4Box("Q10Solid", Q10Size.x()/2,
                            Q10Size.y()/2, Q10Size.z()/2);
  auto Q10LV = new G4LogicalVolume(Q10Solid,
                                   m_material_map["Vacuum"],
                                   "Q10LV");
  Q10LV->SetVisAttributes(ORANGE);
  x = (D4Rho*(1. - std::cos(D4BendAngle)) +
       (Q10Size.z()/2 + VC1Size.z() + Q11Size.z() + VC2Size.z()) *
       std::sin(D4BendAngle));
  y = 0.*m;
  z = (- D4Rho*std::sin(D4BendAngle)
       - (Q10Size.z()/2 + VC1Size.z() + Q11Size.z() + VC2Size.z()) *
       std::cos(D4BendAngle)
       - driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6);
  auto rotQ10 = new G4RotationMatrix;
  rotQ10->rotateY(D4BendAngle);
  new G4PVPlacement(rotQ10, G4ThreeVector(x, y, z),
                    Q10LV, "Q10PV", m_world_lv, false, 0);
  MagnetInfo Q10Info("Q10");
  Q10Info.type = MagnetInfo::kQuadrupole;
  Q10Info.b0 = Q10B0;
  Q10Info.a0 = Q10a0;
  Q10Info.pos = G4ThreeVector(x, y, z);
  Q10Info.size = Q10Size*0.5;
  Q10Info.ra1 = -D4BendAngle;
  m_field->AddMagnetInfo(Q10Info);
  // Q11 magnet
  auto Q11Solid = new G4Box("Q11Solid", Q11Size.x()/2,
                            Q11Size.y()/2, Q11Size.z()/2);
  auto Q11LV = new G4LogicalVolume(Q11Solid,
                                   m_material_map["Vacuum"],
                                   "Q11LV");
  Q11LV->SetVisAttributes(ORANGE);
  x = (D4Rho*(1. - std::cos(D4BendAngle)) +
       (VC2Size.z() + Q11Size.z()/2) * std::sin(D4BendAngle));
  y = 0.*m;
  z = (- D4Rho*std::sin(D4BendAngle)
       - (VC2Size.z() + Q11Size.z()/2) * std::cos(D4BendAngle)
       - driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6);
  auto rotQ11 = new G4RotationMatrix;
  rotQ11->rotateY(D4BendAngle);
  new G4PVPlacement(rotQ11, G4ThreeVector(x, y, z),
                    Q11LV, "Q11PV", m_world_lv, false, 0);
  MagnetInfo Q11Info("Q11");
  Q11Info.type = MagnetInfo::kQuadrupole;
  Q11Info.b0 = Q11B0;
  Q11Info.a0 = Q11a0;
  Q11Info.pos = G4ThreeVector(x, y, z);
  Q11Info.size = Q11Size*0.5;
  Q11Info.ra1 = -D4BendAngle;
  m_field->AddMagnetInfo(Q11Info);
  // D4 magnet
  G4cout << "   D4Rho = " << D4Rho << G4endl;
  G4cout << "   D4r1 = " << D4r1/m  << ", D4r2 = " << D4r2/m  << G4endl;
  auto D4Solid = new G4Tubs("D4Solid", D4r1, D4r2, D4FieldSize.y()/2,
                            0.*deg, D4BendAngle);
  auto D4OutSolid = new G4Tubs("D4OutSolid", D4Rho - D4Width1,
                               D4Rho + D4Width1, Q13Size.y()/2,
                               0.*deg, D4BendAngle);
  auto D4PoleSolid = new G4SubtractionSolid("D4PoleSolid", D4OutSolid,
                                            D4Solid);
  auto D4PoleLV = new G4LogicalVolume(D4PoleSolid,
                                      m_material_map["Iron"],
                                      "D4PoleLV");
  D4PoleLV->SetVisAttributes(G4Colour::Green());
  x = D4Rho;
  y = 0.*m;
  z = - driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6;
  auto D4Rot = new G4RotationMatrix;
  D4Rot->rotateX(90.*deg);
  D4Rot->rotateZ(180.*deg + D4BendAngle);
  new G4PVPlacement(D4Rot, G4ThreeVector(x, y, z),
                    D4PoleLV, "D4PolePV", m_world_lv, false, 0);
  auto D4LV = new G4LogicalVolume(D4Solid,
                                  m_material_map["Vacuum"],
                                  "D4LV");
  D4LV->SetVisAttributes(G4Colour::White());
  new G4PVPlacement(D4Rot, G4ThreeVector(x, y, z),
                    D4LV, "D4PV", m_world_lv, false, 0);
#if 0
  auto D4Coil1Solid = new G4Box("D4Coil1Solid", D4Coil1Size.x()/2,
                                D4Coil1Size.y()/2, D4Coil1Size.z()/2);
  auto D4Coil1LV = new G4LogicalVolume(D4Coil1Solid,
                                       m_material_map["Iron"],
                                       "D4Coil1LV");
  D4Coil1LV->SetVisAttributes(G4Colour::Red());
  x = 0.*m;
  y = 0.*m;
  z = (D4Coil1Size.z()/2 - driftL3 - Q12Size.z() - driftL4 - Q13Size.z()
       - driftL5 - driftL6);
  new G4PVPlacement(nullptr, G4ThreeVector(x, y, z),
                    D4Coil1LV, "D4Coil1PV", m_world_lv, false, 0);
#endif
  MagnetInfo D4Info("D4");
  D4Info.type = MagnetInfo::kDipole;
  D4Info.b0 = D4B0;
  D4Info.rho = D4Rho;
  D4Info.alpha = D4alpha;
  D4Info.beta = D4beta;
  D4Info.pos = G4ThreeVector(x, y, z);
  D4Info.size = D4FieldSize*0.5;
  D4Info.ra1 = -D4BendAngle;
  D4Info.bend = D4BendAngle/deg;
  m_field->AddMagnetInfo(D4Info);
  // Q12 magnet
  auto Q12Solid = new G4Box("Q12Solid", Q12Size.x()/2,
                            Q12Size.y()/2, Q12Size.z()/2);
  auto Q12LV = new G4LogicalVolume(Q12Solid,
                                   m_material_map["Vacuum"],
                                   "Q12LV");
  Q12LV->SetVisAttributes(ORANGE);
  x = 0.*m;
  y = 0.*m;
  z = - Q12Size.z()/2 - driftL4 - Q13Size.z() - driftL5 - driftL6;
  new G4PVPlacement(nullptr, G4ThreeVector(x, y, z),
                    Q12LV, "Q12PV", m_world_lv, false, 0);
  MagnetInfo Q12Info("Q12");
  Q12Info.type = MagnetInfo::kQuadrupole;
  Q12Info.b0 = Q12B0;
  Q12Info.a0 = Q12a0;
  Q12Info.pos = G4ThreeVector(x, y, z);
  Q12Info.size = Q12Size*0.5;
  Q12Info.ra1 = -D4BendAngle;
  m_field->AddMagnetInfo(Q12Info);
  // Q13 magnet
  auto Q13Solid = new G4Box("Q13Solid", Q13Size.x()/2,
                            Q13Size.y()/2, Q13Size.z()/2);
  auto Q13LV = new G4LogicalVolume(Q13Solid,
                                   m_material_map["Vacuum"],
                                   "Q13LV");
  Q13LV->SetVisAttributes(ORANGE);
  x = 0.*m;
  y = 0.*m;
  z =  - Q13Size.z()/2 - driftL5 - driftL6;
  new G4PVPlacement(nullptr, G4ThreeVector(x, y, z),
                    Q13LV, "Q13PV", m_world_lv, false, 0);
  MagnetInfo Q13Info("Q13");
  Q13Info.type = MagnetInfo::kQuadrupole;
  Q13Info.b0 = Q13B0;
  Q13Info.a0 = Q13a0;
  Q13Info.pos = G4ThreeVector(x, y, z);
  Q13Info.size = Q13Size*0.5;
  Q13Info.ra1 = -D4BendAngle;
  m_field->AddMagnetInfo(Q13Info);
  // Vacuum Chamber 1
  auto VC1Solid = new G4Box("VC1Solid", VC1Size.x()/2,
                            VC1Size.y()/2, VC1Size.z()/2);
  auto VC1LV = new G4LogicalVolume(VC1Solid,
                                   m_material_map["Vacuum"],
                                   "VC1LV");
  x = (D4Rho*(1. - std::cos(D4BendAngle)) +
       (VC1Size.z()/2 + Q11Size.z() + VC2Size.z()) * std::sin(D4BendAngle));
  y = 0.*m;
  z = (- D4Rho*std::sin(D4BendAngle)
       - (VC1Size.z()/2 + Q11Size.z() + VC2Size.z()) * std::cos(D4BendAngle)
       - driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6);
  auto rotVC1 = new G4RotationMatrix;
  rotVC1->rotateY(D4BendAngle);
  new G4PVPlacement(rotVC1, G4ThreeVector(x, y, z),
                    VC1LV, "VC1PV", m_world_lv, false, 0);
  // Vacuum Chamber 2
  auto VC2Solid = new G4Box("slidVC2", VC2Size.x()/2,
                            VC2Size.y()/2, VC2Size.z()/2);
  auto VC2LV = new G4LogicalVolume(VC2Solid,
                                   m_material_map["Vacuum"],
                                   "VC2LV");
  x = (D4Rho*(1. - std::cos(D4BendAngle)) +
       VC2Size.z()/2 * std::sin(D4BendAngle));
  y = 0.*m;
  z = (- D4Rho*std::sin(D4BendAngle)
       - VC2Size.z()/2 * std::cos(D4BendAngle)
       - driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6);
  auto rotVC2 = new G4RotationMatrix;
  rotVC2->rotateY(D4BendAngle);
  new G4PVPlacement(rotVC2, G4ThreeVector(x, y, z),
                    VC2LV, "VC2PV", m_world_lv, false, 0);
  // Vacuum Chamber 3
  auto VC3Solid = new G4Box("VC3Solid", VC3Size.x()/2,
                            VC3Size.y()/2, VC3Size.z()/2);
  auto VC3LV = new G4LogicalVolume(VC3Solid,
                                   m_material_map["Vacuum"],
                                   "VC3LV");
  x = 0.*m;
  y = 0.*m;
  z = -driftL3/2 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6;
  new G4PVPlacement(nullptr, G4ThreeVector(x, y, z),
                    VC3LV, "VC3PV", m_world_lv, false, 0);
  // Vacuum Chamber 4
  auto VC4Solid = new G4Box("VC4Solid", VC4Size.x()/2, VC4Size.y()/2, VC4Size.z()/2);
  auto VC4LV = new G4LogicalVolume(VC4Solid,
                                   m_material_map["Vacuum"],
                                   "VC4LV");
  x = 0.*m;
  y = 0.*m;
  z = - driftL4/2 - Q13Size.z() - driftL5 - driftL6;
  new G4PVPlacement(nullptr, G4ThreeVector(x, y, z),
                    VC4LV, "VC4PV", m_world_lv, false, 0);
  m_field->SetStatusK18Field(true);
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructKuramaMagnet()
{
  using CLHEP::m;
  using CLHEP::mm;
  using CLHEP::deg;
  auto vp_sd = new VPSD("/VP");
  G4SDManager::GetSDMpointer()->AddNewDetector(vp_sd);
  const auto field_size = gSize.GetSize("KuramaField") * 0.5 * mm;
  const auto coil_color = G4Colour::Red();
  // size
  const G4ThreeVector coil1_size(900.*mm/2, 193.*mm/2, 280.*mm/2);
  const G4ThreeVector coil2_size(193.*mm/2, 117.*mm/2, 280.*mm/2);
  const G4ThreeVector coil3_size(193.*mm/2, 137.*mm/2, 280.*mm/2);
  const G4ThreeVector coil4_size(193.*mm/2, 280.*mm/2, 740.*mm/2);
  const G4ThreeVector coil5_size(900.*mm/2, 214.*mm/2, 280.*mm/2);
  const G4ThreeVector yoke_inner_size(field_size.x() + 193.*mm,
                                      field_size.y(), field_size.z());
  const G4ThreeVector yoke_outer_size(1800.*mm/2, 2200.*mm/2,
                                      field_size.z());
  const G4ThreeVector yoke_ud_size(2200.*mm/2, 370.*mm/2, field_size.z());
  const G4ThreeVector yoke_lr_size(200.*mm/2, field_size.y(),
                                   field_size.z());
  const G4ThreeVector uguard_inner_size(600.*mm/2, 300.*mm/2, 100.*mm/2);
  const G4ThreeVector uguard_inner2_size(50.*mm/2, 50.*mm/2, 100.*mm/2);
  const G4ThreeVector uguard_outer_size(1600.*mm/2, 2200.*mm/2, 100.*mm/2);
  const G4ThreeVector dguard_inner_size(1100.*mm/2, 1100.*mm/2, 100.*mm/2);
  const G4ThreeVector dguard_outer_size(1600.*mm/2, 2200.*mm/2, 100.*mm/2);
  // position
  const G4ThreeVector field_pos = gGeom.GetGlobalPosition("KURAMA");
  const G4ThreeVector yoke_pos(field_pos);
  const G4ThreeVector uguard_pos(field_pos.x(), field_pos.y(),
                                 field_pos.z() - (820. - 50.)*mm);
  const G4ThreeVector dguard_pos(field_pos.x(), field_pos.y(),
                                 field_pos.z() + (820. - 50.)*mm);
  const G4ThreeVector coil1_pos =
    { field_pos.x(),
      field_size.y() + yoke_ud_size.y()*2. - (coil1_size.y() + 20.*mm),
      field_pos.z() - field_size.z() - (coil1_size.z() + 20.*mm) };
  const G4ThreeVector coil4l_pos =
    { field_pos.x() + field_size.x() + coil4_size.x(),
      field_size.y()/2,
      field_pos.z() };
  const G4ThreeVector coil4r_pos =
    { field_pos.x() - field_size.x() - coil4_size.x(),
      field_size.y()/2,
      field_pos.z() };
  const G4ThreeVector coil2l_pos =
    { field_pos.x() + field_size.x() + coil4_size.x(),
      (coil4l_pos.y() + coil1_pos.y() + coil4_size.y() - coil1_size.y())/2,
      field_pos.z() - field_size.z() - coil1_size.z() + 20.*mm };
  const G4ThreeVector coil2r_pos =
    { field_pos.x() - field_size.x() - coil4_size.x(),
      (coil4r_pos.y() + coil1_pos.y() + coil4_size.y() - coil1_size.y())/2,
      field_pos.z() - field_size.z() - coil1_size.z() + 20.*mm };
  const G4ThreeVector coil5_pos =
    { field_pos.x(),
      field_size.y() + yoke_ud_size.y()*2. - coil5_size.y(),
      field_pos.z() + field_size.z() + coil5_size.z() + 21.*mm };
  const G4ThreeVector coil3l_pos =
    { field_pos.x() + field_size.x() + coil4_size.x(),
      (coil4l_pos.y() + coil5_pos.y() + coil4_size.y() - coil5_size.y())/2,
      field_pos.z() + field_size.z() + coil5_size.z() + 21.*mm };
  const G4ThreeVector coil3r_pos =
    { field_pos.x() - field_size.x() - coil4_size.x(),
      (coil4r_pos.y() + coil5_pos.y() + coil4_size.y() - coil5_size.y())/2,
      field_pos.z() + field_size.z() + coil5_size.z() + 21.*mm };
  const G4ThreeVector yoke_u_pos(field_pos.x(),
                                 field_size.y() + yoke_ud_size.y(),
                                 field_pos.z());
  const G4ThreeVector yoke_d_pos(field_pos.x(),
                                 - field_size.y() - yoke_ud_size.y(),
                                 field_pos.z());
  const G4ThreeVector yoke_l_pos =
    { field_pos.x() + field_size.x() + yoke_lr_size.x() + 200.*mm,
      0.*mm,
      field_pos.z() };
  const G4ThreeVector yoke_r_pos =
    { field_pos.x() - field_size.x() + yoke_lr_size.x() - 200.*mm,
      0.*mm,
      field_pos.z() };
  // Construct KURAMA Magnet
  // auto kurama_solid = new G4Box("KuramaSolid", 4.*m/2, 3.*m/2, 4.*m/2);
  // auto kurama_lv = new G4LogicalVolume(kurama_solid, m_material_map["Air"],
  // 					"KuramaLV");
  // kurama_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // auto kurama_pv = new G4PVPlacement(m_rotation_angle, G4ThreeVector(), "KuramaPV",
  // 				      kurama_lv, m_world_pv, false, 0);

  // Coil1
  auto coil1_solid = new G4Box("Coil1Solid", coil1_size.x(),
                               coil1_size.y(), coil1_size.z());
  auto coil1_lv = new G4LogicalVolume(coil1_solid, m_material_map["Copper"],
                                      "Coil1LV");
  coil1_lv->SetVisAttributes(coil_color);
  auto coil1u_pos = coil1_pos;
  coil1u_pos.rotateY(m_rotation_angle);
  new G4PVPlacement(m_rotation_matrix, coil1u_pos, coil1_lv,
                    "Coil1UPV", m_world_lv, false, 0);
  auto coil1d_pos = G4ThreeVector(coil1_pos.x(),
                                  -coil1_pos.y(),
                                  coil1_pos.z());
  coil1d_pos.rotateY(m_rotation_angle);
  new G4PVPlacement(m_rotation_matrix,	coil1d_pos, coil1_lv,
                    "Coil1DPV", m_world_lv, false, 0);

  // Coil4
  auto coil4_solid = new G4Box("Coil4Solid", coil4_size.x(),
                               coil4_size.y(), coil4_size.z());
  auto coil4_lv = new G4LogicalVolume(coil4_solid, m_material_map["Copper"],
                                      "Coil4LV");
  coil4_lv->SetVisAttributes(coil_color);
  auto coil4ur_pos = coil4l_pos;
  coil4ur_pos.rotateY(m_rotation_angle);
  new G4PVPlacement(m_rotation_matrix, coil4ur_pos, coil4_lv,
                    "Coil4URPV", m_world_lv, false, 0);
  auto coil4ul_pos = coil4r_pos;
  coil4ul_pos.rotateY(m_rotation_angle);
  new G4PVPlacement(m_rotation_matrix, coil4ul_pos, coil4_lv,
                    "Coil4ULPV", m_world_lv, false, 0);
  auto coil4dr_pos = G4ThreeVector(coil4l_pos.x(),
                                   -coil4l_pos.y(),
                                   coil4l_pos.z()).rotateY(m_rotation_angle);
  new G4PVPlacement(m_rotation_matrix, coil4dr_pos, coil4_lv,
                    "Coil4DRPV", m_world_lv, false, 0);
  auto coil4dl_pos = G4ThreeVector(coil4r_pos.x(),
                                   -coil4r_pos.y(),
                                   coil4r_pos.z()).rotateY(m_rotation_angle);
  new G4PVPlacement(m_rotation_matrix, coil4ul_pos, coil4_lv,
                    "Coil4DLPV", m_world_lv, false, 0);

  // Coil5
  auto coil5_solid = new G4Box("Coil5Solid", coil5_size.x(),
                               coil5_size.y(), coil5_size.z());
  auto coil5_lv = new G4LogicalVolume(coil5_solid, m_material_map["Copper"],
                                      "Coil5LV");
  coil5_lv->SetVisAttributes(coil_color);
  auto coil5u_pos = G4ThreeVector(coil5_pos.x(),
                                  coil5_pos.y(),
                                  coil5_pos.z()).rotateY(m_rotation_angle);
  new G4PVPlacement(m_rotation_matrix,	coil5u_pos, coil5_lv,
                    "Coil5UPV", m_world_lv, false, 0);
  auto coil5d_pos = G4ThreeVector(coil5_pos.x(),
                                  -coil5_pos.y(),
                                  coil5_pos.z()).rotateY(m_rotation_angle);
  new G4PVPlacement(m_rotation_matrix, coil5d_pos, coil5_lv,
                    "Coil5DPV", m_world_lv, false, 0);

  // Coil6
  G4double size_COIL6[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL6[0] = 50.0*mm;
  size_COIL6[1] = coil1_size.y()*2+size_COIL6[0];
  //  size_COIL6[1] = 330.*mm;
  size_COIL6[2] = 280.*mm/2;
  size_COIL6[3] = 90.*deg;

  G4double pos_COIL6LU[3];
  G4double pos_COIL6RU[3];
  G4double pos_COIL6LD[3];
  G4double pos_COIL6RD[3];
  //LU
  pos_COIL6LU[ThreeVector::X] = field_pos.x() +field_size.x()-size_COIL6[0];
  pos_COIL6LU[ThreeVector::Y] = coil1_pos.y()  -(size_COIL6[0]+coil1_size.y());
  pos_COIL6LU[ThreeVector::Z] = field_pos.z() - field_size.z()-(size_COIL6[ThreeVector::Z]+21.*mm);
  //RU
  pos_COIL6RU[ThreeVector::X] = field_pos.x() -field_size.x()+size_COIL6[0];
  pos_COIL6RU[ThreeVector::Y] = coil1_pos.y()  -(size_COIL6[0]+coil1_size.y());
  pos_COIL6RU[ThreeVector::Z] = field_pos.z() - field_size.z()-(size_COIL6[ThreeVector::Z]+21.*mm);
  //LD
  pos_COIL6LD[ThreeVector::X] = field_pos.x()  +field_size.x()-size_COIL6[0];
  pos_COIL6LD[ThreeVector::Y] = -coil1_pos.y() +(size_COIL6[0]+coil1_size.y());
  pos_COIL6LD[ThreeVector::Z] = field_pos.z()  - field_size.z()-(size_COIL6[ThreeVector::Z]+21.*mm);
  //RD
  pos_COIL6RD[ThreeVector::X] = field_pos.x()  -field_size.x()+size_COIL6[0];
  pos_COIL6RD[ThreeVector::Y] = -coil1_pos.y() +(size_COIL6[0]+coil1_size.y());
  pos_COIL6RD[ThreeVector::Z] = field_pos.z()  - field_size.z()-(size_COIL6[ThreeVector::Z]+21.*mm);


  G4Tubs* Coil6_tub = new G4Tubs("Coil6_tubs",
				 size_COIL6[0],size_COIL6[1],size_COIL6[2],0.,size_COIL6[3]);
  G4LogicalVolume*  Coil6_log = new G4LogicalVolume(Coil6_tub, m_material_map["Copper"], "Coil6_log",0,0,0);
  Coil6_log->SetVisAttributes(coil_color);
  G4RotationMatrix *rotcoil6lu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6ru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6ld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6rd = new G4RotationMatrix();
  rotcoil6lu->rotateY(- m_rotation_angle);
  rotcoil6ru->rotateY(- m_rotation_angle);
  rotcoil6ld->rotateY(- m_rotation_angle);
  rotcoil6rd->rotateY(- m_rotation_angle);

  rotcoil6lu->rotateZ(0.*deg);
  new G4PVPlacement(rotcoil6lu,
                    G4ThreeVector(pos_COIL6LU[ThreeVector::X],
                                  pos_COIL6LU[ThreeVector::Y],
                                  pos_COIL6LU[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil6_log,
                    "Coil6LU_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil6ru->rotateZ(-90.*deg);
  new G4PVPlacement(rotcoil6ru,
                    G4ThreeVector(pos_COIL6RU[ThreeVector::X],
                                  pos_COIL6RU[ThreeVector::Y],
                                  pos_COIL6RU[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil6_log,
                    "Coil6RU_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil6ld->rotateZ(-180.*deg);
  new G4PVPlacement(rotcoil6ld,
                    G4ThreeVector(pos_COIL6RD[ThreeVector::X],
                                  pos_COIL6RD[ThreeVector::Y],
                                  pos_COIL6RD[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil6_log,
                    "Coil6RD_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil6rd->rotateZ(-270.*deg);
  new G4PVPlacement(rotcoil6rd,
                    G4ThreeVector(pos_COIL6LD[ThreeVector::X],
                                  pos_COIL6LD[ThreeVector::Y],
                                  pos_COIL6LD[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil6_log,
                    "Coil6LD_phys",
                    m_world_lv,
                    false,
                    0);
  //  m_rotation_matrix->rotateZ(90.*deg);
  //  m_rotation_matrix->rotateZ(90.*deg);

  //////////////coil8RLUD
  G4double size_COIL8[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL8[0] = 50.0*mm;
  size_COIL8[1] = coil5_size.y()*2.+size_COIL8[0];
  size_COIL8[2] = 280.*mm/2;
  size_COIL8[3] = 90.*deg;

  G4double pos_COIL8LU[3];
  G4double pos_COIL8RU[3];
  G4double pos_COIL8LD[3];
  G4double pos_COIL8RD[3];
  //LU
  pos_COIL8LU[ThreeVector::X] = field_pos.x() +field_size.x()-size_COIL8[0];
  pos_COIL8LU[ThreeVector::Y] = coil5_pos.y()  -(size_COIL8[0]+coil5_size.y());
  pos_COIL8LU[ThreeVector::Z] = field_pos.z() + field_size.z()+(size_COIL8[ThreeVector::Z]+21.*mm);
  //RU
  pos_COIL8RU[ThreeVector::X] = field_pos.x() -field_size.x()+size_COIL8[0];
  pos_COIL8RU[ThreeVector::Y] = coil5_pos.y()  -(size_COIL8[0]+coil5_size.y());
  pos_COIL8RU[ThreeVector::Z] = field_pos.z() + field_size.z()+(size_COIL8[ThreeVector::Z]+21.*mm);
  //LD
  pos_COIL8LD[ThreeVector::X] = field_pos.x()  +field_size.x()-size_COIL8[0];
  pos_COIL8LD[ThreeVector::Y] = -coil5_pos.y() +(size_COIL8[0]+coil5_size.y());
  pos_COIL8LD[ThreeVector::Z] = field_pos.z()  + field_size.z()+(size_COIL8[ThreeVector::Z]+21.*mm);
  //RD
  pos_COIL8RD[ThreeVector::X] = field_pos.x()  -field_size.x()+size_COIL8[0];
  pos_COIL8RD[ThreeVector::Y] = -coil5_pos.y() +(size_COIL8[0]+coil5_size.y());
  pos_COIL8RD[ThreeVector::Z] = field_pos.z()  + field_size.z()+(size_COIL8[ThreeVector::Z]+21.*mm);

  G4Tubs* Coil8_tub = new G4Tubs("Coil8_tubs",
				 size_COIL8[0],size_COIL8[1],size_COIL8[2],0.,size_COIL8[3]);
  G4LogicalVolume*  Coil8_log = new G4LogicalVolume(Coil8_tub, m_material_map["Copper"], "Coil8_log",0,0,0);
  Coil8_log->SetVisAttributes(coil_color);
  G4RotationMatrix *rotcoil8lu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8ru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8ld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8rd = new G4RotationMatrix();
  rotcoil8lu->rotateY(- m_rotation_angle);
  rotcoil8ru->rotateY(- m_rotation_angle);
  rotcoil8ld->rotateY(- m_rotation_angle);
  rotcoil8rd->rotateY(- m_rotation_angle);

  rotcoil8lu->rotateZ(0.*deg);
  new G4PVPlacement(rotcoil8lu,
                    G4ThreeVector(pos_COIL8LU[ThreeVector::X],
                                  pos_COIL8LU[ThreeVector::Y],
                                  pos_COIL8LU[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil8_log,
                    "Coil8LU_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil8ru->rotateZ(-90.*deg);
  new G4PVPlacement(rotcoil8ru,
                    G4ThreeVector(pos_COIL8RU[ThreeVector::X],
                                  pos_COIL8RU[ThreeVector::Y],
                                  pos_COIL8RU[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil8_log,
                    "Coil8RU_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil8ld->rotateZ(-180.*deg);
  new G4PVPlacement(rotcoil8ld,
                    G4ThreeVector(pos_COIL8RD[ThreeVector::X],
                                  pos_COIL8RD[ThreeVector::Y],
                                  pos_COIL8RD[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil8_log,
                    "Coil8RD_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil8rd->rotateZ(-270.*deg);
  new G4PVPlacement(rotcoil8rd,
                    G4ThreeVector(pos_COIL8LD[ThreeVector::X],
                                  pos_COIL8LD[ThreeVector::Y],
                                  pos_COIL8LD[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil8_log,
                    "Coil8LD_phys",
                    m_world_lv,
                    false,
                    0);
  //  m_rotation_matrix->rotateZ(90.*deg);
  //  m_rotation_matrix->rotateZ(90.*deg);

  //////////////coil7RLUD
  G4double size_COIL7[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL7[0] = 50.0*mm;
  size_COIL7[1] = coil4_size.y()*2.+size_COIL7[0];
  size_COIL7[2] = coil4_size.x();
  size_COIL7[3] = 90.*deg;

  G4double pos_COIL7ULU[3];
  G4double pos_COIL7URU[3];
  G4double pos_COIL7ULD[3];
  G4double pos_COIL7URD[3];
  G4double pos_COIL7DLU[3];
  G4double pos_COIL7DRU[3];
  G4double pos_COIL7DLD[3];
  G4double pos_COIL7DRD[3];
  //ULU
  // pos_COIL7ULU[ThreeVector::X] = field_pos.x() +field_size.x()+size_COIL7[0];
  pos_COIL7ULU[ThreeVector::X] = coil4l_pos.x();
  pos_COIL7ULU[ThreeVector::Y] = coil4l_pos.y() +(size_COIL7[0]+coil4_size.y());
  pos_COIL7ULU[ThreeVector::Z] = field_pos.z() - coil4_size.z();
  //URU
  pos_COIL7URU[ThreeVector::X] = coil4r_pos.x();
  pos_COIL7URU[ThreeVector::Y] = coil4r_pos.y() +(size_COIL7[0]+coil4_size.y());
  pos_COIL7URU[ThreeVector::Z] = field_pos.z() - coil4_size.z();
  //ULD
  pos_COIL7ULD[ThreeVector::X] = coil4l_pos.x();
  pos_COIL7ULD[ThreeVector::Y] = -coil4l_pos.y() -(size_COIL7[0]+coil4_size.y());
  pos_COIL7ULD[ThreeVector::Z] = field_pos.z()  - coil4_size.z();
  //URD
  pos_COIL7URD[ThreeVector::X] = coil4r_pos.x();
  pos_COIL7URD[ThreeVector::Y] = -coil4r_pos.y() -(size_COIL7[0]+coil4_size.y());
  pos_COIL7URD[ThreeVector::Z] = field_pos.z() - coil4_size.z();


  //DLU
  //  pos_COIL7ULU[ThreeVector::X] = field_pos.x() +field_size.x()+size_COIL7[0];
  pos_COIL7DLU[ThreeVector::X] = coil4l_pos.x();
  pos_COIL7DLU[ThreeVector::Y] = coil4l_pos.y() +(size_COIL7[0]+coil4_size.y());
  pos_COIL7DLU[ThreeVector::Z] = field_pos.z() + coil4_size.z();
  //DRU
  pos_COIL7DRU[ThreeVector::X] = coil4r_pos.x();
  pos_COIL7DRU[ThreeVector::Y] = coil4r_pos.y() +(size_COIL7[0]+coil4_size.y());
  pos_COIL7DRU[ThreeVector::Z] = field_pos.z() + coil4_size.z();
  //DLD
  pos_COIL7DLD[ThreeVector::X] = coil4l_pos.x();
  pos_COIL7DLD[ThreeVector::Y] = -coil4l_pos.y() -(size_COIL7[0]+coil4_size.y());
  pos_COIL7DLD[ThreeVector::Z] = field_pos.z()  + coil4_size.z();
  //DRD
  pos_COIL7DRD[ThreeVector::X] = coil4r_pos.x();
  pos_COIL7DRD[ThreeVector::Y] = -coil4r_pos.y() -(size_COIL7[0]+coil4_size.y());
  pos_COIL7DRD[ThreeVector::Z] = field_pos.z() + coil4_size.z();


  G4Tubs* Coil7_tub = new G4Tubs("Coil7_tubs",
				 size_COIL7[0],size_COIL7[1],size_COIL7[2],0.,size_COIL7[3]);
  G4LogicalVolume*  Coil7_log = new G4LogicalVolume(Coil7_tub, m_material_map["Copper"], "Coil7_log",0,0,0);
  Coil7_log->SetVisAttributes(coil_color);
  G4RotationMatrix *rotcoil7ulu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7uru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7uld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7urd = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dlu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7drd = new G4RotationMatrix();
  rotcoil7ulu->rotateY(- m_rotation_angle);
  rotcoil7uru->rotateY(- m_rotation_angle);
  rotcoil7uld->rotateY(- m_rotation_angle);
  rotcoil7urd->rotateY(- m_rotation_angle);
  rotcoil7dlu->rotateY(- m_rotation_angle);
  rotcoil7dru->rotateY(- m_rotation_angle);
  rotcoil7dld->rotateY(- m_rotation_angle);
  rotcoil7drd->rotateY(- m_rotation_angle);

  rotcoil7ulu->rotateZ(0.*deg);
  rotcoil7ulu->rotateX(180.*deg);
  rotcoil7ulu->rotateY(90.*deg);
  new G4PVPlacement(rotcoil7ulu,
                    G4ThreeVector(pos_COIL7ULU[ThreeVector::X],
                                  pos_COIL7ULU[ThreeVector::Y],
                                  pos_COIL7ULU[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil7_log,
                    "Coil7ULU_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil7uru->rotateZ(0.*deg);
  rotcoil7uru->rotateX(180.*deg);
  rotcoil7uru->rotateY(90.*deg);
  new G4PVPlacement(rotcoil7uru,
                    G4ThreeVector(pos_COIL7URU[ThreeVector::X],
                                  pos_COIL7URU[ThreeVector::Y],
                                  pos_COIL7URU[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil7_log,
                    "Coil7URU_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil7urd->rotateZ(0.*deg);
  rotcoil7urd->rotateX(90.*deg);
  rotcoil7urd->rotateY(90.*deg);
  new G4PVPlacement(rotcoil7urd,
                    G4ThreeVector(pos_COIL7URD[ThreeVector::X],
                                  pos_COIL7URD[ThreeVector::Y],
                                  pos_COIL7URD[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil7_log,
                    "Coil7URD_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil7uld->rotateZ(0.*deg);
  rotcoil7uld->rotateX(90.*deg);
  rotcoil7uld->rotateY(90.*deg);
  new G4PVPlacement(rotcoil7uld,
                    G4ThreeVector(pos_COIL7ULD[ThreeVector::X],
                                  pos_COIL7ULD[ThreeVector::Y],
                                  pos_COIL7ULD[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil7_log,
                    "Coil7ULD_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil7dlu->rotateZ(0.*deg);
  rotcoil7dlu->rotateX(-90.*deg);
  rotcoil7dlu->rotateY(90.*deg);
  new G4PVPlacement(rotcoil7dlu,
                    G4ThreeVector(pos_COIL7DLU[ThreeVector::X],
                                  pos_COIL7DLU[ThreeVector::Y],
                                  pos_COIL7DLU[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil7_log,
                    "Coil7DLU_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil7dru->rotateZ(0.*deg);
  rotcoil7dru->rotateX(-90.*deg);
  rotcoil7dru->rotateY(90.*deg);
  new G4PVPlacement(rotcoil7dru,
                    G4ThreeVector(pos_COIL7DRU[ThreeVector::X],
                                  pos_COIL7DRU[ThreeVector::Y],
                                  pos_COIL7DRU[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil7_log,
                    "Coil7DRU_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil7drd->rotateZ(0.*deg);
  rotcoil7drd->rotateX(0.*deg);
  rotcoil7drd->rotateY(90.*deg);
  new G4PVPlacement(rotcoil7drd,
                    G4ThreeVector(pos_COIL7DRD[ThreeVector::X],
                                  pos_COIL7DRD[ThreeVector::Y],
                                  pos_COIL7DRD[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil7_log,
                    "Coil7DRD_phys",
                    m_world_lv,
                    false,
                    0);
  rotcoil7dld->rotateZ(0.*deg);
  rotcoil7dld->rotateX(0.*deg);
  rotcoil7dld->rotateY(90.*deg);
  new G4PVPlacement(rotcoil7dld,
                    G4ThreeVector(pos_COIL7DLD[ThreeVector::X],
                                  pos_COIL7DLD[ThreeVector::Y],
                                  pos_COIL7DLD[ThreeVector::Z]).rotateY(m_rotation_angle),
                    Coil7_log,
                    "Coil7DLD_phys",
                    m_world_lv,
                    false,
                    0);

  ///coil2
  //////////////coil2
  G4Box* Coil2_box = new G4Box("Coil2_box",
                               coil2_size.x(),
                               coil2_size.y(),
                               coil2_size.z());
  G4LogicalVolume*  Coil2_log = new G4LogicalVolume(Coil2_box, m_material_map["Copper"], "Coil2_log",0,0,0);
  Coil2_log->SetVisAttributes(coil_color);
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(coil2l_pos.x(),
                                  coil2l_pos.y(),
                                  coil2l_pos.z()).rotateY(m_rotation_angle),
                    Coil2_log,
                    "Coil2UL_phys",
                    m_world_lv,
                    false,
                    0);
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(coil2r_pos.x(),
                                  coil2r_pos.y(),
                                  coil2r_pos.z()).rotateY(m_rotation_angle),
                    Coil2_log,
                    "Coil2UR_phys",
                    m_world_lv,
                    false,
                    0);
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(coil2l_pos.x(),
                                  -coil2l_pos.y(),
                                  coil2l_pos.z()).rotateY(m_rotation_angle),
                    Coil2_log,
                    "Coil2DL_phys",
                    m_world_lv,
                    false,
                    0);
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(coil2r_pos.x(),
                                  -coil2r_pos.y(),
                                  coil2r_pos.z()).rotateY(m_rotation_angle),
                    Coil2_log,
                    "Coil2DR_phys",
                    m_world_lv,
                    false,
                    0);

  ///coil3
  //////////////coil3
  G4Box* Coil3_box = new G4Box("Coil3_box",
			       coil3_size.x(), coil3_size.y(), coil3_size.z());
  G4LogicalVolume*  Coil3_log = new G4LogicalVolume(Coil3_box, m_material_map["Copper"], "Coil3_log",0,0,0);
  Coil3_log->SetVisAttributes(coil_color);
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(coil3l_pos.x(),
                                  coil3l_pos.y(),
                                  coil3l_pos.z()).rotateY(m_rotation_angle),
                    Coil3_log,
                    "Coil3UL_phys",
                    m_world_lv,
                    false,
                    0);
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(coil3r_pos.x(),
                                  coil3r_pos.y(),
                                  coil3r_pos.z()).rotateY(m_rotation_angle),
                    Coil3_log,
                    "Coil3UR_phys",
                    m_world_lv,
                    false,
                    0);
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(coil3l_pos.x(),
                                  -coil3l_pos.y(),
                                  coil3l_pos.z()).rotateY(m_rotation_angle),
                    Coil3_log,
                    "Coil3DL_phys",
                    m_world_lv,
                    false,
                    0);
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(coil3r_pos.x(),
                                  -coil3r_pos.y(),
                                  coil3r_pos.z()).rotateY(m_rotation_angle),
                    Coil3_log,
                    "Coil3DR_phys",
                    m_world_lv,
                    false,
                    0);
  // Upstream End Guard
  auto uguard_inner_solid = new G4Box("UGuardInnerSolid",
                                      uguard_inner_size.x(),
                                      uguard_inner_size.y(),
                                      uguard_inner_size.z());
  auto uguard_inner2_solid = new G4Box("UGuardInner2Solid",
                                       uguard_inner2_size.x(),
                                       uguard_inner2_size.y(),
                                       uguard_inner2_size.z());
  auto uguard_outer_solid = new G4Box("UGuardOuterSolid",
                                      uguard_outer_size.x(),
                                      uguard_outer_size.y(),
                                      uguard_outer_size.z());
  auto sub_solid = new G4SubtractionSolid("SubSolid",
                                          uguard_outer_solid,
                                          uguard_inner_solid);
  G4ThreeVector uhole(-uguard_inner_size.x()-uguard_inner2_size.x(), 0., 0.);
  auto uguard_solid = new G4SubtractionSolid("UGuardSolid",
                                             sub_solid, uguard_inner2_solid,
                                             nullptr, uhole);
  auto uguard_lv = new G4LogicalVolume(uguard_solid, m_material_map["Iron"],
                                       "UGuardLV");
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(uguard_pos).rotateY(m_rotation_angle),
                    uguard_lv, "UGuardPV", m_world_lv, false, 0);
  // Yoke
  auto yoke_inner_solid = new G4Box("YokeInnerSolid",
                                    yoke_inner_size.x(),
                                    yoke_inner_size.y(),
                                    yoke_inner_size.z());
  auto yoke_outer_solid = new G4Box("YokeOuterSolid",
                                    yoke_outer_size.x(),
                                    yoke_outer_size.y(),
                                    yoke_outer_size.z());
  auto yoke_solid = new G4SubtractionSolid("YokeSolid",
                                           yoke_outer_solid,
                                           yoke_inner_solid);
  auto yoke_lv = new G4LogicalVolume(yoke_solid, m_material_map["Iron"],
                                     "YokeLV");
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(yoke_pos).rotateY(m_rotation_angle),
                    yoke_lv, "YokePV", m_world_lv, false, 0);
  // Downstream End Guard
  auto dguard_inner_solid = new G4Box("DGuardInnerSolid",
                                      dguard_inner_size.x(),
                                      dguard_inner_size.y(),
                                      dguard_inner_size.z());
  auto dguard_outer_solid = new G4Box("DGuardOuterSolid",
                                      dguard_outer_size.x(),
                                      dguard_outer_size.y(),
                                      dguard_outer_size.z());
  auto dguard_solid = new G4SubtractionSolid("DGuardSolid",
                                             dguard_outer_solid,
                                             dguard_inner_solid);
  auto dguard_lv = new G4LogicalVolume(dguard_solid, m_material_map["Iron"],
                                       "DGuardLV");
  new G4PVPlacement(m_rotation_matrix,
                    G4ThreeVector(dguard_pos).rotateY(m_rotation_angle),
                    dguard_lv, "DGuardPV", m_world_lv, false, 0);
  // Virtual Plane
  auto vp_solid = new G4Box("VPSolid", 2.*m/2, 2.*m/2, 0.001*mm/2);
  auto vp_lv = new G4LogicalVolume(vp_solid, m_material_map["Air"], "VPLV");
  vp_lv->SetSensitiveDetector(vp_sd);
  vp_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  for(G4int i=0; i<NumOfSegVP; ++i){
    auto vp_pos = (gGeom.GetGlobalPosition("KURAMA") +
                   gGeom.GetGlobalPosition(Form("VP%d", i+1)));
    new G4PVPlacement(m_rotation_matrix, vp_pos, vp_lv,
                      Form("VP%dPV", i+1), m_world_lv, false, i);
  }
  m_field->SetStatusKuramaField(true);
  m_field->SetKuramaFieldMap(gConf.Get<G4String>("KURAMAFLDMAP"));
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructLAC()
{
  using CLHEP::mm;
  using CLHEP::deg;
  auto lac_sd = new LACSD("/LAC");
  lac_sd->SetRefractiveIndex(1.05);
  G4SDManager::GetSDMpointer()->AddNewDetector(lac_sd);
  const auto& ra2 = gGeom.GetRotAngle2("LAC") * deg;
  const auto& frame_size = gSize.GetSize("LacFrame") * 0.5 * mm;
  const auto& radiator_size = gSize.GetSize("LacRadiator") * 0.5 * mm;
  // Mother
  auto mother_solid = new G4Box("LacMotherSolid",
                                frame_size.x() + 50.*mm,
                                frame_size.y() + 50.*mm,
                                frame_size.z() + 50.*mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "LacMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(- ra2 - m_rotation_angle);
  auto pos = (gGeom.GetGlobalPosition("KURAMA") +
              gGeom.GetGlobalPosition("LAC"));
  pos.rotateY(m_rotation_angle);
  new G4PVPlacement(rot, pos, mother_lv,
                    "LacMotherPV", m_world_lv, false, 0);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Frame
  auto frame_solid = new G4Box("LacFrameSolid", frame_size.x(),
                               frame_size.y(), frame_size.z());
  auto frame_lv = new G4LogicalVolume(frame_solid,
                                      m_material_map["Air"],
                                      "LacFrameLV");
  pos.setMag(0.);
  new G4PVPlacement(nullptr, pos, frame_lv,
                    "LacFramePV", mother_lv, false, 0);
  // Radiator
  auto radiator_solid = new G4Box("LacRadiatorSolid", radiator_size.x(),
                                  radiator_size.y(), radiator_size.z());
  auto radiator_lv = new G4LogicalVolume(radiator_solid,
                                         m_material_map["SilicaAerogelLAC"],
                                         "LacRadiatorLV");
  radiator_lv->SetSensitiveDetector(lac_sd);
  radiator_lv->SetVisAttributes(G4Color::Magenta());
  pos.set(0., 0., -frame_size.z() + radiator_size.z() + 10.*mm);
  new G4PVPlacement(nullptr, pos, radiator_lv,
                    "LacRadiatorPV", frame_lv, false, 0);
  // Mirror
  const G4double mirror_thickness = 1.*mm/2.;
  const G4double mirror_space = 20.*mm;
  const G4ThreeVector triangle_size(710.59*mm, frame_size.y(), 500.*mm);
  const G4double mirror_angle = std::atan2(triangle_size.z(),
                                           triangle_size.x());
  const G4ThreeVector mirror1_size((frame_size.x() - triangle_size.x())/2.,
                                   triangle_size.y(), mirror_thickness);
  const G4ThreeVector mirror2_size(std::hypot(triangle_size.x(),
                                              triangle_size.z())/2.,
                                   triangle_size.y(),
                                   mirror_thickness);
  auto mirror1_solid = new G4Box("LacMirror1Solid", mirror1_size.x(),
                                 mirror1_size.y(), mirror1_size.z());
  auto mirror1_lv = new G4LogicalVolume(mirror1_solid,
                                        m_material_map["Aluminum"],
					"LacMirror1LV");
  auto mirror2_solid = new G4Box("LacMirror2Solid", mirror2_size.x(),
                                 mirror2_size.y(), mirror2_size.z());
  auto mirror2_lv = new G4LogicalVolume(mirror2_solid,
                                        m_material_map["Aluminum"],
					"LacMirror2LV");
  for(G4int i=0; i<2; ++i){
    pos.set((triangle_size.x() + mirror1_size.x()) * (i*2 - 1),
            0., frame_size.z() - mirror_space);
    new G4PVPlacement(nullptr, pos, mirror1_lv,
                      "LacMirrorPV", frame_lv, false, 0);
    pos.set(triangle_size.x()/2 * (i*2 - 1),
            0., frame_size.z() - triangle_size.z()/2 - mirror_space);
    rot = new G4RotationMatrix;
    rot->rotateY(mirror_angle * (i*2 - 1));
    new G4PVPlacement(rot, pos, mirror2_lv,
                      "LacMirrorPV", frame_lv, false, 0);
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructSCH()
{
  using CLHEP::mm;
  const auto& sch_pos = (gGeom.GetGlobalPosition("KURAMA") +
                         gGeom.GetGlobalPosition("SCH"));
  const auto& sch_size = gSize.GetSize("SchSeg") * 0.5 * mm;
  const G4double dXdW = gGeom.GetWirePitch("SCH") * mm;
  G4LogicalVolume* sch_lv[NumOfSegSCH];
  auto sch_solid = new G4Box("SchSolid", sch_size.x(),
                             sch_size.y(), sch_size.z());
  for(G4int i=0; i<NumOfSegSCH; ++i){
    sch_lv[i] = new G4LogicalVolume(sch_solid, m_material_map["Scintillator"],
                                    Form("SchSeg%dLV", i), 0, 0, 0);
    sch_lv[i]->SetVisAttributes(G4Colour::Cyan());
    G4double ipos_x = dXdW * (i - (NumOfSegSCH - 1)/2.);
    G4ThreeVector pos(sch_pos.x() + ipos_x,
                      sch_pos.y(),
                      sch_pos.z() + 1.*mm*(1- 2*(i%2)));
    // pos.rotateY(m_rotation_angle);
    new G4PVPlacement(m_rotation_matrix, pos, sch_lv[i],
                      Form("SchSeg%dPV", i), m_world_lv, false, i);
  }
  auto schSD = new SCHSD("/SCH");
  G4SDManager::GetSDMpointer()->AddNewDetector(schSD);
  for(G4int i = 0; i<NumOfSegSCH; ++i){
    sch_lv[i]->SetSensitiveDetector(schSD);
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructSDC1()
{
  using CLHEP::mm;
  if(!sdc_sd){
    sdc_sd = new SDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector(sdc_sd);
  }
  const auto& sdc1_pos = (gGeom.GetGlobalPosition("KURAMA") +
                          (gGeom.GetGlobalPosition("SDC1-V1") +
                           gGeom.GetGlobalPosition("SDC1-U2")) * 0.5);
  const auto& frame_size = gSize.GetSize("Sdc1Frame") * 0.5 * mm;
  const auto& drift_size = gSize.GetSize("Sdc1Drift") * 0.5 * mm;
  auto sdc1_solid = new G4Box("Sdc1Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc1_lv = new G4LogicalVolume(sdc1_solid, m_material_map["Argon"],
                                     "Sdc1LV", 0, 0, 0);
  sdc1_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(m_rotation_matrix, sdc1_pos,
                    sdc1_lv, "Sdc1PV", m_world_lv, false, 0);
  auto sdc1pl_solid = new G4Box("Sdc1PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc1V1", "Sdc1V2", "Sdc1X1",
			    "Sdc1X2", "Sdc1U1", "Sdc1U2" };
  for(G4int i=0; i<NumOfLayersSDC1; ++i){
    G4ThreeVector pos;
    switch (i) {
    case 0:
      pos.setZ(-22.5985*mm);
      break;
    case 1:
      pos.setZ(-17.4015*mm);
      break;
    case 2:
      pos.setZ(-2.5985*mm);
      break;
    case 3:
      pos.setZ(2.5985*mm);
      break;
    case 4:
      pos.setZ(17.4015*mm);
      break;
    case 5:
      pos.setZ(22.5985*mm);
      break;
    }
    auto sdc1pl_lv = new G4LogicalVolume(sdc1pl_solid,
                                         m_material_map["Argon"],
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc1pl_lv->SetSensitiveDetector(sdc_sd);
    new G4PVPlacement(nullptr, pos, sdc1pl_lv, plane_name[i] + "PV",
                      sdc1_lv, false, 101+i);
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructSDC3()
{
  using CLHEP::mm;
  if(!sdc_sd){
    sdc_sd = new SDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector(sdc_sd);
  }
  const auto& sdc3_pos = (gGeom.GetGlobalPosition("KURAMA") +
                          (gGeom.GetGlobalPosition("SDC3-X1") +
                           gGeom.GetGlobalPosition("SDC3-Y2")) * 0.5);
  const auto& frame_size = gSize.GetSize("Sdc3Frame") * 0.5 * mm;
  const auto& drift_size = gSize.GetSize("Sdc3Drift") * 0.5 * mm;
  auto sdc3_solid = new G4Box("Sdc3Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc3_lv = new G4LogicalVolume(sdc3_solid, m_material_map["Argon"],
                                     "Sdc3LV", 0, 0, 0);
  sdc3_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(m_rotation_matrix, sdc3_pos,
                    sdc3_lv, "Sdc3PV", m_world_lv, false, 0);
  auto sdc3pl_solid = new G4Box("Sdc3PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc3X1", "Sdc3X2", "Sdc3Y1", "Sdc3Y2" };
  for(G4int i=0; i<NumOfLayersSDC3; ++i){
    G4ThreeVector pos;
    switch (i) {
    case 0:
      pos.setZ(-16.0*mm);
      break;
    case 1:
      pos.setZ(-8.206*mm);
      break;
    case 2:
      pos.setZ(8.206*mm);
      break;
    case 3:
      pos.setZ(16.0*mm);
      break;
    }
    auto sdc3pl_lv = new G4LogicalVolume(sdc3pl_solid,
                                         m_material_map["Argon"],
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc3pl_lv->SetSensitiveDetector(sdc_sd);
    new G4PVPlacement(nullptr, pos, sdc3pl_lv, plane_name[i] + "PV",
                      sdc3_lv, false, 301+i);
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructSDC4()
{
  using CLHEP::mm;
  if(!sdc_sd){
    sdc_sd = new SDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector(sdc_sd);
  }
  const auto& sdc4_pos = (gGeom.GetGlobalPosition("KURAMA") +
                          (gGeom.GetGlobalPosition("SDC4-Y1") +
                           gGeom.GetGlobalPosition("SDC4-X2")) * 0.5);
  const auto& frame_size = gSize.GetSize("Sdc4Frame") * 0.5 * mm;
  const auto& drift_size = gSize.GetSize("Sdc4Drift") * 0.5 * mm;
  auto sdc4_solid = new G4Box("Sdc4Solid", frame_size.x(),
                              frame_size.y(), frame_size.z());
  auto sdc4_lv = new G4LogicalVolume(sdc4_solid, m_material_map["Argon"],
                                     "Sdc4LV", 0, 0, 0);
  sdc4_lv->SetVisAttributes(G4Colour::Green());
  new G4PVPlacement(m_rotation_matrix, sdc4_pos,
                    sdc4_lv, "Sdc4PV", m_world_lv, false, 0);
  auto sdc4pl_solid = new G4Box("Sdc4PlSolid", drift_size.x(),
                                drift_size.y(), drift_size.z());
  G4String plane_name[] = { "Sdc4Y1", "Sdc4Y2", "Sdc4X1", "Sdc4X2" };
  for(G4int i=0; i<NumOfLayersSDC4; ++i){
    G4ThreeVector pos;
    switch (i) {
    case 0:
      pos.setZ(-34.868*mm);
      break;
    case 1:
      pos.setZ(-17.547*mm);
      break;
    case 2:
      pos.setZ(17.547*mm);
      break;
    case 3:
      pos.setZ(34.868*mm);
      break;
    }
    auto sdc4pl_lv = new G4LogicalVolume(sdc4pl_solid,
                                         m_material_map["Argon"],
                                         plane_name[i] + "LV", 0, 0, 0);
    sdc4pl_lv->SetSensitiveDetector(sdc_sd);
    new G4PVPlacement(nullptr, pos, sdc4pl_lv, plane_name[i] + "PV",
                      sdc4_lv, false, 401+i);
  }
}

//_____________________________________________________________________________
void
DetectorConstruction::ConstructWC()
{
  using CLHEP::mm;
  using CLHEP::deg;
  const auto& ra2 = gGeom.GetRotAngle2("WC") * deg;
  const auto& half_size_In = gSize.GetSize("WcSegIn") * 0.5 * mm;
  const auto& half_size_Out = gSize.GetSize("WcSegOut") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("WC");
  auto wcSD = new WCSD("/WC");
  wcSD->SetRefractiveIndex(1.33);
  G4SDManager::GetSDMpointer()->AddNewDetector(wcSD);
  // Mother
  auto mother_solid = new G4Box("WcMotherSolid",
                                half_size_Out.x()*NumOfSegWC + 200.*mm,
                                half_size_Out.y() + 200.*mm,
                                half_size_Out.z()*2 + 200.*mm);
  // half_size_Out.x()*NumOfSegWC + 50.*mm,
  // half_size_Out.y() + 50.*mm,
  // half_size_Out.z()*2 + 50.*mm);

  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "WcMotherLV");
  auto rot = new G4RotationMatrix;
  rot->rotateY(- ra2 - m_rotation_angle);
  auto pos = (gGeom.GetGlobalPosition("KURAMA") +
              gGeom.GetGlobalPosition("WC"));
  pos.rotateY(m_rotation_angle);
  new G4PVPlacement(rot, pos, mother_lv,
                    "WcMotherPV", m_world_lv, false, 0);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Segment
  auto segment_solid = new G4Box("WcSegmentSolid", half_size_In.x(),
                                 half_size_In.y(), half_size_In.z());
  auto segment_lv = new G4LogicalVolume(segment_solid,
                                        m_material_map["Water"],
                                        "WcSegmentLV");

  auto WCContainer     = new G4Box("WCContainer",
				   half_size_Out.x(),
				   half_size_Out.y(),
				   half_size_Out.z());
  auto WCContainer_gap = new G4Box("WCContainer_gap",
				   half_size_In.x(),
				   half_size_In.y(),
				   half_size_In.z());
  //G4RotationMatrix* rot_wccontainer_gap;
  auto rot_wccontainer_gap = new G4RotationMatrix;
  G4ThreeVector pos_wccontainer_gap(0.0, 0.0 ,0.0);
  auto solid_WCContainer
    = new G4SubtractionSolid("solid_WCContainer",
   			     WCContainer, WCContainer_gap,
   			     rot_wccontainer_gap,
			     pos_wccontainer_gap);

  auto logWCContainer = new G4LogicalVolume(solid_WCContainer,
					    m_material_map["Acrylic"],
					    "logWCContainer");

  for(G4int i=0; i<NumOfSegWC; ++i){


    pos = G4ThreeVector((-NumOfSegWC/2 + i)*pitch,
                        0.0,
                        2.*(- i%2 + 0.5)*half_size_Out.z());
    //for Vessel
    //    logWCContainer->SetVisAttributes(G4Colour::White());
    logWCContainer->SetVisAttributes(G4Colour::Cyan());
    new G4PVPlacement(nullptr, pos, logWCContainer,
                      "WcSegmentContainerPV", mother_lv, false, i);
    //for Water
    segment_lv->SetVisAttributes(G4Colour::Cyan());
    segment_lv->SetSensitiveDetector(wcSD);
    new G4PVPlacement(nullptr, pos, segment_lv,
                      "WcSegmentPV", mother_lv, false, i);

  }
  /*
  //Temporary!!!!!!
  //WC frame (test)
  //double Alcut_z =0.*mm;
  double Alcut_z =(1000.-(250.+403.))*mm;

  auto Alframe1     = new G4Box("Alframe1",
  80.*mm/2.,
  80.*mm/2.,
  2000.*mm/2. - Alcut_z/2.);

  auto Alframe2     = new G4Box("Alframe2",
  80.*mm/2.,
  80.*mm/2.,
  750.*mm/2.);



  auto Alframe1_lv = new G4LogicalVolume(Alframe1,
  m_material_map["Aluminum"],
  "Alframe1_LV");
  auto Alframe2_lv = new G4LogicalVolume(Alframe2,
  m_material_map["Aluminum"],
  "Alframe2_LV");

  auto pos_frame1_1 = G4ThreeVector(-pitch/2.-3000.*mm/2.,
  -2000.*mm,
  -1.*Alcut_z/2.);

  auto pos_frame1_2 = G4ThreeVector(-pitch/2.+3000.*mm/2.,
  -2000.*mm,
  -1.*Alcut_z/2.);

  Alframe1_lv->SetVisAttributes(G4Colour::White());
  new G4PVPlacement(nullptr, pos_frame1_1, Alframe1_lv,
  "WcFramePV", mother_lv, false, 0);

  new G4PVPlacement(nullptr, pos_frame1_2, Alframe1_lv,
  "WcFramePV", mother_lv, false, 1);


  auto pos_frame2_1 = G4ThreeVector(-pitch/2.-3000.*mm/2.,
  -1000.*mm,
  750/2.);

  auto pos_frame2_2 = G4ThreeVector(-pitch/2.+3000.*mm/2.,
  -1000.*mm,
  750/2.);

  Alframe2_lv->SetVisAttributes(G4Colour::Red());
  new G4PVPlacement(nullptr, pos_frame2_1, Alframe2_lv,
  "WcFramePV", mother_lv, false, 3);

  new G4PVPlacement(nullptr, pos_frame2_2, Alframe2_lv,
  "WcFramePV", mother_lv, false, 4);

  //Temporary!!!!!!
  */


}
