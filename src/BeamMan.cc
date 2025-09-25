// -*- C++ -*-

#include "BeamMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TFile.h>
#include <TTree.h>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
void
BeamInfo::Print() const
{
  PrintHelper helper(4, std::ios::fixed, G4cout);
  G4cout << "   BeamInfo : "
         << "pos=" << pos << ",  "
         << "mom=" << mom << G4endl;
}

//_____________________________________________________________________________
BeamMan::BeamMan()
  : m_is_ready(false),
    m_file_name(),
    m_file(),
    m_param_array(),
    m_n_param()
{
}

//_____________________________________________________________________________
BeamMan::~BeamMan()
{
}

//_____________________________________________________________________________
G4bool
BeamMan::Initialize()
{
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const G4double p0 = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
  const G4double x0 = gGeom.GetGlobalPosition("BH2").x()*CLHEP::mm;
  const G4double z0 = gGeom.GetGlobalPosition("BH2").z()*CLHEP::mm;

  if(m_file_name.empty())
    return true;

  m_file = new TFile(m_file_name);
  TTree* tree = dynamic_cast<TTree*>(m_file->Get("tr"));

  if(!m_file->IsOpen() || !tree)
    return false;

  m_param_array.clear();
  BeamInfo beam;
  tree->SetBranchAddress("pointInx", &beam.x);
  tree->SetBranchAddress("pointIny", &beam.y);
  tree->SetBranchAddress("pointInz", &beam.z);
  tree->SetBranchAddress("pInx", &beam.px);
  tree->SetBranchAddress("pIny", &beam.py);
  tree->SetBranchAddress("pInz", &beam.pz);

// === [추가] 빔 위치 오프셋 (conf에서 읽음, 기본값 0) =======================
const G4double dx = gConf.Get<G4double>("BeamShiftX_mm") * CLHEP::mm;
const G4double dy = gConf.Get<G4double>("BeamShiftY_mm") * CLHEP::mm;
const G4double dz = gConf.Get<G4double>("BeamShiftZ_mm") * CLHEP::mm;
// ===========================================================================


  for(Long64_t i=0, n=tree->GetEntries(); i<n; ++i){
    tree->GetEntry(i);
    //beam.pos.set(beam.x + x0, beam.y, z0);
     //beam.pos.set(beam.x + x0, beam.y, beam.z);
    beam.pos.set(beam.x, beam.y, beam.z);
    beam.pos *= CLHEP::mm;
     // 오프셋 적용
    beam.pos += G4ThreeVector(dx, dy, dz);

    G4double scale = p0/0.907;
    beam.mom.set(beam.px * scale, beam.py * scale, beam.pz * scale);
    m_param_array.push_back(beam);
  }

  m_file->Close();
  m_n_param = m_param_array.size();
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
BeamMan::Initialize(const G4String& filename)
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
const BeamInfo&
BeamMan::Get() const
{
  if (m_is_ready && m_n_param > 0){
    return m_param_array.at(G4RandFlat::shootInt(m_n_param));
  } else {
    static BeamInfo nullinfo;
    return nullinfo;
  }
}

//_____________________________________________________________________________
void
BeamMan::Print() const
{
  PrintHelper helper(4, std::ios::fixed, G4cout);
  G4cout << FUNC_NAME << G4endl;
  for(const auto& b : m_param_array){
    G4cout << "   "
	   << "pos=" << b.pos << ",  "
	   << "mom=" << b.mom << G4endl;
  }
  G4cout << "   nparam = " << m_param_array.size() << G4endl;
}
