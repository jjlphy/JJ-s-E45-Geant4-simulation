// -*- C++ -*-

#include "JamMan.hh"

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

#define SKIP_NP0 1

//_____________________________________________________________________________
void
JamInfo::Print() const
{
  PrintHelper helper(4, std::ios::fixed, G4cout);
  const G4int w = 8;
  G4cout << "   np=" << np << G4endl;
  for(G4int i=0; i<np; ++i){
    G4cout << "   - " << i << " "
	   << "pid=" << std::setw(w) << pid[i] << " "
	   << "px=" << std::setw(w) << px[i] << " "
	   << "py=" << std::setw(w) << py[i] << " "
	   << "pz=" << std::setw(w) << pz[i] << G4endl;
  }
}

//_____________________________________________________________________________
JamMan::JamMan()
  : m_is_ready(false),
    m_file_name(),
    m_file(),
    m_tree(),
    m_event(new JamInfo),
    m_n_event()
{
}

//_____________________________________________________________________________
JamMan::~JamMan()
{
  if(m_file && m_file->IsOpen())
    m_file->Close();
}

//_____________________________________________________________________________
G4bool
JamMan::Initialize()
{
  if(m_file_name.empty())
    return true;

  m_file = new TFile(m_file_name);
  m_tree = dynamic_cast<TTree*>(m_file->Get("tree"));

  if(!m_file->IsOpen() || !m_tree)
    return false;

  m_event = new JamInfo;
  m_file = new TFile(m_file_name);
  m_tree = dynamic_cast<TTree*>(m_file->Get("tree"));
  m_tree->SetBranchAddress("np", &m_event->np);
  m_tree->SetBranchAddress("pid", m_event->pid);
  m_tree->SetBranchAddress("px", m_event->px);
  m_tree->SetBranchAddress("py", m_event->py);
  m_tree->SetBranchAddress("pz", m_event->pz);
  m_tree->SetBranchStatus("*", false);
  m_tree->SetBranchStatus("np", true);
  m_tree->SetBranchStatus("pid", true);
  m_tree->SetBranchStatus("px", true);
  m_tree->SetBranchStatus("py", true);
  m_tree->SetBranchStatus("pz", true);
  m_n_event = m_tree->GetEntries();
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
JamMan::Initialize(const G4String& filename)
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
JamInfo*
JamMan::Get() const
{
  return Get(G4RandFlat::shootInt(m_n_event));
}

//_____________________________________________________________________________
JamInfo*
JamMan::Get(Int_t i) const
{
  auto curr_file = gFile;
  m_file->cd();
  m_tree->GetEntry(i);
  if(curr_file)
    curr_file->cd();
  return m_event;
}

//_____________________________________________________________________________
void
JamMan::Print() const
{
  PrintHelper helper(4, std::ios::fixed, G4cout);

  G4cout << FUNC_NAME << G4endl;
  if(m_event)
    m_event->Print();
  G4cout << "   n_event = " << m_n_event << G4endl;
}
