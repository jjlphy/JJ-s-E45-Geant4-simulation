// -*- C++ -*-

#include "DiffCrossSectionMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TSpline.h>

//_____________________________________________________________________________
DiffCrossSectionMan::DiffCrossSectionMan()
  : m_is_ready(false),
    m_file_name(),
    m_file(),
    m_n_spline(),
    m_spline_container(),
    m_range_min(),
    m_range_max()
{
}

//_____________________________________________________________________________
DiffCrossSectionMan::~DiffCrossSectionMan()
{
}

//_____________________________________________________________________________
G4bool
DiffCrossSectionMan::Initialize()
{
  if(m_file_name.empty())
    return true;

  m_file = new TFile(m_file_name);

  if(!m_file->IsOpen())
    return false;

  m_n_spline = m_file->GetListOfKeys()->GetSize();
  for (G4int i = 0; i < m_n_spline; i++) {
    G4String spline_name = "A" + G4String(std::to_string(i)) + "_spline";
    TSpline3* spline = dynamic_cast<TSpline3*>(m_file->Get(spline_name.c_str()));
    if (spline) {
      m_spline_container.push_back(spline);
    } else {
      G4cerr << "Error: Failed to retrieve TSpline3: " << spline_name << G4endl;
    }
  }
  m_range_min = m_spline_container[0]->GetXmin();
  m_range_max = m_spline_container[0]->GetXmax();

  m_file->Close();
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
DiffCrossSectionMan::Initialize(const G4String& filename)
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
G4double
DiffCrossSectionMan::GetCoeff(G4int order, G4double mom_kaon) const
{
  if (order < 0 || order >= m_n_spline || !m_spline_container[order]) {
    G4cerr << "Error: Invalid spline index: " << order << G4endl;
    return 0.0;
  }
  
  if (mom_kaon < m_range_min) {
    return m_spline_container[order]->Eval(m_range_min);
  } else if (mom_kaon > m_range_max) {
    return m_spline_container[order]->Eval(m_range_max);
  } else {
    return m_spline_container[order]->Eval(mom_kaon);
  }
}
