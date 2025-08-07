// -*- C++ -*-

#ifndef DIFF_CROSS_SECTION_MAN_HH
#define DIFF_CROSS_SECTION_MAN_HH

#include <string>
#include <vector>

#include <globals.hh>

#include <TSpline.h>

#include <Rtypes.h>

class TFile;

//_____________________________________________________________________________
class DiffCrossSectionMan
{
public:
  static G4String ClassName();
  static DiffCrossSectionMan& GetInstance();
  ~DiffCrossSectionMan();

private:
  DiffCrossSectionMan();
  DiffCrossSectionMan(const DiffCrossSectionMan& );
  DiffCrossSectionMan& operator =(const DiffCrossSectionMan&);

private:
  G4bool        m_is_ready;
  G4String      m_file_name;
  TFile*        m_file;
  G4int         m_n_spline;
  std::vector<TSpline3*> m_spline_container;
  G4double      m_range_min;
  G4double      m_range_max;
  
public:
  G4bool               Initialize();
  G4bool               Initialize(const G4String& filename);
  G4bool               IsReady() const { return m_is_ready; }
  G4double             GetCoeff(G4int order, G4double mom_kaon) const;
  G4int                GetNumSpline() const { return m_n_spline; }
  
};

//_____________________________________________________________________________
inline G4String
DiffCrossSectionMan::ClassName()
{
  static G4String s_name("DiffCrossSectionMan");
  return s_name;
}

//_____________________________________________________________________________
inline DiffCrossSectionMan&
DiffCrossSectionMan::GetInstance()
{
  static DiffCrossSectionMan s_instance;
  return s_instance;
}

#endif
