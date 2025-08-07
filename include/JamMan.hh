// -*- C++ -*-

#ifndef JAM_MAN_HH
#define JAM_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <globals.hh>
#include <G4ThreeVector.hh>

#include <Rtypes.h>

#include "AnaManager.hh"

class TFile;
class TTree;

//_____________________________________________________________________________
struct JamInfo
{
  Int_t np;
  Int_t pid[MaxHits];
  Double_t px[MaxHits]; // [GeV/c]
  Double_t py[MaxHits]; // [GeV/c]
  Double_t pz[MaxHits]; // [GeV/c]
  void Print() const;
};

//_____________________________________________________________________________
class JamMan
{
public:
  static G4String ClassName();
  static JamMan& GetInstance();
  ~JamMan();

private:
  JamMan();
  JamMan(const JamMan& );
  JamMan& operator =(const JamMan&);

private:
  G4bool     m_is_ready;
  G4String   m_file_name;
  TFile*     m_file;
  TTree*     m_tree;
  JamInfo*   m_event;
  G4int      m_n_event;

public:
  JamInfo* Get() const;
  JamInfo* Get(Int_t i) const;
  G4bool   Initialize();
  G4bool   Initialize(const G4String& filename);
  G4bool   IsReady() const { return m_is_ready; }
  void     Print() const;
};

//_____________________________________________________________________________
inline G4String
JamMan::ClassName()
{
  static G4String s_name("JamMan");
  return s_name;
}

//_____________________________________________________________________________
inline JamMan&
JamMan::GetInstance()
{
  static JamMan s_instance;
  return s_instance;
}

#endif
