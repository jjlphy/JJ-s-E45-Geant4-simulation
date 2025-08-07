// -*- C++ -*-

#ifndef DET_SIZE_MAN_HH
#define DET_SIZE_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <globals.hh>
#include <G4ThreeVector.hh>

//_____________________________________________________________________________
class DetSizeMan
{
public:
  static G4String    ClassName();
  static DetSizeMan& GetInstance();
  ~DetSizeMan();

private:
  DetSizeMan();
  DetSizeMan(const DetSizeMan& );
  DetSizeMan& operator =(const DetSizeMan&);

private:
  typedef std::vector<G4double>          ParamArray;
  typedef std::map<G4String, ParamArray> ParamMap;
  typedef ParamMap::const_iterator       PIterator;
  G4bool   m_is_ready;
  G4String m_file_name;
  ParamMap m_param_map;

public:
  G4bool        Initialize();
  G4bool        Initialize(const G4String& filename);
  G4bool        IsReady() const { return m_is_ready; }
  G4double      Get(const G4String& key, G4int i=0) const;
  G4ThreeVector GetSize(const G4String& key) const;
  void          Print() const;
};

//_____________________________________________________________________________
inline G4String
DetSizeMan::ClassName()
{
  static G4String s_name("DetSizeMan");
  return s_name;
}

//_____________________________________________________________________________
inline DetSizeMan&
DetSizeMan::GetInstance()
{
  static DetSizeMan s_instance;
  return s_instance;
}

#endif
