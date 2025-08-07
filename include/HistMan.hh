// -*- C++ -*-

#ifndef HIST_MAN_HH
#define HIST_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <G4String.hh>

#include <TString.h>

//_____________________________________________________________________________
class HistMan
{
public:
  static G4String ClassName();
  static HistMan& GetInstance();
  ~HistMan();

private:
  HistMan();
  HistMan(const HistMan& );
  HistMan& operator =(const HistMan&);

private:
  typedef std::vector<Double_t>         ParamArray;
  typedef std::map<TString, ParamArray> ParamMap;
  typedef ParamMap::const_iterator      PIterator;
  Bool_t   m_is_ready;
  TString  m_file_name;
  ParamMap m_param_map;

public:
  const std::vector<Double_t>& Get(const TString& key) const;
  Bool_t        Initialize();
  Bool_t        Initialize(const TString& filename);
  Bool_t        IsReady() const { return m_is_ready; }
  void          Print() const;
};

//_____________________________________________________________________________
inline G4String
HistMan::ClassName()
{
  static G4String s_name("HistMan");
  return s_name;
}

//_____________________________________________________________________________
inline HistMan&
HistMan::GetInstance()
{
  static HistMan s_instance;
  return s_instance;
}

#endif
