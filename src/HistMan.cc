// -*- C++ -*-

#include "HistMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <G4ThreeVector.hh>

#include <TSystem.h>

#include "FuncName.hh"

//_____________________________________________________________________________
HistMan::HistMan()
  : m_is_ready(false),
    m_file_name(),
    m_param_map()
{
}

//_____________________________________________________________________________
HistMan::~HistMan()
{
}

//_____________________________________________________________________________
Bool_t
HistMan::Initialize()
{
  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    G4cerr << "#E " << FUNC_NAME << " "
           << "No such parameter file : " << m_file_name << std::endl;
    return false;
  }

  G4String line;
  while(ifs.good() && std::getline(ifs, line)){
    if(line.empty() || line[0]=='#') continue;
    std::istringstream input_line(line);

    TString first_param;
    input_line >> first_param;

    TString key = first_param;
    ParamArray param_array;
    Double_t   param;
    while(input_line >> param){
      param_array.push_back(param);
    }
    m_param_map[key] = param_array;
  }

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
HistMan::Initialize(const TString& filename)
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
const std::vector<Double_t>&
HistMan::Get(const TString& key) const
{
  PIterator itr = m_param_map.find(key);
  if(itr==m_param_map.end()){
    Print();
    TString msg(FUNC_NAME+" No such key : "+key);
    G4cout << msg << G4endl;
    gSystem->Exit(1);
    // throw std::invalid_argument(msg);
  }
  return itr->second;
}

//_____________________________________________________________________________
void
HistMan::Print() const
{
  G4cout << FUNC_NAME << std::endl;

  const G4int w = 20;
  PIterator itr, end=m_param_map.end();
  for(itr=m_param_map.begin(); itr!=end; ++itr){
    G4cout << " key = " << std::setw(w) << std::left
           << itr->first << itr->second.size() << " : ";
    for(G4int i=0, n=itr->second.size(); i<n; ++i){
      G4cout << std::setw(5) << std::right
             << itr->second.at(i) << " ";
    }
    G4cout << std::endl;
  }
}
