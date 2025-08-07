// -*- C++ -*-

#include "DetSizeMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <G4ThreeVector.hh>

#include "FuncName.hh"

//_____________________________________________________________________________
DetSizeMan::DetSizeMan()
  : m_is_ready(false),
    m_file_name(),
    m_param_map()
{
}

//_____________________________________________________________________________
DetSizeMan::~DetSizeMan()
{
}

//_____________________________________________________________________________
G4bool
DetSizeMan::Initialize()
{
  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    G4cerr << "#E " << FUNC_NAME << " "
           << "No such parameter file : " << m_file_name << std::endl;
    return false;
  }

  G4String line;
  while(ifs.good() && std::getline(ifs, line)){
    if(line[0]=='#') continue;
    std::istringstream input_line(line);

    G4String first_param;
    input_line >> first_param;

    G4String key = first_param;
    ParamArray param_array;
    G4double   param;
    while(input_line >> param){
      param_array.push_back(param);
    }
    m_param_map[key] = param_array;
  }

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
DetSizeMan::Initialize(const G4String& filename)
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
G4double
DetSizeMan::Get(const G4String& key, G4int i) const
{
  std::stringstream param;
  param << key << "(" << i << ")";

  PIterator itr = m_param_map.find(key);

  if(itr==m_param_map.end() ||
     i+1 > (G4int)itr->second.size()){
    Print();
    G4String msg(FUNC_NAME+" No such key : "+key);
    msg += "(i=" + std::to_string(i) + ")";
    throw std::invalid_argument(msg);
  }

  return itr->second.at(i);
}

//_____________________________________________________________________________
G4ThreeVector
DetSizeMan::GetSize(const G4String& key) const
{
  return G4ThreeVector(Get(key, G4ThreeVector::X),
                       Get(key, G4ThreeVector::Y),
                       Get(key, G4ThreeVector::Z));
}

//_____________________________________________________________________________
void
DetSizeMan::Print() const
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
