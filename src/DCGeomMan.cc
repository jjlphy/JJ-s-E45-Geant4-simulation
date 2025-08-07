// -*- C++ -*-

#include "DCGeomMan.hh"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

#include "DCGeomRecord.hh"
#include "FuncName.hh"

//_____________________________________________________________________________
DCGeomMan::DCGeomMan()
  : m_is_ready(false),
    m_file_name(),
    m_container(),
    m_detector_id_map(),
    m_global_z_map(),
    m_local_z_map()
{
}

//_____________________________________________________________________________
DCGeomMan::~DCGeomMan()
{
  Clear();
}

//_____________________________________________________________________________
Int_t
DCGeomMan::CalcWireNumber(Int_t lnum, Double_t pos) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record){
    return record->WireNumber(pos);
  }
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
Int_t
DCGeomMan::CalcWireNumber(const TString& key, Double_t wire) const
{
  return CalcWireNumber(GetDetectorId(key), wire);
}

//_____________________________________________________________________________
Double_t
DCGeomMan::CalcWirePosition(Int_t lnum, Double_t wire) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record){
    return record->WirePos(wire);
  }
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
Double_t
DCGeomMan::CalcWirePosition(const TString& key, Double_t wire) const
{
  return CalcWirePosition(GetDetectorId(key), wire);
}

//_____________________________________________________________________________
void
DCGeomMan::Clear()
{
  for(auto&& itr : m_container){
    delete itr.second;
  }
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetLocalZ(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record){
    return record->Length();
  }else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetLocalZ(const TString& key) const
{
  return GetLocalZ(GetDetectorId(key));
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetResolution(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->Resolution();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetResolution(const TString& key) const
{
  return GetResolution(GetDetectorId(key));
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetRotAngle1(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->RotationAngle1();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetRotAngle1(const TString& key) const
{
  return GetRotAngle1(GetDetectorId(key));
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetRotAngle2(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->RotationAngle2();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetRotAngle2(const TString& key) const
{
  return GetRotAngle2(GetDetectorId(key));
}

//_____________________________________________________________________________
const ThreeVector&
DCGeomMan::GetGlobalPosition(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->Pos();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
const ThreeVector&
DCGeomMan::GetGlobalPosition(const TString& key) const
{
  return GetGlobalPosition(GetDetectorId(key));
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetTiltAngle(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->TiltAngle();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetTiltAngle(const TString& key) const
{
  return GetTiltAngle(GetDetectorId(key));
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetWirePitch(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->WirePitch();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
Double_t
DCGeomMan::GetWirePitch(const TString& key) const
{
  return GetWirePitch(GetDetectorId(key));
}

//_____________________________________________________________________________
Bool_t
DCGeomMan::Initialize(const TString& file_name)
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::NormalVector(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->NormalVector();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::NormalVector(const TString& key) const
{
  return NormalVector(GetDetectorId(key));
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::UnitVector(Int_t lnum) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(record) return record->UnitVector();
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::UnitVector(const TString& key) const
{
  return UnitVector(GetDetectorId(key));
}

//_____________________________________________________________________________
const DCGeomRecord*
DCGeomMan::GetRecord(Int_t lnum) const
{
  DCGeomIterator itr = m_container.find(lnum);
  DCGeomIterator end = m_container.end();
  DCGeomRecord* record = 0;
  if(itr!=end)
    record = itr->second;
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
  return record;
}

//_____________________________________________________________________________
const DCGeomRecord*
DCGeomMan::GetRecord(const TString& key) const
{
  return GetRecord(GetDetectorId(key));
}

//_____________________________________________________________________________
Bool_t
DCGeomMan::Initialize()
{
  if(m_is_ready){
    std::cerr << "#W " << FUNC_NAME
		<< " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    std::cerr << "#E " << FUNC_NAME
		<< " file open fail : " << m_file_name << std::endl;
    return false;
  }

  Clear();

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    if(line[0]=='#') continue;
    std::istringstream iss(line.Data());
    Int_t id; TString name;
    Double_t gx, gy, gz, ta, ra1, ra2, l, res, w0, dd, ofs;
    if(iss >> id >> name >> gx >> gy >> gz >> ta >> ra1 >> ra2
	>> l >> res >> w0 >> dd >> ofs){
      DCGeomRecord* record =
	new DCGeomRecord(id, name, gx, gy, gz, ta, ra1, ra2,
			  l, res, w0, dd, ofs);
      if(m_container[id]){
	std::cerr << "#W " << FUNC_NAME << " "
		    << "duplicated key is deleted : " << id << std::endl;
	m_container[id]->Print();
	delete m_container[id];
      }
      m_container[id] = record;

      m_detector_id_map[name] = id;
      m_global_z_map[name]    = gz;
      m_local_z_map[name]     = l;
    }else{
      std::cerr << "#E " << FUNC_NAME << " "
		  << "invalid format : " << line << std::endl;
    }
  }

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
std::vector<Int_t>
DCGeomMan::GetDetectorIDList() const
{
  std::vector<Int_t> vlist;
  vlist.reserve(m_container.size());
  DCGeomIterator itr, end=m_container.end();
  for(itr=m_container.begin(); itr!=end; ++itr){
    vlist.push_back(itr->first);
  }

  return vlist;
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::Local2GlobalPos(Int_t lnum, const ThreeVector& in) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }

  Double_t x = record->dxds()*in.x() + record->dxdt()*in.y()
    + record->dxdu()*in.z() + record->Pos().x();
  Double_t y = record->dyds()*in.x() + record->dydt()*in.y()
    + record->dydu()*in.z() + record->Pos().y();
  Double_t z = record->dzds()*in.x() + record->dzdt()*in.y()
    + record->dzdu()*in.z() + record->Pos().z();

  return ThreeVector(x, y, z);
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::Local2GlobalPos(const TString& key, const ThreeVector& in) const
{
  return Local2GlobalPos(GetDetectorId(key), in);
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::Global2LocalPos(Int_t lnum, const ThreeVector& in) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }

  Double_t x
    = record->dsdx()*(in.x()-record->Pos().x())
    + record->dsdy()*(in.y()-record->Pos().y())
    + record->dsdz()*(in.z()-record->Pos().z());
  Double_t y
    = record->dtdx()*(in.x()-record->Pos().x())
    + record->dtdy()*(in.y()-record->Pos().y())
    + record->dtdz()*(in.z()-record->Pos().z());
  Double_t z
    = record->dudx()*(in.x()-record->Pos().x())
    + record->dudy()*(in.y()-record->Pos().y())
    + record->dudz()*(in.z()-record->Pos().z());

  return ThreeVector(x, y, z);
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::Global2LocalPos(const TString& key, const ThreeVector& in) const
{
  return Global2LocalPos(GetDetectorId(key), in);
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::Local2GlobalDir(Int_t lnum, const ThreeVector& in) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }

  Double_t x = record->dxds()*in.x() + record->dxdt()*in.y()
    + record->dxdu()*in.z();
  Double_t y = record->dyds()*in.x() + record->dydt()*in.y()
    + record->dydu()*in.z();
  Double_t z = record->dzds()*in.x() + record->dzdt()*in.y()
    + record->dzdu()*in.z();

  return ThreeVector(x, y, z);
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::Local2GlobalDir(const TString& key, const ThreeVector& in) const
{
  return Local2GlobalDir(GetDetectorId(key), in);
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::Global2LocalDir(Int_t lnum, const ThreeVector& in) const
{
  const DCGeomRecord* record = GetRecord(lnum);
  if(!record){
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }

  Double_t x = record->dsdx()*in.x() + record->dsdy()*in.y()
    + record->dsdz()*in.z();
  Double_t y = record->dtdx()*in.x() + record->dtdy()*in.y()
    + record->dtdz()*in.z();
  Double_t z = record->dudx()*in.x() + record->dudy()*in.y()
    + record->dudz()*in.z();

  return ThreeVector(x, y, z);
}

//_____________________________________________________________________________
ThreeVector
DCGeomMan::Global2LocalDir(const TString& key, const ThreeVector& in) const
{
  return Global2LocalDir(GetDetectorId(key), in);
}

//_____________________________________________________________________________
Int_t
DCGeomMan::GetDetectorId(const TString& key) const
{
  DCGeomIterator itr, end = m_container.end();
  for(itr=m_container.begin(); itr!=end; ++itr){
    if (itr->second->Name() == key)
      return itr->second->Id();
  }
  G4String e(FUNC_NAME + Form(" No record : %s", key.Data()));
  throw std::invalid_argument(e);
}

//_____________________________________________________________________________
void
DCGeomMan::SetFileName(const TString& file_name)
{
  m_file_name = file_name;
}

//_____________________________________________________________________________
void
DCGeomMan::SetResolution(Int_t lnum, Double_t res)
{
  DCGeomRecord* record = m_container[lnum];
  if(record) {
    record->SetResolution(res);
  }
  else{
    G4String e(FUNC_NAME + Form(" No record : %d", lnum));
    throw std::invalid_argument(e);
  }
}

//_____________________________________________________________________________
void
DCGeomMan::SetResolution(const TString& key, Double_t res)
{
  SetResolution(GetDetectorId(key), res);
}
