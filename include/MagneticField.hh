// -*- C++ -*-

#ifndef MAGNETIC_FIELD_HH
#define MAGNETIC_FIELD_HH

#include <vector>
#include <map>
#include <globals.hh>
#include <G4MagneticField.hh>
#include <G4String.hh>
#include <G4ThreeVector.hh>

class FieldMap;

//_____________________________________________________________________________
class MagnetInfo
{
public:
  static G4String ClassName();
  MagnetInfo(G4String n="")
    : name(n)
  {}
  enum EMagnetType { kDipole, kQuadrupole, kShs, kKurama, NMagnetType };

public:
  G4int    type;
  G4String name;
  G4double b0;        // magnetic field
  G4double l;         // length
  G4double ra1;       // rotation angle around z-axis
  G4ThreeVector pos;  // center position
  G4ThreeVector size; // half width, height, and thickness
  // for Dipole
  G4double rho;   // bending radius
  G4double bend;  // bending angle
  G4double alpha; // angle between conetral track and entrance pole
  G4double beta;  // angle between conetral track and exit pole
  // for Quadrupole
  G4double a0; // aperture

public:
  G4bool CalcK18Field(const G4ThreeVector& point, G4double* bfield) const;
};

//_____________________________________________________________________________
inline G4String
MagnetInfo::ClassName()
{
  static G4String s_name("MagnetInfo");
  return s_name;
}

//_____________________________________________________________________________
class MagneticField : public G4MagneticField
{
public:
  static G4String ClassName();
  MagneticField();
  ~MagneticField();

private:
  MagneticField(const MagneticField&);
  MagneticField& operator =(const MagneticField&);

private:
  typedef std::map<G4String, MagnetInfo> MagnetMap;
  typedef std::vector< std::vector< std::vector<G4ThreeVector> > > Field;
  G4bool    m_is_ready;
  G4bool    m_k18_status;
  G4bool    m_kurama_status;
  G4bool    m_shs_status;
  FieldMap* m_kurama_field_map;
  FieldMap* m_shs_field_map;
  MagnetMap m_magnet_map;

public:
  virtual void GetFieldValue(const G4double Point[4], G4double* Bfield) const;

public:
  void   AddMagnetInfo(const MagnetInfo& mag);
  const G4ThreeVector& GetSizeShsField() const;
  G4bool GetStatusK18Field() const { return m_k18_status; }
  G4bool GetStatusKuramaField() const { return m_kurama_status; }
  G4bool GetStatusShsField() const { return m_shs_status; }
  G4bool Initialize();
  G4bool IsReady() const { return m_is_ready; }
  void   SetKuramaFieldMap(G4String map);
  void   SetShsFieldMap(G4String map);
  void   SetStatusK18Field(G4bool flag=true){ m_k18_status = flag; }
  void   SetStatusKuramaField(G4bool flag=true){ m_kurama_status = flag; }
  void   SetStatusShsField(G4bool flag=true){ m_shs_status = flag; }

};

//_____________________________________________________________________________
inline G4String
MagneticField::ClassName()
{
  static G4String s_name("MagneticField");
  return s_name;
}

#endif
