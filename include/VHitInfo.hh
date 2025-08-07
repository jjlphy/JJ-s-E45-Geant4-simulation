// -*- C++ -*-

#ifndef VHIT_INFO_HH
#define VHIT_INFO_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>

class G4Step;
class TParticle;

//_____________________________________________________________________________
class VHitInfo
{
public:
  static G4String ClassName();
  VHitInfo(const G4String& name, G4Step* step=nullptr,
           G4bool use_center_point=false);
  ~VHitInfo();

private:
  G4String      m_detector_name;
  G4String      m_particle_name;
  G4ThreeVector m_position;
  G4ThreeVector m_momentum;
  G4double      m_time;
  G4double      m_energy_deposit;
  G4int         m_track_id;
  G4int         m_pdg_encoding;
  G4int         m_copy_number;
  G4double      m_step_length;
  G4double      m_mass;
  G4double      m_charge;
  G4int         m_parent_id;
  G4double      m_track_length;
  G4ThreeVector m_vertex_position;
  G4ThreeVector m_vertex_momentum;
  G4double      m_vertex_kinetic_energy; // sqrt(p^2+m^2)-m
  TParticle*    m_particle; // for tree

public:
  void     AddEnergyDeposit(G4double de){ m_energy_deposit += de; }
  const G4String& GetDetectorName() const { return m_detector_name; }
  const G4String& GetParticleName() const { return m_particle_name; }
  const G4ThreeVector& GetPosition() const { return m_position; }
  const G4ThreeVector& GetMomentum() const { return m_momentum; }
  G4double GetTime() const { return m_time; }
  G4double GetEnergyDeposit() const { return m_energy_deposit; }
  G4int    GetTrackID() const { return m_track_id; }
  G4int    GetPDGEncoding() const { return m_pdg_encoding; }
  G4int    GetCopyNumber() const { return m_copy_number; }
  G4double GetMass() const { return m_mass; }
  G4double GetCharge() const { return m_charge; }
  G4int    GetParentID() const { return m_parent_id; }
  G4double GetTrackLength() const { return m_track_length; }
  const G4ThreeVector& GetVertexPosition() const
  { return m_vertex_position; }
  const G4ThreeVector& GetVertexMomentum() const
  { return m_vertex_momentum; }
  G4double      GetVertexKineticEnergy() const
  { return m_vertex_kinetic_energy; }
  void          Print() const;
  TParticle* GetParticle() const { return m_particle; }
};

//_____________________________________________________________________________
inline G4String
VHitInfo::ClassName()
{
  static G4String s_name("VHitInfo");
  return s_name;
}

#endif
