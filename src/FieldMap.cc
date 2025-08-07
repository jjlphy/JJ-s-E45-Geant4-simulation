// -*- C++ -*-

#include "FieldMap.hh"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <CLHEP/Units/SystemOfUnits.h>

#include "ConfMan.hh"
#include "FuncName.hh"

//______________________________________________________________________________
FieldMap::FieldMap(const G4String& file_name)
  : m_is_ready(false),
    m_file_name(file_name),
    m_b(),
    m_nx(0), m_ny(0), m_nz(0),
    m_xmin(0.), m_ymin(0.), m_zmin(0.),
    m_xmax(0.), m_ymax(0.), m_zmax(0.),
    m_dx(1.), m_dy(1.), m_dz(1.),
    m_value_calc(1.),
    m_value_nmr(1.),
    m_field_size()
{
}

//______________________________________________________________________________
FieldMap::~FieldMap()
{
  ClearField();
}

//______________________________________________________________________________
G4bool
FieldMap::Initialize()
{
  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    G4cerr << FUNC_NAME << G4endl
	   << " file open fail : " << m_file_name << G4endl;
    return false;
  }

  if(m_is_ready){
    G4cerr << FUNC_NAME << G4endl
	   << " * already initialied" << G4endl;
    return false;
  }

  ClearField();

  G4cerr << "   open : " << m_file_name << G4endl;

  if(!(ifs >> m_nx >> m_ny >> m_nz >> m_xmin >> m_ymin >> m_zmin >> m_dx >> m_dy >> m_dz)){
    G4cerr << "#E " << FUNC_NAME << " invalid format" << G4endl;
    return false;
  }

  if(m_nx<0 || m_ny<0 || m_nz<0){
    G4cerr << "#E " << FUNC_NAME << " Nx, Ny, Nz must be positive" << G4endl;
    return false;
  }

  m_xmax = m_xmin + (m_nx - 1) * m_dx;
  m_ymax = m_ymin + (m_ny - 1) * m_dy;
  m_zmax = m_zmin + (m_nz - 1) * m_dz;
  m_field_size.set(m_xmax - m_xmin,
                   m_ymax - m_ymin,
                   m_zmax - m_zmin);
  m_field_size *= CLHEP::cm;

  G4cout << "   x = (" << m_xmin << ", " << m_xmax << ")"
	 << "  y = (" << m_ymin << ", " << m_ymax << ")"
	 << "  z = (" << m_zmin << ", " << m_zmax << ")" << std::endl;

  m_b.resize(m_nx);
  for(G4int ix=0; ix<m_nx; ++ix){
    m_b[ix].resize(m_ny);
    for(G4int iy=0; iy<m_ny; ++iy){
      m_b[ix][iy].resize(m_nz);
    }
  }

  if(m_value_calc == 0. || !std::isfinite(m_value_calc) ||
     m_value_nmr == 0. || !std::isfinite(m_value_nmr)){
    return true;
  }
  const G4double factor = m_value_nmr / m_value_calc;
  G4cout << "   Calc = " << m_value_calc
	 << ", NMR = " << m_value_nmr
	 << ", factor = " << factor << G4endl;

  G4double x, y, z, bx, by, bz;

  G4cout << "   reading fieldmap " << std::flush;

  G4int line = 0;
  while(ifs.good()){
    if(line++%1000000==0)
      G4cout << "." << std::flush;
    ifs >> x >> y >> z >> bx >> by >> bz;
    G4int ix = G4int((x-m_xmin+0.1*m_dx)/m_dx);
    G4int iy = G4int((y-m_ymin+0.1*m_dy)/m_dy);
    G4int iz = G4int((z-m_zmin+0.1*m_dz)/m_dz);
    if(ix>=0 && ix<m_nx && iy>=0 && iy<m_ny && iz>=0 && iz<m_nz){
      m_b[ix][iy][iz].set(bx*factor, by*factor, bz*factor);
    }
  }

  G4cout << " done" << G4endl;
  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
G4bool
FieldMap::GetFieldValue(const G4double point[3] /* cm */,
                        G4double *bfield /* tesla */) const
{
  G4double xt = point[0];
  G4double yt = point[1];
  G4double zt = point[2];

  G4int ix1, ix2, iy1, iy2, iz1, iz2;
  ix1 = G4int((xt-m_xmin)/m_dx);
  iy1 = G4int((yt-m_ymin)/m_dy);
  iz1 = G4int((zt-m_zmin)/m_dz);

  G4double wx1, wx2, wy1, wy2, wz1, wz2;
  if(ix1<0) return false; // { ix1=ix2=0; wx1=1.; wx2=0.; }
  else if(ix1>=m_nx-1) return false; // { ix1=ix2=m_nx-1; wx1=1.; wx2=0.; }
  else { ix2=ix1+1; wx1=(m_xmin+m_dx*ix2-xt)/m_dx; wx2=1.-wx1; }

  if(iy1<0) return false; // { iy1=iy2=0; wy1=1.; wy2=0.; }
  else if(iy1>=m_ny-1) return false; // { iy1=iy2=m_ny-1; wy1=1.; wy2=0.; }
  else { iy2=iy1+1; wy1=(m_ymin+m_dy*iy2-yt)/m_dy; wy2=1.-wy1; }

  if(iz1<0) return false; // { iz1=iz2=0; wz1=1.; wz2=0.; }
  else if(iz1>=m_nz-1) return false; // { iz1=iz2=m_nz-1; wz1=1.; wz2=0.; }
  else { iz2=iz1+1; wz1=(m_zmin+m_dz*iz2-zt)/m_dz; wz2=1.-wz1; }

  G4double bx1 = wx1*wy1*m_b[ix1][iy1][iz1].x() + wx1*wy2*m_b[ix1][iy2][iz1].x()
    + wx2*wy1*m_b[ix2][iy1][iz1].x() + wx2*wy2*m_b[ix2][iy2][iz1].x();
  G4double bx2 = wx1*wy1*m_b[ix1][iy1][iz2].x() + wx1*wy2*m_b[ix1][iy2][iz2].x()
    + wx2*wy1*m_b[ix2][iy1][iz2].x() + wx2*wy2*m_b[ix2][iy2][iz2].x();
  G4double bx  = wz1*bx1 + wz2*bx2;

  G4double by1 = wx1*wy1*m_b[ix1][iy1][iz1].y() + wx1*wy2*m_b[ix1][iy2][iz1].y()
    + wx2*wy1*m_b[ix2][iy1][iz1].y() + wx2*wy2*m_b[ix2][iy2][iz1].y();
  G4double by2 = wx1*wy1*m_b[ix1][iy1][iz2].y() + wx1*wy2*m_b[ix1][iy2][iz2].y()
    + wx2*wy1*m_b[ix2][iy1][iz2].y() + wx2*wy2*m_b[ix2][iy2][iz2].y();
  G4double by  = wz1*by1 + wz2*by2;

  G4double bz1 = wx1*wy1*m_b[ix1][iy1][iz1].z() + wx1*wy2*m_b[ix1][iy2][iz1].z()
    + wx2*wy1*m_b[ix2][iy1][iz1].z() + wx2*wy2*m_b[ix2][iy2][iz1].z();
  G4double bz2 = wx1*wy1*m_b[ix1][iy1][iz2].z() + wx1*wy2*m_b[ix1][iy2][iz2].z()
    + wx2*wy1*m_b[ix2][iy1][iz2].z() + wx2*wy2*m_b[ix2][iy2][iz2].z();
  G4double bz  = wz1*bz1 + wz2*bz2;

  bfield[0] += bx*CLHEP::tesla;
  bfield[1] += by*CLHEP::tesla;
  bfield[2] += bz*CLHEP::tesla;

  return true;
}

//______________________________________________________________________________
G4bool
FieldMap::IsInsideField(G4double* pos) const
{
  return (m_xmin < pos[0] && pos[0] < m_xmax &&
          m_ymin < pos[1] && pos[1] < m_ymax &&
          m_zmin < pos[2] && pos[2] < m_zmax);
}

//______________________________________________________________________________
G4bool
FieldMap::IsInsideField(const G4ThreeVector& pos) const
{
  return (m_xmin < pos.x() && pos.x() < m_xmax &&
          m_ymin < pos.y() && pos.y() < m_ymax &&
          m_zmin < pos.z() && pos.z() < m_zmax);
}

//______________________________________________________________________________
void
FieldMap::ClearField()
{
  for(G4int ix=0; ix<m_nx; ++ix){
    for(G4int iy=0; iy<m_ny; ++iy){
      m_b[ix][iy].clear();
    }
    m_b[ix].clear();
  }
  m_b.clear();
}
