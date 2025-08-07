// -*- C++ -*-

#ifndef FUNK_NAME_HH
#define FUNK_NAME_HH

#include <G4String.hh>

#define FUNC_NAME G4String("["+G4String(ClassName())+"::"+__func__+"()]")

#endif
