// -*- C++ -*-

#ifndef LINK_DEF_H
#define LINK_DEF_H

#include <vector>

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class std::vector<std::vector<Int_t>>+;
#pragma link C++ class std::vector<std::vector<Bool_t>>+;
#pragma link C++ class std::vector<std::vector<Char_t>>+;
#pragma link C++ class std::vector<std::vector<Short_t>>+;
#pragma link C++ class std::vector<std::vector<Long_t>>+;
#pragma link C++ class std::vector<std::vector<UChar_t>>+;
#pragma link C++ class std::vector<std::vector<UShort_t>>+;
#pragma link C++ class std::vector<std::vector<UInt_t>>+;
#pragma link C++ class std::vector<std::vector<ULong_t>>+;
#pragma link C++ class std::vector<std::vector<Double_t>>+;
#pragma link C++ class std::vector<std::vector<Char_t*>>+;
#pragma link C++ class std::vector<std::vector<const Char_t*>>+;
#pragma link C++ class std::vector<std::vector<Long64_t>>+;
#pragma link C++ class std::vector<std::vector<ULong64_t>>+;
#pragma link C++ class std::vector<std::vector<void*>>+;
#pragma link C++ class std::vector<TParticle>+;
#pragma link C++ class std::vector<TVector3>+;
#endif

#endif
