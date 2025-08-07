// -*- C++ -*-

#include "PhysicsList.hh"

#include <iomanip>
#include <CLHEP/Units/SystemOfUnits.h>

#include <globals.hh>
#include <G4ios.hh>

#include <G4DecayPhysics.hh>
#include <G4EmStandardPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4IonPhysics.hh>
#include <G4StoppingPhysics.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4NeutronTrackingCut.hh>
#include <G4HadronPhysicsQGSP_BERT.hh>
#include <G4VPhysicsConstructor.hh>

// ionization process and model of generic ions 
#include <G4ionIonisation.hh>
#include <G4BraggIonGasModel.hh>
#include <G4BetheBlochIonGasModel.hh>
#include <G4IonFluctuations.hh>
#include <G4UniversalFluctuation.hh>
#include <G4EmParameters.hh>
// Multipple scattering process and model of generic ions
#include <G4hMultipleScattering.hh>
#include <G4UrbanMscModel.hh>

// Single coulomb scattering process and model of generic ions
#include <G4CoulombScattering.hh>
#include <G4IonCoulombScatteringModel.hh>

// Nuclear stopping process and model of generic ions
#include <G4NuclearStopping.hh>
#include <G4ICRU49NuclearStoppingModel.hh>

#include "ConfMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"

#define G4MT_physicsVector                                                    \
  ((G4VMPLsubInstanceManager.offset[g4vmplInstanceID]).physicsVector)

namespace
{
  const auto& confMan = ConfMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
PhysicsList::PhysicsList(G4int ver)
  : G4VUserPhysicsList()
{
  verboseLevel = ver;
  defaultCutValue = 0.7*CLHEP::mm;
  // confMan.Get<G4double>("DefaultCutValue");

  // EM Physics
  RegisterPhysics(new G4EmStandardPhysics(ver));
  // Synchrotron Radiation & GN Physics
  RegisterPhysics(new G4EmExtraPhysics(ver));

  // Decays
  RegisterPhysics(new G4DecayPhysics(ver));

  // Hadron Elastic scattering
  RegisterPhysics(new G4HadronElasticPhysics(ver));
  // Hadron Physics
  RegisterPhysics(new G4HadronPhysicsQGSP_BERT(ver));

  ///// Others
  // Stopping Physics
  RegisterPhysics(new G4StoppingPhysics(ver));
  // Ion Physics
  RegisterPhysics(new G4IonPhysics(ver));
  // Neutron tracking cut
  RegisterPhysics(new G4NeutronTrackingCut(ver));
}

//_____________________________________________________________________________
void
PhysicsList::ConstructParticle()
{
  for(auto itr = G4MT_physicsVector->cbegin();
      itr != G4MT_physicsVector->cend(); ++itr)
  {
    (*itr)->ConstructParticle();
  }
}

//_____________________________________________________________________________
void
PhysicsList::ConstructProcess()
{

  const G4int pad_configure = gSize.Get("TpcPadConfigure");

  
  // G4AutoLock l(&constructProcessMutex);
  AddTransportation();

  if(pad_configure ==4)
    AddIonGasProcess(); 


  for(auto itr = G4MT_physicsVector->cbegin();
      itr != G4MT_physicsVector->cend(); ++itr)
  {
    auto name = (*itr)->GetPhysicsName();
    if(name == "G4EmStandard" || name == "G4GammaLeptoNuclearPhys"){
      if(!confMan.Get<G4bool>("EM")) continue;
    }else if(name == "Decay"){
      if(!confMan.Get<G4bool>("DECAY")) continue;
    }else{
      if(!confMan.Get<G4bool>("HADRON")) continue;
    }

    if(verboseLevel > 0)
      G4cout << FUNC_NAME << " Construct " << name << G4endl;

    (*itr)->ConstructProcess();


    if(pad_configure ==4){
      G4EmParameters* emParameters = G4EmParameters::Instance();
      emParameters->SetMinEnergy(10*CLHEP::eV);
      emParameters->SetMaxEnergy(2.*CLHEP::GeV);
      emParameters->SetNumberOfBinsPerDecade(100);
    }

  }
}

//_____________________________________________________________________________
void PhysicsList::AddIonGasProcess()
{
    auto ph = G4PhysicsListHelper::GetPhysicsListHelper();
    auto pIterator = GetParticleIterator();
    pIterator->reset();
    while((*pIterator)())
    {
        G4ParticleDefinition *pDefinition = pIterator->value();
        G4String pName = pDefinition->GetParticleName();
        if(pName == "proton" || pName == "kaon" || pName == "pion" || pName == "muon")
        {
            // effective charge and energy loss model of ion
            G4ionIonisation *iIon = new G4ionIonisation();
            G4BraggIonGasModel *bIgm = new G4BraggIonGasModel();
            G4BetheBlochIonGasModel *bbIgm = new G4BetheBlochIonGasModel();
	    
	    bIgm->SetActivationHighEnergyLimit(2.*CLHEP::MeV*pDefinition->GetPDGMass()/CLHEP::proton_mass_c2);
	    bbIgm->SetActivationLowEnergyLimit(2.*CLHEP::MeV*pDefinition->GetPDGMass()/CLHEP::proton_mass_c2);

            iIon->AddEmModel(0, bIgm, new G4IonFluctuations);
            iIon->AddEmModel(0, bbIgm, new G4UniversalFluctuation);
	    
            // no delta ray
            iIon->ActivateSecondaryBiasing("World", 1e-10, 100*CLHEP::TeV);
	    
            G4hMultipleScattering *hMsc = new G4hMultipleScattering();
            hMsc->AddEmModel(0, new G4UrbanMscModel());
            G4CoulombScattering *csc = new G4CoulombScattering();
            csc->AddEmModel(0, new G4IonCoulombScatteringModel());
            G4NuclearStopping *nsp = new G4NuclearStopping();
            nsp->AddEmModel(0, new G4ICRU49NuclearStoppingModel());
            ph->RegisterProcess(iIon, pDefinition);
            ph->RegisterProcess(hMsc, pDefinition);
            ph->RegisterProcess(csc, pDefinition);
            ph->RegisterProcess(nsp, pDefinition);
        }
    }
}

//_____________________________________________________________________________
void
PhysicsList::SetCuts()
{
  SetCutsWithDefault();
}
