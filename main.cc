// -*- C++ -*-

#include <G4MTRunManager.hh>
#include <G4RunManager.hh>
#include <G4RunManagerFactory.hh>
#include <G4TrajectoryDrawByCharge.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>
#include <G4UIExecutive.hh>
#include <G4VisExecutive.hh>
#include <QGSP_BERT.hh>

#include "ActionInitialization.hh"
#include "ConfMan.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

//#include "VisManager.hh"

//_____________________________________________________________________________
int
main(int argc, char** argv)
{
  std::vector<G4String> arg(argv, argv + argc);

  if (argc != kArgc-1 && argc != kArgc) {
    G4cout << "Usage: " << arg[kProcess]
	   << " [ConfFile] [OutputName] (G4Macro)" << G4endl;
    return EXIT_SUCCESS;
  }

  auto& gConf = ConfMan::GetInstance();
  if (!gConf.Initialize(arg)) {
    return EXIT_FAILURE;
  }

  auto runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  // G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new ActionInitialization);
  runManager->Initialize();

  auto visManager = new G4VisExecutive;
  visManager->SetVerboseLevel(0);
  visManager->Initialize();

  auto uiManager = G4UImanager::GetUIpointer();

  // interactive session, if no arguments given
  if(argc == kArgc-1) {
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    uiManager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  // batch mode
  else if(argc == kArgc) {
    G4String command("/control/execute ");
    G4String fileName(argv[kG4Macro]);
    uiManager->ApplyCommand(command + fileName);
  }

  delete visManager;
  delete runManager;
  return EXIT_SUCCESS;
}
