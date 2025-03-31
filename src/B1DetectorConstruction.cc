//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VSolid.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
     
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  
G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
G4Box* solidWorld = new G4Box("WorldBox", 100*cm, 100*cm, 100*cm);
G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "LogicalWorld");
G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

//Lead shield
G4Material* matPb = nist->FindOrBuildMaterial("G4_Pb");
G4ThreeVector pos1 = G4ThreeVector(0*cm, 0*cm, 0*cm);
G4Box* PbBox = new G4Box("PbBox", 18*cm, 46*cm, 18*cm);     
G4LogicalVolume* LogPb = new G4LogicalVolume(PbBox, matPb, "LogicmatPb");           
new G4PVPlacement(0, pos1, LogPb, "Pb_layer", logicWorld, false, 0, checkOverlaps);                     
//Strontium source                              
G4Material* matSr = nist->FindOrBuildMaterial("G4_Sr");
G4Box* SrBox = new G4Box("SrBox", 7.5*cm, 4.5*cm, 7.5*cm);     
G4LogicalVolume* LogSr = new G4LogicalVolume(SrBox, matSr, "LogicmatSr");   
for (int j = 1; j < 9; j++){
                      new G4PVPlacement(0, G4ThreeVector(0*cm, -31.5*cm+((j-1)*9)*cm,0*cm), LogSr,"Source", LogPb, false, j, checkOverlaps);         
}
//Aluminum outer shell 
G4Material* matAl = nist->FindOrBuildMaterial("G4_Al");
G4Box* AlLeaf = new G4Box("AlLeaf", 25*cm, 46.0*cm, 1*cm);
G4LogicalVolume* LogAl = new G4LogicalVolume(AlLeaf, matAl, "LogicmatAl");
new G4PVPlacement(0, G4ThreeVector(43*cm, 0*cm, 0*cm), LogAl, "Aluminum_Front", logicWorld, false, 0, checkOverlaps);
new G4PVPlacement(0, G4ThreeVector(-43*cm, 0*cm, 0*cm), LogAl, "Aluminum_Back", logicWorld, false, 1, checkOverlaps);
G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
rotationMatrix->rotateY(90.*deg);
new G4PVPlacement(rotationMatrix, G4ThreeVector(0*cm, 0*cm, 43*cm), LogAl, "Aluminum_Right", logicWorld, false, 2, checkOverlaps);
new G4PVPlacement(rotationMatrix, G4ThreeVector(0*cm, 0*cm, -43*cm), LogAl, "Aluminum_Left", logicWorld, false, 3, checkOverlaps); 
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
