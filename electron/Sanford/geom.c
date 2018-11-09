//---------------------------------------------------------------------------//
//!
//! \file   brem_dosimetry.c
//! \author Luke Kersting
//! \brief  Geometry for 0.75 MeV Electron Bremsstrahlung dosimetry validation problem
//!
//---------------------------------------------------------------------------//

void brem_dosimetry()
{
  // Set up manager of the geometry world
  gSystem->Load( "libGeom" );
  TGeoManager* geom = new TGeoManager(
    "brem_dosimetry",
    "Geometry for the for 0.75 MeV Electron Bremsstrahlung dosimetry test prob.");

//---------------------------------------------------------------------------//
// Material and Mixture Definitions
//---------------------------------------------------------------------------//

  // CaF2 ( weight fraction: F = 0.48668, Ca = 0.51332, density = 3.180 g/cm^3 )
  TGeoMixture* mat_1 = new TGeoMixture("mat_1", 2, 3.180 );
    mat_1->SetUniqueID( 1 );
    mat_1->DefineElement(0, 40.078, 20, 0.51332);
    mat_1->DefineElement(1, 18.998403163, 9, 0.48668);

  // Al ( density = 2.69890 g/cm^3 )
  TGeoMaterial* mat_2 = new TGeoMaterial( "mat_2", 26.9815385, 13, 2.69890 );

  // C ( density = 1.76 g/cm^3 )
  TGeoMaterial* mat_3 = new TGeoMaterial( "mat_3", 12.011, 13, 2.69890 );

  // Dry Air ( weight fraction: C = 0.000124, N = 0.755267, O = 0.231781, Ar = 0.012827
  // density = 9.987075e-4 g/cm^3 )
  TGeoMixture* mat_4 = new TGeoMixture("mat_4", 4, 9.987075e-4 );
    mat_4->SetUniqueID( 4 );
    mat_4->DefineElement(0, 12.011, 6, 0.000124);
    mat_4->DefineElement(1, 14.007, 7, 0.755267);
    mat_4->DefineElement(2, 15.999, 8, 0.231781);
    mat_4->DefineElement(2, 39.948, 18, 0.012827);

  // Be ( density = 1.84800 g/cm^3 )
  TGeoMaterial* mat_5 = new TGeoMaterial( "mat_5", 9.0121831, 4, 1.848 );

  // Void material
  TGeoMaterial* void_mat = new TGeoMaterial( "void", 0, 0, 0 );

  // Graveyard (terminal)
  TGeoMaterial* graveyard_mat = new TGeoMaterial( "graveyard", 0, 0, 0 );

//---------------------------------------------------------------------------//
// Medium Definitions
//---------------------------------------------------------------------------//

  // CaF2 Dosimeter
  TGeoMedium* med_1 = new TGeoMedium( "med_1", 1, mat_1 );

  // Al encapsulator
  TGeoMedium* med_2 = new TGeoMedium( "med_2", 2, mat_2 );

  // Carbon convertor
  TGeoMedium* med_3 = new TGeoMedium( "med_3", 3, mat_3 );

  // Dry Air (Albuquerque)
  TGeoMedium* med_4 = new TGeoMedium( "med_4", 4, mat_4 );

  // Be encasing
  TGeoMedium* med_5 = new TGeoMedium( "med_5", 5, mat_5 );

  // Void
  TGeoMedium* void_med = new TGeoMedium( "void_med", 6, void_mat );

  // Graveyard
  TGeoMedium* graveyard_med = new TGeoMedium( "graveyard", 7, graveyard_mat );

//---------------------------------------------------------------------------//
// Volume Definitions
//---------------------------------------------------------------------------//

  // Create the CaF2 cylinder (R = 60.0 cm, h = 0.089 cm)
  double half_height_1 = 0.0445; // h/2
  TGeoVolume* region_1 =
    geom->MakeTube( "REGION1", med_1, 0.0, 60.0, half_height_1 );

  region_1->SetUniqueID( 1 );

  // Create the Al encapsulators (R = 60.0 cm, h = 0.22 cm)
  double half_height_2 = 0.11;
  TGeoVolume* region_2 =
    geom->MakeTube( "REGION2", med_2, 0.0, 60.0, half_height_2 );

  region_2->SetUniqueID( 2 );

  // Create the C convertor (R = 1.0 cm, h = 0.48 cm)
  double half_height_3 = 0.24
  TGeoVolume* region_3 =
    geom->MakeTube( "REGION3", med_3, 0.0, 1.0, half_height_3 );

  region_3->SetUniqueID( 3 );

  // Create the source void (R = 1.0 cm, h = 0.005 cm)
  double half_height_4 = 0.0025;
  TGeoVolume* region_4 =
    geom->MakeTube( "REGION4", void_med, 0.0, 1.0, half_height_4 );

  region_4->SetUniqueID( 4 );

  // Create the dry air (R = 60.0 cm, h = 30.539 cm)
  double half_height_5 = 15.2695;
  TGeoVolume* region_5 =
    geom->MakeTube( "REGION5", med_4, 0.0, 60.0, half_height_5 );

  region_5->SetUniqueID( 5 );

  // Create the Be encasing (R = 65.0 cm, h = 35.1 cm)
  double half_height_6 = 17.55;
  TGeoVolume* region_6 =
    geom->MakeTube( "REGION6", med_5, 0.0, 65.0, half_height_6 );

  region_5->SetUniqueID( 6 );

  // Create the graveyard volume (R = 70.0 cm, h = 72.0 cm)
  double half_height_graveyard = 36.0;
  TGeoVolume* graveyard_region =
    geom->MakeTube( "GRAVEYARD", graveyard_med, 0.0, 70.0, 36.0 );

  graveyard_region->SetUniqueID( 7 );

//---------------------------------------------------------------------------//
// Heirarchy (Volume) Definitions
//---------------------------------------------------------------------------//

geom->SetTopVolume(ENVL);
 ENVL->AddNode(PRER,1,new TGeoTranslation(0,0,1496.712));
  PRER->AddNode(EHOL,1,gGeoIdentity);
 ENVL->AddNode(CHOL,1,new TGeoTranslation(0,0,1582.346));
 ENVL->AddNode(CAL0,1,new TGeoTranslation(0,-132.36,1574.668));
  TGeoVolume *CAL1 = CAL0->Divide("CAL1",1,44,-242.66,11.03);
   TGeoVolume *CAL2 = CAL1->Divide("CAL2",2,20,-110.3,11.03);
    CAL2->AddNode(CLW1,1,gGeoIdentity);

  // Set the graveyard to be the top volume (rest-of-universe)
  geom->SetTopVolume( graveyard_volume );

  // Place the Be encasing inside of the graveyard (z_min = -0.1, z_max = 35.0)
  double z_trans = half_height_6 - 0.1;
  graveyard_volume->AddNode( region_6, 1, new TGeoTranslation(0, 0, z_trans) );

    // Place the dry air inside of the Be encasing (z_min = -0.005, z_max = 30.534)
    region_6->AddNode( region_5, 1 );

      // Place the dry air inside of the Be encasing
      region_6->AddNode( region_5, 1 );


  region_1->AddNode( )

  // Place the hydrogen sphere 1 inside of hydrogen sphere 2
  h_sphere_volume2->AddNode( h_sphere_volume1, 1 );

  // Place the hydrogen sphere 2 inside of hydrogen sphere 3
  h_sphere_volume3->AddNode( h_sphere_volume2, 1 );

  // Place the hydrogen sphere 3 inside of hydrogen sphere 4
  h_sphere_volume4->AddNode( h_sphere_volume3, 1 );

  // Place the hydrogen sphere 4 inside of hydrogen sphere 5
  h_sphere_volume5->AddNode( h_sphere_volume4, 1 );

  // Place the hydrogen sphere 5 inside of the void sphere
  void_sphere_volume->AddNode( h_sphere_volume5, 1 );

  // Place the void sphere inside of the graveyard
  graveyard_volume->AddNode( void_sphere_volume, 1 );

  // Set the graveyard to be the top volume (rest-of-universe)
  geom->SetTopVolume( graveyard_volume );

//---------------------------------------------------------------------------//
// Export and Drawing Capabilities
//---------------------------------------------------------------------------//

  // Close the geometry
  geom->SetTopVisible();
  geom->CloseGeometry();

  // Draw the geometry
  // h_sphere_volume5->Draw();

  // Export the geometry
  geom->Export( "brem_dosimetry.root" );

  // Finished
  exit(0);
}

//---------------------------------------------------------------------------//
// end h_sphere_1kev.c
//---------------------------------------------------------------------------//
