#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

TFile* CheckFile( TString path_name ){

  if (gSystem -> AccessPathName(path_name)) {
    std::cout << "The " << path_name << " file does not exist. Creating it..." << std::endl;
  } else {
    std::cout << "Opening file: " << path_name << std::endl;
  }

  TFile *file = new TFile(path_name, "UPDATE");
  if ( !(file -> IsOpen()) ) { std::cout << "Unable to open file." << std::endl;}

  return file;
}
