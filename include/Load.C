void Load(){
  
  TString HIPO=gSystem->Getenv("HIPO");
  gSystem->Load(HIPO+"/lib/libhipo4");
  gInterpreter->AddIncludePath(HIPO+"/include");


}
