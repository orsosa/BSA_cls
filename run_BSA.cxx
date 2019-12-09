#include "BSA_survey_cls.h"

int main(int argc, char *argv[]){

  if (argc<2){
    std::cout<<"You need to supply a file or a regexp"<<std::endl;
    return 1;
  }
  TString rg="";
  if (argc==3){
    rg = "rgb";
  }
  
  
  BSA_survey_cls t(argv[1],"binning_info.txt",rg);
  t.Loop();
  std::cout<<"end"<<std::endl;
  return 0;
}

