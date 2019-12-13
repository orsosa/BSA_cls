#include "BSA_survey_cls.h"

int main(int argc, char *argv[]){

  if (argc<2){
    std::cout<<"You need to supply a file or a regexp"<<std::endl;
    return 1;
  }
  TString opt="";
  if (argc>2){
    for (int k = 2;k<argc;k++){
      opt += " ";
      opt += argv[k];
    }
  }

  TString binfo = "binning_info.txt";
  if (opt.Contains("binfo:")){
    Int_t ind = opt.Index("binfo:");
    opt.Tokenize(binfo,ind);
    binfo = binfo.ReplaceAll("binfo:","").Strip();
  }

  std::cout<<"options: "<<opt<<std::endl;
  std::cout<<"binning info file: "<<binfo<<std::endl;
    
  
  BSA_survey_cls t(argv[1],binfo,opt);
  t.Loop();
  std::cout<<"end"<<std::endl;
  return 0;
}

