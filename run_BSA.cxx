#include "BSA_survey_cls.h"


int main(int argc, char *argv[]){

  if (argc!=2){
    std::cout<<"You need to supply a file or a regexp"<<std::endl;
    return 1;
  }
  
  BSA_survey_cls t(argv[1],"binning_info.txt");
  t.Loop();
  return 0;
}
