#include "BSA_survey_cls.h"


int main(){
  BSA_survey_cls t("../JLAB_DATA/data/RGA/mix/pruned_evData_skim4_5204_pippim_*.root");
  t.Loop();
  return 0;
}
