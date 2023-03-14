#include <string>
#include <ctime>
#include <iostream>

using namespace std;

string get_date(){
  time_t now = time(0);
  tm *local_tm = localtime(&now);

  string year = to_string(local_tm->tm_year + 1900);
  string month = to_string(local_tm->tm_mon + 1);
  string day = to_string(local_tm->tm_mday);
  return year+"_"+month+"_"+day;
}

const TVector3 kNuDir(const caf::SRSliceProxy* slc,bool use_true=false){
  TVector3 nu_direction;
  if (use_true){
    if (isnan(slc->truth.position.x)) return TVector3(0.,0.,1.); //assume nu is in z direction
    nu_direction.SetXYZ(slc->truth.position.x-kPrismCentroid[0],
                          slc->truth.position.y-kPrismCentroid[1],
                          slc->truth.position.z+kDistanceFromBNB);
  }
  else{
    if (isnan(slc->vertex.x)) return TVector3(0.,0.,1.); //assume nu is in z direction
    nu_direction.SetXYZ(slc->vertex.x-kPrismCentroid[0],
                          slc->vertex.y-kPrismCentroid[1],
                          slc->vertex.z+kDistanceFromBNB);
  }
  return nu_direction.Unit(); //Direction of neutrino beam
};

//const Var kDoubleToVar([](double w) -> double {return w;}); //Convert double to var
