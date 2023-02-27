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
