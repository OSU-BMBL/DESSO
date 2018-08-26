#include "utilities.h"

int main(){
  std::string upper="abCDEfgH";
  convert_to_uppercase(upper);
  std::cout<<upper<<std::endl;
  std::cout<<std::endl;
  
  string_vector v1;
  parse_string_to_list("  ABCD,BCD;ABC, ;EFG"," ,;",v1);
  for (int i=0;i<v1.size();i++)
    std::cout << i << " " <<v1[i] <<std::endl;
  parse_string_to_list("  ABCD,BCD;ABC, ;EFG  ; "," ,;",v1);
  for (int i=0;i<v1.size();i++)
    std::cout << i << " " <<v1[i] <<std::endl;

  std::cout<<get_exec_path()<<std::endl;

  std::cout<<opposite_strand("ATTGCA")<<std::endl;

  return 0;
}
