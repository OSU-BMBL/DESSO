#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "OptionParser.h"

int main(int argc, char* argv[]){

  OptionParser option;
  option.add_option("-i","input file",true,"");
  option.add_option("-v","Verbose mode",false,"true");
  option.add_option("-fasta","Input file is fasta format",false,"true");
  option.add_option("-l","Input file is consecutive lines",false,"false");
  option.add_option("-s","Give the number of seed",true,"0");
  option.add_option("-d","To test the double",true,"2.5");
  option.print_descrip();

  std::string filename;
  bool verbose,fasta;
  int seed;
  double d;

  
  std::vector<std::string> argument;
  argument.clear();
  for (int i=1;i<argc;i++){
    std::string t(argv[i]);
    argument.push_back(t);
  }

  std::cout << "Number of arguments: "<<argument.size()<<std::endl;

  option.parse(argument);

  option.get_option("-i",filename);
  option.get_option("-v",verbose);
  option.get_option("-s",seed);
  option.get_option("-d",d);



  std::cout << "Filename: " << filename <<std::endl;
  std::cout << "Verbose: " << verbose << std::endl;
  std::cout << "Seed: " << seed << std::endl; 
  std::cout << "Double: " << d << std::endl;





}

