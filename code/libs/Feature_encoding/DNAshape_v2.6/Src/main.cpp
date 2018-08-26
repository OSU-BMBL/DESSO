#include "common.h"
#include "utilities.h"
#include "OptionParser.h"
#include "TableCompiler.h"
#include "Output.h"
#include "ReadData.h"
#include "prediction.h"
#include "SNP.h"

int main(int argc, char* argv[]){

  int initial_time = get_time();
  std::string exec_path,lib_list,inosine_list;     
//lib_list     file containing the list of location of *.minor files, relative path

  exec_path = get_exec_path(argv[0]);
  //lib_list=exec_path+"/SimData/poollist9";
  lib_list = exec_path+"/SimData/nonartifacts.txt";
  inosine_list = exec_path+"/SimData/Inosine.txt";
  std::string querytable_filename = exec_path + "/QueryTable.dat";
  //lib_list = exec_path + "/SimData/dummy.txt";

  OptionParser option;
  option.add_option("-i","Input file containing the DNA sequence.",true,"");
  //  option.add_option("-o","Output file",true,"");
  option.add_option("-v","Verbose mode. Output all useful information.",false,"false");
  option.add_option("-d","Debug mode. Output additional information.",false,"false");
  option.add_option("-sdist","Generate distribution of step parameters (slide,roll,twist) for all 10 dinucleotides.",false,"false");
  option.add_option("-width","Specify the maximum number of values each line of output contain. Default value is 30.",true,"30");
  option.add_option("-SNP","Generate the SNP table.",false,"false");
  option.add_option("-five","Define the minor groove width of a base pair as the average minor groove width of the 5 levels centered around that base pair instead of 3 levels by default.",false,"false");
  option.add_option("-r","Recompile the query table from the simulation data.",false,"false");
  option.add_option("-inosine","Incorporate the inosine simulation data.",false,"false");
  option.add_option("-delimiter","Specify which delimiter to use in formatting the output.(1-comma, 2-semicolon, 3-space)",true,"1");
  option.add_option("-lib","Specify the the lib file.",true,"");
  
  std::string ifilename,libfile;
  int output_width,delimiter_type;
  char delimiter;
  //std::string ofilename;
  bool verbose,debug,sdist,SNP,five,recompile,inosine;
  libfile = "";
  
  if (argc==2){       //  -h     output the description information
    std::string first_arg(argv[1]);
    if (first_arg=="-h")      
      option.print_descrip();
  }  
  if (argc==1){
    option.print_descrip();
  }
  
  std::vector <std::string> argument;
  argument.clear();
  for (int i=1;i<argc;i++){
    std::string t(argv[i]);
    argument.push_back(t);
  }

  option.parse(argument);
  option.get_option("-i", ifilename);
  //option.get_option("-o",ofilename);
  option.get_option("-v",verbose);
  option.get_option("-d",debug);
  option.get_option("-sdist",sdist);
  option.get_option("-width",output_width);
  option.get_option("-SNP",SNP);
  option.get_option("-five",five);
  option.get_option("-r",recompile);
  option.get_option("-inosine",inosine);
  option.get_option("-delimiter",delimiter_type);
  option.get_option("-lib",libfile);

  if (!libfile.empty()) lib_list = libfile;

  if (delimiter_type == 2)
    delimiter = ';';
  else if (delimiter_type == 3)
    delimiter = ' ';
  else
    delimiter = ',';

  //Initialize the map
  DNA_to_properties pentamers_map,inosine_map; 
  pentamers_map.clear();
  //Inosine
  inosine_map.clear();
  if (inosine){
    std::ifstream in_inosine_list(inosine_list.c_str());
    if (!in_inosine_list){
      std::cerr << "Cannot open the inosine_list file: "<<inosine_list<<"\n";
      exit(1);
    }
    std::string line;
    while (getline(in_inosine_list,line)){
	if (!line.empty()){
	  std::string minor_file = exec_path + "/SimData" + line.substr(1,line.size());
	  std::string cclis_file = minor_file;
	  cclis_file.replace(minor_file.size()-6,6,"_cc.lis");
	  std::string basename = minor_file.substr(0,minor_file.find(".minor"));   // without extension	
	  std::string seq_file = basename+".seq";
	  std::string sequence;
	  std::ifstream in_seq_file(seq_file.c_str());
	  if (in_seq_file){    // file exist
	    getline(in_seq_file,sequence);
	    remove_terminal_space(sequence);	  
	    in_seq_file.close();

	    std::vector<double_vector> matrix;	  

	    if ((!sequence.empty()) and (sequence.size()==13)){     //start to input the data to matrix
	      //Information about the Matrix:
	      //  #0 Pos   #1 Level   #3 Average of Minor Groove Width
	      add_groove_width_to_inosine_table(minor_file,inosine_map,sequence,3,"minor",verbose,debug,five);	    	  
	    }	    

	  } // if (in_seq_file)
	  in_seq_file.close();
	}
    }
    DNA_to_properties::const_iterator end=inosine_map.end();
    for (DNA_to_properties::iterator it=inosine_map.begin();it!=end;++it)
      it->second.calc_ave_sd(debug);

  }

  if ((recompile) or (!libfile.empty())){   //Compile the query table from the simulation table or Perform cross validation
    std::ifstream in_lib_list(lib_list.c_str());
    if (!in_lib_list){     //check if the file exist
      std::cerr << "Cannot open the lib_list file: " << lib_list <<"\n";
      exit(1);
    }

    //Generate all 512 unique pentamers
    build_unique_pentamers(pentamers_map);
    //if (debug)
    //  output_pentamers_map(pentamers_map,"minor");

    if (!verbose)
      std::cout << "Compiling Query Table......"<<std::endl;

    std::string line;
    int no_of_lib_files=0;
    while (getline(in_lib_list,line)){
      if (!line.empty()){
	std::string minor_file = exec_path + "/SimData/" + line.substr(1,line.size());
	std::string major_file = minor_file;
	std::string cclis_file = minor_file;
	major_file.replace(minor_file.size()-5,5,"major");
	cclis_file.replace(minor_file.size()-6,6,"_cc.lis");
	std::string basename = minor_file.substr(0,minor_file.find(".minor"));   // without extension	
	std::string seq_file = basename+".seq";
	std::string sequence;
	if (debug){
	  std::cout << "Minor file: " << minor_file << std::endl;
	  std::cout << "Major file: " << major_file << std::endl;
	  std::cout << "seq file: "<< seq_file << std::endl;
	  std::cout << "cclis_file: "<< cclis_file << std::endl;
	}	  
	std::ifstream in_seq_file(seq_file.c_str());
	if (in_seq_file){    // file exist
	  if (debug)
	    std::cout << "File exist\n";

	  getline(in_seq_file,sequence);
	  remove_terminal_space(sequence);	  
	  in_seq_file.close();

	  if (debug){
	    std::cout << "Sequence: "<< sequence << std::endl;
	  }

	  std::vector<double_vector> matrix;	  

	  if ((!sequence.empty()) and (is_DNA(sequence))){     //start to input the data to matrix
	    no_of_lib_files++;	   

	    //Information about the Matrix:
	    //  #0 Pos   #1 Level   #3 Average of Minor Groove Width
	    add_groove_width_to_pentamers_table(minor_file,pentamers_map,sequence,3,"minor",verbose,debug,five);   
	    add_groove_width_to_pentamers_table(major_file,pentamers_map,sequence,3,"major",verbose,debug,five);	    
	    add_step_info_to_pentamers_table(pentamers_map,cclis_file,sequence,verbose,debug);
	  }
	  else if (verbose)
	    std::cout<<"\t Cannot recognize the sequence \n";

	} // if (in_seq_file)
	in_seq_file.close();
      }
    }
    in_lib_list.close();

    if (verbose)
      std::cout << "# of valid lib files: " << no_of_lib_files << std::endl;
    else
      std::cout <<"Query Table Compiled...... "<<std::endl;

    //Calculate the average and standard deviation for each pentamer
    if (debug)
      std::cout << "Call calc_ave_ad() \n";
    DNA_to_properties::const_iterator end=pentamers_map.end();
    for (DNA_to_properties::iterator it=pentamers_map.begin();it!=end;++it)
      it->second.calc_ave_sd(debug);
    if (debug)
      std::cout << "Finished \n";
  }
  else{ // load the query table from file
    process_querytable_file(querytable_filename,pentamers_map);
  }

  if (verbose){
    string_vector output_object_list;
    std::cout << std::endl << "Groove Width:" <<std::endl;
    output_object_list.clear();
    output_object_list.push_back("minor");
    output_object_list.push_back("major");
    output_object_list.push_back("propel");
    output_pentamers_map(pentamers_map,output_object_list,inosine);

    std::cout << std::endl << "Step Information:" <<std::endl;
    output_object_list.clear();
    output_object_list.push_back("slide1");
    output_object_list.push_back("slide2");
    output_object_list.push_back("roll1");
    output_object_list.push_back("roll2");
    output_object_list.push_back("twist1");
    output_object_list.push_back("twist2");
    output_pentamers_map(pentamers_map,output_object_list,inosine);

    output_object_list.clear();
    output_object_list.push_back("minor");
    output_object_list.push_back("major");
    output_object_list.push_back("propel");
    if (inosine){
      std::cout << std::endl << "Inosine:"<<std::endl;
      output_pentamers_map(inosine_map,output_object_list,inosine);
    }
  }

  if (recompile){
    output_pentamers_map_to_querytable(pentamers_map,querytable_filename);
  }  

  if (sdist){
    std::cout << "Roll:"<<std::endl;
    step_parameters_distribution(pentamers_map,"roll");
    std::cout << "Twist:"<<std::endl;
    step_parameters_distribution(pentamers_map,"twist");
  }

  if (SNP){
    SNP_distribution(pentamers_map);
    exit(0);
  }

  std::ifstream in_fstream(ifilename.c_str());
  if (!in_fstream){
    std::cerr << "Cannot open the input file:  "<< ifilename << "\n";
    exit(1);
  }
  
  std::cout << "Reading the input sequence......"<<std::endl;
  string_vector sequence_list;
  string_vector name_list;
  sequence_list.clear();
  name_list.clear();
  read_fasta(in_fstream,sequence_list,name_list,debug);

  std::string minor_filename = ifilename + ".MGW" ;
  //std::string major_filename = ifilename + ".major" ;
  //std::string slide_filename = ifilename + ".slide" ;
  std::string roll_filename = ifilename + ".Roll" ;
  std::string twist_filename = ifilename + ".HelT";
  std::string propel_filename = ifilename + ".ProT";

  if (inosine){
    predict_groove_width_inosine(minor_filename,sequence_list,name_list,pentamers_map,inosine_map);
    exit(0);
  }

  //convert sequence_list to pointers_list
  std::vector <pointers_vector> pointers_matrix;
  std::vector <int_vector> status_matrix;
  std::cout << "Indexing the input sequence......"<<std::endl;
  convert_sequence_list(sequence_list,pointers_matrix,status_matrix,pentamers_map);
  std::cout << "Indexing complete"<<std::endl;


  std::stringstream current_ss;

  std::cout << "Processing......"<<std::endl;

  predict_groove_width(current_ss,pointers_matrix,status_matrix,name_list,debug,pentamers_map,"minor",output_width,delimiter);  
  output_stringstream_to_file(current_ss,minor_filename);
  std::cout << "\r" << "25% Complete"<<std::flush;

  //current_ss.str("");
  //current_ss.clear();
  //predict_groove_width(current_ss,pointers_matrix,status_matrix,name_list,debug,pentamers_map,"major",output_width,delimiter);
  //output_stringstream_to_file(current_ss,major_filename);	\
  //std::cout << "\r" << "22% Complete"<<std::flush;

  current_ss.str("");
  current_ss.clear();
  predict_groove_width(current_ss,pointers_matrix,status_matrix,name_list,debug,pentamers_map,"propel",output_width,delimiter);
  output_stringstream_to_file(current_ss,propel_filename);
  std::cout << "\r" << "50% Complete"<<std::flush;
  
  //current_ss.str("");
  //current_ss.clear();
  //predict_step_parameters(current_ss,pointers_matrix,status_matrix,name_list,debug,pentamers_map,"slide",output_width,delimiter);
  //output_stringstream_to_file(current_ss,slide_filename);
  //std::cout << "\r" << "55% Complete"<<std::flush;

  current_ss.str("");
  current_ss.clear();
  predict_step_parameters(current_ss,pointers_matrix,status_matrix,name_list,debug,pentamers_map,"roll",output_width,delimiter);
  output_stringstream_to_file(current_ss,roll_filename);
  std::cout << "\r" << "75% Complete"<<std::flush;

  current_ss.str("");
  current_ss.clear();
  predict_step_parameters(current_ss,pointers_matrix,status_matrix,name_list,debug,pentamers_map,"twist",output_width,delimiter);
  output_stringstream_to_file(current_ss,twist_filename);
  std::cout << "\r" << "100% Complete"<<std::flush;

  std::cout << std::endl;
  std::cout << "Total running time: "<< (get_time()-initial_time)<<" seconds"<<std::endl;

}
