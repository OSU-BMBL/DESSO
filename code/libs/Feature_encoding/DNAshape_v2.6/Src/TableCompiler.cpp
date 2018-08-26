#include "TableCompiler.h"

void add_groove_width_to_pentamers_table(std::string filename,DNA_to_properties &onemap,std::string sequence,int object_index, \
					 std::string object_name,bool verbose,bool debug,bool five){

  //construct the matrix from the file
  std::vector<double_vector> matrix;
  matrix.clear();
  std::ifstream in_file(filename.c_str());

  if (!in_file){
    if (debug)
      std::cout << "Cannot open the file: " << filename << std::endl;
    in_file.close();
    return;
  }

  std::string line;

  while (getline(in_file,line)){
    double_vector double_row;
    double_row.clear();
    string_vector string_row;
    parse_string_to_list(line,"\t ",string_row);		
    if ((string_row.size()>0) and (string_vector_to_double_vector(string_row,double_row)))
      matrix.push_back(double_row);
  }
  //if the file is empty, then return
  if (matrix.size()==0)
    return;

  //  Matrix definition
  //  0--ResId    1--Levels  3--Average minor/major groove width
  //  8--occurrences

  int no_of_positions = sequence.size();

  if (debug){
    for (int i=0;i<matrix.size();i++){
      for (int j=0;j<matrix[i].size();j++)
	std::cout <<"\t"<<matrix[i][j];
      std::cout<<"\n";
    }
    std::cout<<"No of positions: "<<no_of_positions<<std::endl;
  }

  std::vector<int> num_of_levels;     //num of levels for each ResId
  for (int i=0;i<=no_of_positions;i++)
    num_of_levels.push_back(0);
  
  //std::cout << "Vector size: "<<num_of_levels.size()<<std::endl;
  for (int i=0; i<matrix.size(); i++){
    //std::cout << double_to_int(matrix[i][0])  <<std::endl;
    num_of_levels[double_to_int(matrix[i][0])]++;
  }

  /*  
  int starting_row_no[MAX_SEQ_LEN_IN_LIB]={0};    //starting row no. in matrix for each position, initialize all the elements to 0
  int no_of_levels[MAX_SEQ_LEN_IN_LIB]={0};
  for (int i=0;i<matrix.size();i++){
    if (no_of_levels[double_to_int(matrix[i][0])]==0)
      starting_row_no[double_to_int(matrix[i][0])]=i;
    no_of_levels[double_to_int(matrix[i][0])]++;
  }
  */
  //std::cout << "Mark\n";
 
  int start_pos = 4;  //discard the first 3 positions      //the ResID of the first sequence character is 1
  int end_pos=no_of_positions-4; // discard the last 3
  double current_value,prev_value,next_value,prevprev_value,nextnext_value;
  int current_id;
  std::string current_pentamer;
  for (int i=1;i<(matrix.size()-1);i++){
    if ((double_to_int(matrix[i][0])>=start_pos) and (double_to_int(matrix[i][0])<=end_pos) \
	and (double_to_int(matrix[i][1])==1) and (double_to_int(matrix[i][8])>100000)){   //the first is level is 1
            
      current_id = double_to_int(matrix[i][0]);
      current_pentamer = sequence.substr(current_id-1-2,5);  //1bp shift
      current_value = matrix[i][object_index];
      
      if ((double_to_int(matrix[i+1][0])==current_id) and (double_to_int(matrix[i+1][1])==2) \
	  and (double_to_int(matrix[i+1][8])>100000))
	next_value = matrix[i+1][object_index];
      else
	continue;
      
      if (five){  // Use 5 levels instead of 3 
	if ((i+2<matrix.size()) and (double_to_int(matrix[i+2][0])==current_id)	\
	    and (double_to_int(matrix[i+2][1])==3) and (double_to_int(matrix[i+2][8])>100000))
	  nextnext_value = matrix[i+2][object_index];
	else 
	  continue;
      }	    
      
      int tp = i;
      bool flag = false;
      for (int j=0;j<num_of_levels[current_id-1];j++){
	tp--;
	if ((double_to_int(matrix[tp][1])<=4) and (double_to_int(matrix[tp][8])>100000)){
	    flag = true;
	    prev_value = matrix[tp][object_index];
	    break;
	  }
      }
      if (!flag) continue;

      if (five){
	if ((tp-1>=0) and (double_to_int(matrix[tp-1][0])==current_id-1) and (double_to_int(matrix[tp-1][8])>100000)){
	  prevprev_value = matrix[tp-1][object_index];
	}
	else
	  continue;
      }

      double final_value;
      if (five)
	final_value = (prevprev_value + prev_value + current_value + next_value + nextnext_value)/5;
      else
	final_value = (prev_value + current_value + next_value)/3;

      if (found_str_in_map(current_pentamer,onemap))
	onemap[current_pentamer].push(final_value,object_name);
      else
	onemap[opposite_strand(current_pentamer)].push(final_value,object_name);
    }
  }
  /*
  for (int i=start_pos;i<=end_pos;i++){    
    double current_value =0;
    std::string current_pentamer = sequence.substr(i-1-2,5);   //there is 1bp shift between pos and sequence index       
    if (((starting_row_no[i]+1)<matrix.size()) and ((starting_row_no[i]-1)>=0) and \
	(matrix[starting_row_no[i]+1][0]==matrix[starting_row_no[i]][0]) and \
	(no_of_levels[i-1]>0) and (matrix[starting_row_no[i]][1]==1) and \
	(matrix[starting_row_no[i]-1][0]==(matrix[starting_row_no[i]][0]-1))) {  //make sure it is consecutive, the first level is level 1
      //next level is the same position, and its level 2      
      current_value = matrix[starting_row_no[i]][object_index]+\
	matrix[starting_row_no[i]+1][object_index];           //currentvalue =   3 |||0,1
      if (matrix[starting_row_no[i]-1][1]==5){
	if ((starting_row_no[i]-2<0) or ((matrix[starting_row_no[i]-2][0]+1)!=matrix[starting_row_no[i]][0]))
	  continue;
	current_value+=matrix[starting_row_no[i]-2][object_index];
      }
      else 
	current_value+=matrix[starting_row_no[i]-1][object_index];
      current_value=current_value/3;    
      if (debug)
	std::cout<< current_pentamer << " : "<<current_value <<std::endl;
      if (found_str_in_map(current_pentamer,onemap))
	onemap[current_pentamer].push(current_value,object_name);
      else
	onemap[opposite_strand(current_pentamer)].push(current_value,object_name);
    }
  }
  */
  in_file.close();
}
			 

   
void add_step_info_to_pentamers_table(DNA_to_properties& onemap, std::string filename, std::string sequence,bool verbose, bool debug){
  std::ifstream cclis_inf(filename.c_str());
  if (cclis_inf){
    std::string line;
    std::vector<double_vector> matrix;
    matrix.clear();

    //add propel(omega) informaton
    while (getline(cclis_inf,line)){
      if (line.find("|D| AV")!=std::string::npos){
	//skip 7 lines
	for (int i=0;i<7;i++)
	  getline(cclis_inf,line);

	string_vector string_row;
	double_vector double_row;
	do{
	  getline(cclis_inf,line);
	  parse_string_to_list(line,"\t ",string_row);
	  if (string_row.size()==12){
	    //delete the first 4 elements of the vector
	    for (int j=0;j<4;j++)
	      string_row.erase(string_row.begin());
	    string_vector_to_double_vector(string_row,double_row);
	    matrix.push_back(double_row);
	  }
	  else if (string_row.size()>0){
	    std::cout << "Error: Cannot parse the line: "<<line<<std::endl;
	    return;
	  }
	}while (string_row.size()==8);
	break;
      }
    }
    if (debug){
      std::cout<<"Block |D|"<<std::endl;
      for (int i=0;i<matrix.size();i++){
	for (int j=0;j<matrix[i].size();j++)
	  std::cout << "\t"<<matrix[i][j];
	std::cout << std::endl;
      }
    }
    if (matrix.size()>0){
      add_propel_to_table(onemap,sequence,matrix,4,"propel",verbose,debug);

    matrix.clear();
    //add step information
    while (getline(cclis_inf,line)){
      if (line.find("|H| AV")!=std::string::npos){
	  //skip 7 lines
	  for (int i=0;i<7; i++)
	    getline(cclis_inf,line);

	  string_vector string_row;
	  double_vector double_row;
	  do{
	    getline(cclis_inf,line);
	    parse_string_to_list(line,"\t ",string_row);
	    if (string_row.size()==11){
	      //delete the first 4 elements of the vector
	      for (int j=0;j<4;j++)
		string_row.erase(string_row.begin());
	      string_vector_to_double_vector(string_row,double_row);
	      matrix.push_back(double_row);
	    }
	    else if (string_row.size()>0){
	      std::cout << "Cannot parse the line: "<<line<<std::endl;
	      exit(1);
	    }

	  }while (string_row.size()==7);
	  break;
      }
    }

    if (debug){
      std::cout<<"Block |H|:"<<std::endl;
      for (int i=0;i<matrix.size();i++){
	for (int j=0; j<matrix[i].size(); j++)
	  std::cout << "\t"<<matrix[i][j];
	std::cout <<std::endl;
      }
    }
    
    if (matrix.size()>0){
      add_one_step_info(onemap,sequence,matrix,1,"slide",verbose,debug);
      add_one_step_info(onemap,sequence,matrix,4,"roll",verbose,debug);
      add_one_step_info(onemap,sequence,matrix,5,"twist",verbose,debug);
    }    
  }
  else
    if (debug)
      std::cout<<"Cannot open the cclis file: " <<filename<<std::endl;

    cclis_inf.close();
  }
}



void add_one_step_info(DNA_to_properties& onemap,std::string sequence,std::vector<double_vector> &matrix,int object_index, \
		       std::string object_name, bool verbose, bool debug){
  if (matrix.size()!=sequence.size()-1){
    if (debug)
      std::cout << "Error: Size of step info matrix does not equal to (len(sequence)-1)"<<std::endl;
    return;
  }
  //discard the first 2 step and last 2
  // N - N -(step1) - N - (step2) - N     what is object_name1, what is object_name2
  std::string current_pentamer1,current_pentamer2;
  std::string object_name1 = object_name + "1";
  std::string object_name2 = object_name + "2";
  for (int i=2; i<matrix.size()-2; i++){
    current_pentamer1 = sequence.substr(i-1,5);
    current_pentamer2 = sequence.substr(i-2,5);
    
    if (found_str_in_map(current_pentamer1,onemap))
      onemap[current_pentamer1].push(matrix[i][object_index],object_name1);
    else
      onemap[opposite_strand(current_pentamer1)].push(matrix[i][object_index],object_name2);

    if (found_str_in_map(current_pentamer2,onemap))
      onemap[current_pentamer2].push(matrix[i][object_index],object_name2);
    else
      onemap[opposite_strand(current_pentamer2)].push(matrix[i][object_index],object_name1);
  }

}

void add_propel_to_table(DNA_to_properties& onemap,std::string sequence, std::vector<double_vector> &matrix,int object_index,\
			  std::string object_name,bool verbose,bool debug){
   if (matrix.size()!=sequence.size()){
     if (debug)
       std::cout << "Error: Size of Global Base-Base Parameters does not equal to len(sequence)"<<std::endl;
     return;
   }
   //discard the first 2 liens and last 2 
   std::string current_pentamer;
   for (int i=2; i<matrix.size()-2; i++){
     current_pentamer = sequence.substr(i-2,5);
     if (found_str_in_map(current_pentamer,onemap))
       onemap[current_pentamer].push(matrix[i][object_index],object_name);
     else
       onemap[opposite_strand(current_pentamer)].push(matrix[i][object_index],object_name);
   }
 }
  
void process_querytable_file(std::string querytable_filename,DNA_to_properties& onemap){
  std::ifstream qt_ifstream(querytable_filename.c_str());
  if (!qt_ifstream){     //check if the file exist
    std::cerr << "Cannot open the following file containing query table: " << querytable_filename <<"\n";
    exit(1);
  }
  std::string line,current_pentamer;
  string_vector sv;
  double_vector dv;
  while (getline(qt_ifstream,line)){
    if (line.size()>0){
      parse_string_to_list(line," ",sv);
      if (sv.size()==28){
	current_pentamer = sv[0];
	sv.erase(sv.begin());
	string_vector_to_double_vector(sv,dv);
	properties p = properties();
	p.load_data_from_vector(dv);
	onemap[current_pentamer] = p;	
      }
      else {
	std::cerr << "Cannot parse the following line:\n"<<line<<std::endl;
	exit(1);
      }      
    }
  }
}

void add_groove_width_to_inosine_table(std::string filename,DNA_to_properties &onemap,std::string sequence,int object_index, \
					 std::string object_name,bool verbose,bool debug,bool five){

  //construct the matrix from the file
  std::vector<double_vector> matrix;
  matrix.clear();
  std::ifstream in_file(filename.c_str());

  if (!in_file){
    if (debug)
      std::cout << "Cannot open the file: " << filename << std::endl;
    in_file.close();
    return;
  }

  std::string line;

  while (getline(in_file,line)){
    double_vector double_row;
    double_row.clear();
    string_vector string_row;
    parse_string_to_list(line,"\t ",string_row);		
    if ((string_row.size()>0) and (string_vector_to_double_vector(string_row,double_row)))
      matrix.push_back(double_row);
  }
  //if the file is empty, then return
  if (matrix.size()==0)
    return;

  //  Matrix definition
  //  0--ResId    1--Levels  3--Average minor/major groove width
  //  8--occurrences

  int no_of_positions = sequence.size();

  if (debug){
    for (int i=0;i<matrix.size();i++){
      for (int j=0;j<matrix[i].size();j++)
	std::cout <<"\t"<<matrix[i][j];
      std::cout<<"\n";
    }
    std::cout<<"No of positions: "<<no_of_positions<<std::endl;
  }

  std::vector<int> num_of_levels;     //num of levels for each ResId
  for (int i=0;i<=no_of_positions;i++)
    num_of_levels.push_back(0);
  
  //std::cout << "Vector size: "<<num_of_levels.size()<<std::endl;
  for (int i=0; i<matrix.size(); i++){
    //std::cout << double_to_int(matrix[i][0])  <<std::endl;
    num_of_levels[double_to_int(matrix[i][0])]++;
  }

 
  int start_pos = 4;  //discard the first 3 positions      //the ResID of the first sequence character is 1
  int end_pos=no_of_positions-4; // discard the last 3
  double current_value,prev_value,next_value,prevprev_value,nextnext_value;
  int current_id;
  std::string current_pentamer;
  for (int i=1;i<(matrix.size()-1);i++){
    if ((double_to_int(matrix[i][0])>=start_pos) and (double_to_int(matrix[i][0])<=end_pos) \
	and (double_to_int(matrix[i][1])==1) and (double_to_int(matrix[i][8])>100000)   //the first is level is 1
	and (double_to_int(matrix[i][0])==7)){   //CGCG NNNNN CGCG  , only utilize the information of central pentamer
            
      current_id = double_to_int(matrix[i][0]);
      current_pentamer = sequence.substr(current_id-1-2,5);  //1bp shift
      current_value = matrix[i][object_index];
      
      if ((double_to_int(matrix[i+1][0])==current_id) and (double_to_int(matrix[i+1][1])==2) \
	  and (double_to_int(matrix[i+1][8])>100000))
	next_value = matrix[i+1][object_index];
      else
	continue;
      
      if (five){  // Use 5 levels instead of 3 
	if ((i+2<matrix.size()) and (double_to_int(matrix[i+2][0])==current_id)	\
	    and (double_to_int(matrix[i+2][1])==3) and (double_to_int(matrix[i+2][8])>100000))
	  nextnext_value = matrix[i+2][object_index];
	else 
	  continue;
      }	    
      
      int tp = i;
      bool flag = false;
      for (int j=0;j<num_of_levels[current_id-1];j++){
	tp--;
	if ((double_to_int(matrix[tp][1])<=4) and (double_to_int(matrix[tp][8])>100000)){
	    flag = true;
	    prev_value = matrix[tp][object_index];
	    break;
	  }
      }
      if (!flag) continue;

      if (five){
	if ((tp-1>=0) and (double_to_int(matrix[tp-1][0])==current_id-1) and (double_to_int(matrix[tp-1][8])>100000)){
	  prevprev_value = matrix[tp-1][object_index];
	}
	else
	  continue;
      }

      double final_value;
      if (five)
	final_value = (prevprev_value + prev_value + current_value + next_value + nextnext_value)/5;
      else
	final_value = (prev_value + current_value + next_value)/3;

      //std::cout << "Final value: "<<final_value << std::endl;
      if (!found_str_in_map(current_pentamer,onemap)){
	properties p=properties();
	onemap[current_pentamer]=p;
      }
      
      onemap[current_pentamer].push(final_value,object_name);
      
    }
  }
  in_file.close();
}

