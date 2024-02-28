/* 
MIT License Copyright 2024 Francesco Savelli - see LICENSE file.

Generic functions used elsewhere in the codebase. 

*/

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "util.hpp"

using namespace std;

vector<string> readlines(const string& namefile){

  vector<string> result;
  string line;
  ifstream infile (namefile.c_str());

  while (getline(infile,line,'\n'))
    result.push_back(line);
  
  infile.close();
  return result;
}

vector<string> splitstring(const string& line, char sep){

  stringstream ss(line);
  string word;
  vector<string> words;
  while(std::getline(ss, word, sep))
     words.push_back(word);
  return words;
}

string itos(long long int number){
	stringstream ss;
	ss<<number;
	return ss.str();
}

string dtos(double number){
	stringstream ss;
	ss<<number;
	return ss.str();
}

/*
string_matrix_t& readstringparams(const string& namefile){

  string_matrix_t& result;

  ifstream infile (namefile.c_str());
  string line;
  while (getline(infile, line, '\n')){

    stringstream ss(line);
    string word;
    vector<string> words;
    while(getline(ss, word, " "))
       words.push_back(word);
    
    result.push_back(words);
  }
  
  infile.close();

  return result
}

void readlines(const string& namefile, vector<string>& result){

  //vector<string> lines;
  string line;
  result.clear();
  ifstream infile (namefile.c_str());

  while (getline(infile,line,'\n'))
    result.push_back(line);
  
  infile.close();
}


vector<string> splitstring(const string& input, const string& separator){

  stringstream ss(input);
  string word;
  vector<string> words;

  while(getline(ss, word, separator))
     words.push_back(segment);
  
  return words;
}



void splitline(const string& line, const string& sep,  vector<string>& words){

  words.clear();

  string::size_type current_pos=0;
  string::size_type next_pos=line.find(sep);
  while(next_pos!=string::npos){ 
    //npos is a special number returned when no sep is found in the string
    //string word=line.substr(0, next_pos);
    //string word=line.substr(current_pos, next_pos);
    //string::size_type begin=word.find_first_not_of(" ");
    string::size_type begin=line.find_first_not_of(" ", current_pos);
    string::size_type end=line.find_last_not_of(" ", next_pos-1);
    if ((begin==string::npos)&&(end==string::npos)&&((end-begin+1)<=0))
      cerr<<"Error in reading line: "<<line<<endl;
    else
      words.push_back(string(line, begin, end-begin+1));
    current_pos=next_pos+1;
    next_pos=line.find(sep, current_pos);
  }

  string::size_type begin=line.find_first_not_of(" ", current_pos);
  string::size_type end=line.find_last_not_of(" ");
  if ((begin==string::npos)&&(end==string::npos)&&((end-begin+1)<=0))
    cerr<<"Error in reading line: "<<line<<endl;
  else
    words.push_back(string(line, begin, end-begin+1));
}

void readstringmatrix(const string& namefile, const string& sep, 
                      string_matrix_t& result){

  result.clear();
  vector<string> lines;
  readlines(namefile, lines);
  vector<string> words;
  for (unsigned int line_ind=0; line_ind<lines.size(); line_ind++){
    splitline(lines[line_ind], sep, words);
    result.push_back(words); 
    // push_back makes a deep copy, so do not worry about words being rewritten 
  }
}
*/


