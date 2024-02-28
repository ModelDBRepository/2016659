/* 
MIT License Copyright 2024 Francesco Savelli - see LICENSE file.

Generic functions used elsewhere in the codebase. 

*/

#ifndef UTIL_HPP
#define UTIL_HPP

#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

/*
typedef vector<vector<string> > string_matrix_t;

string_matrix_t& readstringparams(const string& namefile);

void readlines(const string& namefile, vector<string>& result);

void splitline(const string& line, const string& sep,  vector<string>& words);

void readstringmatrix(const string& namefile, const string& sep, 
                      string_matrix_t& result);
*/

// IO utilities
vector<string> readlines(const string& namefile);

vector<string> splitstring(const string& line, char sep);

// Conversions
string itos(long long int number);

string dtos(double number);

#endif

