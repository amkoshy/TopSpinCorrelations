#include "ttmd_fileReader.h"


std::vector<std::string> split(const std::string& s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

// bool BothAreSpaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }
