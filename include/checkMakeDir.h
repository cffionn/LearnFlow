#ifndef CHECKMAKEDIR_H
#define CHECKMAKEDIR_H

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

struct stat st;

bool checkDir(const std::string inPath)
{
  if(stat(inPath.c_str(), &st) != 0) return false;

  if(st.st_mode & S_IFDIR) return true;
  else return false;
}

bool checkFile(const std::string inFile)
{
  if(stat(inFile.c_str(), &st) != 0) return false;

  if(st.st_mode & S_IFREG) return true;
  else return false;
}

bool checkFileExt(const std::string inFile, const std::string ext)
{
  if(!checkFile(inFile)) return false;
  
  if(ext.size() == 0){
    std::cout << "Given extension \'" << ext << "\' is invalid. return false" << std::endl;
    return false;
  }
  else if(inFile.size() < ext.size()){
    std::cout << "Given file \'" << inFile << "\' does not have given extension \'" << ext << "\'. return false" << std::endl;
    return false;
  }
  else if(inFile.substr(inFile.size() - ext.size(), ext.size()).find(ext.c_str()) == std::string::npos){
    std::cout << "Given file \'" << inFile << "\' does not have given extension \'" << ext << "\'. return false" << std::endl;
    return false;
  }

  return true;
}


bool checkMakeDir(const std::string inPath)
{
  bool dirIsMade = false;
  
  if(checkFile(inPath)){
    std::cout << "Path \'" << inPath << "\' already exists as file. Return false." << std::endl;
  }
  else if(checkDir(inPath)) dirIsMade = true;
  else if(!checkDir(inPath)){
    mkdir(inPath.c_str(), 0700);
    dirIsMade = true;
  }

  return dirIsMade;
}

void invalidFileMessage(const std::string inFile)
{
  std::cout << "Input \'" << inFile << "\' is not a valid file. return" << std::endl;
  return;
}

#endif
