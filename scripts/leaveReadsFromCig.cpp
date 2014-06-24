#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
using namespace std;

map <string, int> MAPA;

int main(int argc, char *argv[]) {
  
  if (argc != 3) {
    cout << "leaveReadsFromCig <file.cig> <reads.fq>" << endl;
    return -1;  
  }

  ifstream cig_in(argv[1]);
  string line;

  while (getline(cig_in, line)) {
    
    istringstream iss(line);
    string col;

    vector <string> cols;
    while (iss >> col) {
      cols.push_back(col);  
    }

    MAPA[ "@" + cols[0] ] = 1;

  }

  ifstream reads_in(argv[2]);
  string line1, line2, line3, line4;

  while (getline(reads_in, line1)) {
    
    getline(reads_in, line2);
    getline(reads_in, line3);
    getline(reads_in, line4);
  
    if (MAPA[line1]) {
      cout << line1 << endl;
      cout << line2 << endl;  
      cout << line3 << endl;  
      cout << line4 << endl;  
    }

  }

  return 0;
}
