#include <vector>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
using namespace std;

vector <string> readIds;

map <string, pair <int, int > > CIG;

int is_paired_end = 0;

void compareAlns(string cigFile, string bamFile) {

  ifstream cig_in(cigFile.c_str());
  ifstream bam_in(bamFile.c_str());
  string line;
  int ll = 0;

  while (std::getline(cig_in, line)) {

    istringstream iss(line);
    string col;
    vector <string> cols;
    while(iss >> col) {
      cols.push_back(col);
    }

    int first, last;
    if (cols[4][cols[4].size() - 1] != ',') {
      string aln1 = cols[4];
      for (int i = 0; i < aln1.size(); i++) {
        if (aln1[i] == '-') {
          first = atoi( (aln1.substr(0, i)).c_str() );
          last = atoi( (aln1.substr(i + 1)).c_str() );
          break;
        }  
      }  
    } else {
      string aln1 = cols[4];
      for (int i = 0; i < aln1.size(); i++) {
        if (aln1[i] == '-') {
          first = atoi( (aln1.substr(0, i)).c_str() );
          break;
        }  
      }
      int k = 4;
      while (cols[k][cols[k].size() - 1] == ',') {
        k++;  
      }
      aln1 = cols[k];
      for (int i = 0; i < aln1.size(); i++) {
        if (aln1[i] == '-') {
          last = atoi( (aln1.substr(i + 1)).c_str() );
          break;
        }  
      }
    }

    CIG[ cols[0] ] = make_pair(first, last);
  
    //cout << cols[0] << " " << first << " " << last << endl;
  
    ll++;
  }

  int hit = 0, reported_reads = 0;
  while (std::getline(bam_in, line)) {

    istringstream iss(line);
    string col;
    vector <string> cols;
    while(iss >> col) {
      cols.push_back(col);
    }

    reported_reads++;

    pair <int, int> par = CIG[cols[3]];
    if (par.first - 1 == atoi(cols[1].c_str()) || par.second == atoi(cols[2].c_str())) {
      hit++;  
    } else {
      //cout << par.first << " " << atoi(cols[1].c_str()) << " " << par.second << " " << atoi(cols[2].c_str()) << " " << cols[3] << endl;  
    } 
  
  }
  
  if (!is_paired_end) {
    ll /= 2;  
  }

  cout << "OUT OF " << reported_reads << " ALL READS " << hit << " ALIGNED CORRECTLY ";
  cout << (double)hit/reported_reads * 100 << "%" << endl;

  cout << "OUT OF " << ll << " ALL READS IN CIG FILE " << hit << " ALIGNED CORRECTLY ";
  cout << (double)hit/ll * 100 << "%" << endl;
}

int main(int argc, char *argv[]) {

  if (argc < 3) {
    cout << "usage: bamVsBeers <beers.cig> <aln.bam>" << endl;
    cout << "options:" << endl;
    cout << "--paired-end" << endl;
    return -1;
  }

  if (argc >= 4 && strcmp(argv[3], "--paired-end") == 0) {
    is_paired_end = 1;
  }

  compareAlns(argv[1], argv[2]);

  return 0;
}
