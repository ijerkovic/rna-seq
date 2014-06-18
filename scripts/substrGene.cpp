#include <vector>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;

string whole_genome;
/*
 * Extracts the whole genome from genome file and stores it into a string
 */
void extractGenomeFromFile(string genomeFile) {

  whole_genome = "";

  FILE *inGenome = fopen(genomeFile.c_str(), "r");
  char buff[256];
  while (fgets(buff, 255, inGenome)) {
    // ignore lines that start with genome desc, they start with '>'
    if (buff[0] != '>') {
      string tmp = buff;
      tmp.resize(tmp.size() - 1);  // remove endl
      whole_genome += tmp;
    }
  }
  // cout << "WHOLE " << whole_genome.size() << endl;
}

map <string, int> STATS;

void getDonorAcceptorStats(string juncsFilename) {
  ifstream juncsFile(juncsFilename.c_str());
  string line;

  while(getline(juncsFile, line)) {
    int pos1, pos2;
    for (int i = 0; i < line.size(); i++) {
      if (line[i] == ':') pos1 = i;
      if (line[i] == '-') pos2 = i;
    }
    pos1++;
    string leftnum = line.substr(pos1, pos2 - pos1);
    string rightnum = line.substr(pos2 + 1);
  
    int start = atoi(leftnum.c_str());
    int end = atoi(rightnum.c_str()) - 2;

    string intron = whole_genome.substr( start, end - start + 1 );
    
    //cout << intron << endl;
    
    string acc_don = intron.substr(0, 2) + "-" + intron.substr(intron.size() - 2);
    std::transform(acc_don.begin(), acc_don.end(), acc_don.begin(), ::toupper);
    //cout << acc_don << endl;
    STATS[acc_don]++;
  }

  typedef std::map<string, int>::iterator it_type;
 
  vector < pair<int, string> > acc_don_stats;
  int all = 0;
  for(it_type iterator = STATS.begin(); iterator != STATS.end(); iterator++) {
    //cout << iterator->first << " ";
    //cout << iterator->second << endl;
    all += iterator->second;
    acc_don_stats.push_back(make_pair(iterator->second, iterator->first));
  }

  sort(acc_don_stats.rbegin(), acc_don_stats.rend());

  cout << "<table>" << endl;
  for (int i = 0; i < acc_don_stats.size(); i++) {
    cout << "<tr>" << endl;
    cout << "<td>" << acc_don_stats[i].second << "</td><td>" << acc_don_stats[i].first << "</td><td>" << (double)acc_don_stats[i].first / all * 100. << "%</td>" << endl;
    cout << "</tr>" << endl;
  }
  cout << "</table>" << endl;

}

int main(int argc, char *argv[]) {

  if (argc < 4) {
    cout << "usage: substrGene <gene.fa> <start> <end> [--junc-stats junc-reads-file]" << endl;
    return -1;
  }
  
  extractGenomeFromFile(argv[1]);
  
  string juncsFile = "";

  for (int i = 0; i < argc; i++) { 
    if (strcmp(argv[i], "--junc-stats") == 0) {
      juncsFile = argv[i + 1];
    }
  }

  if (juncsFile != "") {
    getDonorAcceptorStats(juncsFile);
    return 0;  
  }

  int start = atoi(argv[2]) - 1;
  int end = atoi(argv[3]) - 1;

  cout << whole_genome.substr(start, end - start + 1) << endl;

  return 0; 
}
