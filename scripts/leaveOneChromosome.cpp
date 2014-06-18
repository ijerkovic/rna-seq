#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

vector <string> readIds;

void leaveOnlyOneChrFromFile(string inputFile, string chr) {

  ifstream input_file(inputFile.c_str());
  string line;

  while (std::getline(input_file, line)) {

    istringstream iss(line);
    string col;
    vector <string> cols;
    while(iss >> col) {
      cols.push_back(col);
    }

    if (cols[1] == chr) {
      readIds.push_back(cols[0]);
      cout << line << endl;  
    }

  }

}

void leaveMatchedReadIds(string readFile) {
    
  ifstream input_file(readFile.c_str());
  string line;

  int num = 0;

  FILE *tmp1 = fopen("sim_1.tmp", "w");
  FILE *tmp2 = fopen("sim_2.tmp", "w");

  while (std::getline(input_file, line)) {

    if (num >= readIds.size()) break;

    string line2;
    std::getline(input_file, line2);

    if (line.substr(1, readIds[num].size()) == readIds[num]) {
      //cout << line << endl;
      //cout << line2 << endl;
      num++;
      if (line[line.size() - 1] == 'a') {
        fprintf(tmp1, "%s\n", line.c_str());
        fprintf(tmp1, "%s\n", line2.c_str());
      } else {
        fprintf(tmp2, "%s\n", line.c_str());
        fprintf(tmp2, "%s\n", line2.c_str());
      }
    }

  }

  fclose(tmp1);
  fclose(tmp2);

  // create fastq files

  string sys1 = "perl ~/rna-seq/scripts/fasta_to_fastq.pl sim_1.tmp > ";
  sys1 += readFile + "_1.fq";
  system(sys1.c_str());

  string sys2 = "perl ~/rna-seq/scripts/fasta_to_fastq.pl sim_1.tmp > ";
  sys2 += readFile + "_2.fq";
  system(sys2.c_str());

  system("rm sim_1.tmp sim_2.tmp");
}

int main(int argc, char *argv[]) {

  if (argc != 4) {
    cout << "usage: leaveOneChr <reads.cig> <chr-name eg. chr19> <reads.fa>" << endl;
    return -1;
  }

  leaveOnlyOneChrFromFile(argv[1], argv[2]);
  leaveMatchedReadIds(argv[3]);

  return 0;
}
