#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <vector>
using namespace std;

int main(int argc, char *argv[] ) {

  if(argc != 4) {
    cout << "Usage: ./get_exon_reads genome-file.fa exon-regions-file read_length" << endl;
    return -1;
  }

  FILE *genomeFile = fopen(argv[1], "r");
  FILE *exonRegionsFile = fopen(argv[2], "r");

  int read_length = atoi(argv[3]);

  if (read_length > 1000) {
    read_length = 1000;  
  }

  if (read_length < 20) {
    read_length = 20;  
  }

  vector < pair < int, int > > exon_regions;

  int left, right;
  while(fscanf(exonRegionsFile, "%d %d", &left, &right) != EOF) {
    cout << left << " " << right << endl;
    exon_regions.push_back(make_pair(left, right));
  }

  int wholeGenomeLength = 0;
  char line[256];

  FILE *reads = fopen(strcat(argv[1], "-reads.fq"), "w+");

  string allGenome = "";
  while(fgets(line, 255, genomeFile)) {
    // ignore lines that are genome desc, they start with '>'
    if(line[0] != '>') {
      string tmp = line;
      tmp.resize(tmp.size() - 1);
      allGenome += tmp;
    }
  }

  string fq_4row = string(read_length, 'E');

  int num_reads = 0;
  string tmp = "";
  for(int i = 0; i < exon_regions.size(); i++) {
    for(int j = exon_regions[ i ].first; j <= exon_regions[ i ].second; j++) {
      if( j % 10000 == 0) cout << j << endl;
      if(tmp.size() < read_length) {
        tmp += allGenome[ j ];
      }
      if(tmp.size() == read_length) {
        fprintf(reads, "@test_mRNA_%d_%d\n", i, j);
        fprintf(reads, "%s\n", tmp.c_str());
        fprintf(reads, "+\n%s\n", fq_4row.c_str());
        tmp = tmp.substr(1);
        num_reads++;
      }
    }  
  }

  cout << "NUMBER OF READS GENERATED: " << num_reads << endl;

  return 0;
}
