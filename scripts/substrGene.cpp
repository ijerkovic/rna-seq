#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
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

int main(int argc, char *argv[]) {

  if (argc != 4) {
    cout << "usage: substrGene <gene.fa> <start> <end>" << endl;
  }

  extractGenomeFromFile(argv[1]);

  int start = atoi(argv[2]) - 1;
  int end = atoi(argv[3]) - 1;

  cout << whole_genome.substr(start, end - start + 1) << endl;

  return 0; 
}
