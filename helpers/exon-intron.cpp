#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
using namespace std;

int main(int argc, char *argv[] ) {

  if(argc != 2) {
    cout << "Usage: ./exon-intron genome-file.fa" << endl;
    return -1;
  }

  FILE *in = fopen(argv[1], "r");

  int wholeGenomeLength = 0;
  char line[256];

  while(fgets(line, 255, in)) {
    // ignore lines that start with genome desc, they start with '>'
    if(line[0] != '>') {
      wholeGenomeLength += strlen(line) - 1;
    }
  }

  srand(time(NULL));

  FILE *exonRegions = fopen(strcat(argv[1], "-exon-regions"), "w+");

  for(int i = rand() % 1000; i < wholeGenomeLength; ) {
    int left = i;
    int right = left + rand() % 1000 + 100;
    if (right >= wholeGenomeLength) {
      break;  
    }
    fprintf(exonRegions, "%d %d\n", left, right);
    i += (right - left) * (rand() % 900 + 100);
  }

  return 0;
}
