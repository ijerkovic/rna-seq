  // Copyright 2014 Igor Jerkovic

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include <string>
using namespace std;

#define MAX(A, B) (((A) > (B)) ? (A) : (B))
#define MIN(A, B) (((A) < (B)) ? (A) : (B))

#define GET_DATA_FOR_HISTOGRAM 1

  // set to true if there is no need to do the build-index phase again
  int skip_build_index = 0;

  int read_length = 50;

  int coverage_limit = 20;

  struct Read {
    string id;
    string value;
    int map_position;
    Read(string id_, string value_, int map_position_) {
      id = id_; value = value_; map_position = map_position_;  
    }

    bool friend operator < (const Read &A, const Read &B) {
      return A.map_position < B.map_position;  
    }
  };

#define MAX_EXON_LEN 8192

  struct Exon {
    int start;
    int left;
    int right;
    int coverage[MAX_EXON_LEN];
    int c;

    Exon() {}
    Exon(int left_, int right_) {
      start = left_;
      left = left_;
      right = right_;
      c = 0;
      for (int i = 0; i < MAX_EXON_LEN; i++)
        coverage[i] = 0;
    }
  };

  vector <Exon> exon_regions;
  vector <Exon> possible_exon_regions;

  void addCoverage(Exon &e, int start, int end) {
    e.c++;
    if (end >= MAX_EXON_LEN) cout << "YES";
    for (int i = start; i <= MIN(end, MAX_EXON_LEN - 1); i++) {
      e.coverage[i]++;
    }
  }

  double calcCoverage(Exon &e) {
    double cov = 0;
    for (int i = 0; i < MAX_EXON_LEN; i++) {
      cov += e.coverage[i];
    }
    return cov / (double)(e.right - e.left);
  }

  int newE = 0, addedE = 0;

  void insertIntoExonRegions(int left) {
    int found = 0;
    for (int i = 0; i < exon_regions.size(); i++) {
      if (left >= exon_regions[i].left && left <= exon_regions[i].right) {
        exon_regions[i].right = MAX(left + read_length, exon_regions[i].right);
        found = 1;
        addedE++;
        addCoverage(exon_regions[i],
                    left - exon_regions[i].start,
                    left + read_length - exon_regions[i].start);
      }
    }
    if (!found) {
      newE++;
      exon_regions.push_back(Exon(left, left + read_length));
      addCoverage(exon_regions[exon_regions.size() - 1], 0, read_length);
    }
  }


  vector <Read> matched_reads;
  vector <Read> unmatched_reads;

  /**
   * builds index that is later on used for genome alignment
   * genomeFileName - filename of genome for which index is being built
   * indexPrefix - some aligners need a prefix name for additional files
   * aligner - for now bwa, bowtie, bowtie2 and lisa are supported
   */
  void buildIndexForGenomeAlignment(string genomeFileName, string indexPrefix, string aligner) {
    if (aligner == "bwa") {
      string build_index = "bwa index ";
      build_index += genomeFileName;
      system(build_index.c_str());
    } else if (aligner == "bowtie" || aligner == "bowtie2") {
      string build_index = aligner;
      build_index += "-build ";
      build_index += genomeFileName;
      build_index += " ";
      build_index += indexPrefix;
      cout << build_index << endl;
      system(build_index.c_str());
    }
  }

  /**
   * aligns read against the prebuilt index
   * indexPrefix - prefix name used by some aligners
   * readName - filename of the read that is going to be aligned
   * aligner - bwa, bowtie, bowtie2 or lisa
   * resultFile - filename for the result alignment
   */
  void alignReadAgainstGenome(string indexPrefix, string readName, string aligner, string resultFile) {
    if (aligner == "bwa") {
      string mapper = "bwa mem ";
      mapper += indexPrefix;
      mapper += " ";
      mapper += readName;
      mapper += " > ";
      mapper += resultFile;
      system(mapper.c_str());
    } else if (aligner == "bowtie" || aligner == "bowtie2") {
      string mapper = aligner;
      mapper += " ";
      mapper += indexPrefix;
      mapper += " ";
      mapper += readName;
      mapper += " > ";
      mapper += resultFile;
      cout << mapper << endl;
      system(mapper.c_str());  
    }
  }

  /*
   * checks input arguments passed from main
   */
  int checkArguments(int argc, char *argv[]) {

    if (argc < 5) {
      cout << "Usage: rna-seq <aligner> <genome.fa> <index_dir> <reads.fq> <result_file> <optional-gtf-file> --cov-limit [coverage-limit-number] [--skip-building-index]" << endl;
      return -1;
    }

    for(int i = 0; i < argc; i++) {
      if(strcmp(argv[i], "--skip-building-index") == 0) {
        skip_build_index = 1;
      }
    }

    string alignerName = argv[1];
    if((alignerName != "bwa") && (alignerName != "bowtie") && (alignerName != "bowtie2") && 
       (alignerName != "lisa")) {
      cout << "supported aligners are 'bwa', 'bowtie', 'bowtie2' and 'lisa'" << endl;
      return -1;  
    }

  }

  /*
   * Gets read alignement data that is used in other phases of alignment
   * samFilename - filename in sam format that contains read alignment data
   */
  void extractReadAlignmentDataFromSamFile(string samFilename) {
    ifstream samFile(samFilename.c_str());
    string line;
    
    int num_reads = 0;
    int num_unmapped = 0, num_mapped = 0;

    while(getline(samFile, line)) {
      if(line[0] == '@') {
        // skip header line from sam file
        continue;  
      }
      num_reads++;

      istringstream iss(line);
      string col;
      vector <string> cols;
      while(iss >> col) {
        cols.push_back(col);
      }

      if(cols[3] == "0") {
        num_unmapped++;
        unmatched_reads.push_back(Read(cols[0], cols[9], 0));
      } else {
        num_mapped++;
        matched_reads.push_back(Read(cols[0], cols[9], atoi(cols[3].c_str())));
      }
    }

    cout << "number of reads: " << num_reads << endl;
    cout << "number of mapped reads: " << num_mapped << " (" << (double)num_mapped/num_reads*100 << ")" << endl;
    cout << "number of unmapped reads: " << num_unmapped << " (" << (double)num_unmapped/num_reads*100 << ")" << endl;
  }

#define HIST_SIZE 60000000
  int hist[HIST_SIZE] = {0};

  void prepareDataForHistogramDrawing() {
    int in = 0, out = 0;
    for (int i = 0; i < matched_reads.size(); i++) {
      if(matched_reads[i].map_position >= 3000000 && matched_reads[i].map_position <= 3500000) in++;
      else out++;
      for (int j = matched_reads[i].map_position; j < matched_reads[i].map_position + read_length; j++) {
        hist[j]++;  
      }  
    }

    FILE *hist_f = fopen("hist.in", "w");
    cout << "JURO " << matched_reads[0].map_position << endl;
    for (int i = matched_reads[0].map_position; i < HIST_SIZE; i++)
      if (i >= 3000000 && i <= 3500000 ) {
        fprintf(hist_f, "%d\n", hist[i]);
      }
    cout << in << " vs " << out << endl;

    cout << matched_reads[0].map_position << " " << matched_reads[matched_reads.size() - 1].map_position << endl;
  }

  void calculatePossibleExonRegions() {
    
    int total = matched_reads.size() + unmatched_reads.size();
  
  sort(matched_reads.begin(), matched_reads.end());

  if (GET_DATA_FOR_HISTOGRAM) {
    prepareDataForHistogramDrawing();  
  }

  for (int i = 0; i < matched_reads.size(); i++) {
    insertIntoExonRegions(matched_reads[i].map_position);
  }

  int all = 0;
  for (int i = 0; i < exon_regions.size(); i++) {
    if (calcCoverage(exon_regions[i]) > coverage_limit) {
      all += exon_regions[i].c;
      possible_exon_regions.push_back(exon_regions[i]);
    }
  }

  for (int i = 0; i < matched_reads.size(); i++) {
    int l = matched_reads[i].map_position;
    int r = l + read_length - 1;

    int found = 0;
    // this should be done in binary search
    for (int j = 0; j < possible_exon_regions.size(); j++) {
      if (l >= possible_exon_regions[ j ].left && r <=
          possible_exon_regions[ j ].right) {
        found = 1;
        break;
      }
    }

    if (!found) {
      unmatched_reads.push_back(matched_reads[i]);
    }
  }

  cout << "AFTER INITIAL EXON COVERAGE SEARCH" << endl;
  cout << "MAPPED ON EXONS: " << total - unmatched_reads.size() << " (";
  cout <<  (total - unmatched_reads.size()) / (double)total * 100 << ")" << endl;
  cout << "UNMAPPED READS: " << unmatched_reads.size() << " (";
  cout <<  unmatched_reads.size() / (double)total * 100 << ")" << endl;
}

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
  cout << "WHOLE " << whole_genome.size() << endl;
 }

vector <int> spliceSites;
string exon_genome;
/*
 * Creates one big string that contains bases from original genome that are calculated
 * to be exon regions
 */
void connectExonRegions() {
  exon_genome = "";
  spliceSites.clear();

  cout << "w2 " << whole_genome.size() << endl;
  
  for (int i = 0; i < possible_exon_regions.size(); i++) {
    cout << possible_exon_regions[i].left << " " << possible_exon_regions[i].right << ", LENGTH= " << possible_exon_regions[i].right - possible_exon_regions[i].left << endl;
    /*cout << 
      whole_genome.substr(
        possible_exon_regions[i].left,
        possible_exon_regions[i].right - possible_exon_regions[i].left + 1) << endl;
        cout << endl;*/

    exon_genome +=
      whole_genome.substr(
        possible_exon_regions[i].left,
        possible_exon_regions[i].right - possible_exon_regions[i].left + 1);

    spliceSites.push_back(exon_genome.size());
  }

  if (spliceSites.size()) spliceSites.pop_back();
 
  cout << "original genome length: " << whole_genome.length() << endl;
  cout << "exon genome length: " << exon_genome.length() << " (";
  cout << exon_genome.length() / (double)whole_genome.length() * 100 << "% of original size)" << endl;
}

//#define INDEX_SIZE (1LL<<32)
//vector <int> exon_index[INDEX_SIZE];

map <string, vector<int> > exon_index;

#define SLIDING_WINDOW_SIZE 13

/*
 * calculates hash based on input string
 * characters of key can be 'A', 'C', 'T', 'G' or 'N'
 * when 'A' 0 is assigned, when 'N' assuming it is 'A' since we don't know the exact value
 * other chars get 1, 2 and 3 assigned respectivly
 */
unsigned int hash(string key) {

  assert(key.size() == SLIDING_WINDOW_SIZE);
  
  unsigned int h = 0;
  unsigned int pot = 1;
  for (int i = 0; i < key.size(); i++) {
    char c = toupper(key[i]);
    if (c == 'C') {
      h += 1 * pot;
    } else if (c == 'T') {
      h += 2 * pot;  
    } else if (c == 'G') {
      h += 3 * pot;  
    } // if A or N: h += 0
    pot *= 4;
  }
  
  return h;
}

/*
 * builds an index for the whole exon genome
 * exon genome is considerably smaller than original genome since exons are much smaller
 * than introns and introns are separated out
 */
void createIndex() {
  
  for (int i = 0; i < exon_genome.size() - SLIDING_WINDOW_SIZE + 1; i++) {
    string window = exon_genome.substr(i, SLIDING_WINDOW_SIZE);
    exon_index[window].push_back(i);
    //exon_index[ hash(window) ].push_back(i);
  }

}

void printIndexValue(string hashValue) {
  cout << "(";
  vector <int> v = exon_index[hashValue];
  for (int i = 0; i < v.size(); i++ ) {
    cout << v[i] << "," ; 
  }
  cout << ")";
}

/*
 * Matches the rest of reads against the exon index
 */
void unmatchedReadsAgainstIndex() {
  
  for (int i = 0; i < unmatched_reads.size(); i++) {
    
    string read = unmatched_reads[i].value;
    cout << read << endl;
    for (int j = 0; j < read.size() - SLIDING_WINDOW_SIZE + 1; j++) {
      //printIndexValue(hash(read.substr(j, SLIDING_WINDOW_SIZE)));
      printIndexValue(read.substr(j, SLIDING_WINDOW_SIZE));
    }
    cout << endl;
    
  }

}

/*
 * Main method of rna-seq
 */
int main(int argc, char *argv[]) {

  if (checkArguments(argc, argv) == -1) return -1;

  string alignerName = argv[1];
  string refGenomeName = argv[2];
  string indexPrefix = argv[3];
  string readInputName = argv[4];
  string resultFile = argv[5];

  cout << skip_build_index << endl;

  if (!skip_build_index) {
    buildIndexForGenomeAlignment(refGenomeName, indexPrefix, alignerName);
  }
  
  alignReadAgainstGenome(indexPrefix, readInputName, alignerName, resultFile);

  extractReadAlignmentDataFromSamFile(resultFile);

  calculatePossibleExonRegions();

  extractGenomeFromFile(refGenomeName);

  connectExonRegions();

  createIndex();

  unmatchedReadsAgainstIndex();

  FILE *all_reads_file = fopen(readInputName.c_str(), "r");
  FILE *unmatched_reads_file = fopen("unmatched.fastq", "w");
  
#define MAX_READ_LENGTH 10000
  
  char buff1[MAX_READ_LENGTH];
  char buff2[MAX_READ_LENGTH];
  char buff3[MAX_READ_LENGTH];
  char buff4[MAX_READ_LENGTH];
  char buff5[MAX_READ_LENGTH];

  while (fgets(buff1, MAX_READ_LENGTH, all_reads_file)) {
    fgets(buff2, MAX_READ_LENGTH, all_reads_file);

    string read_id = buff1;
    sscanf(buff1, "%s", buff5);
    read_id = buff5;
    read_id = read_id.substr(1);

    string read_val = buff2;
    read_val.resize(read_val.size() - 1);

    fgets(buff3, MAX_READ_LENGTH, all_reads_file);
    fgets(buff4, MAX_READ_LENGTH, all_reads_file);
    for (int i = 0; i < unmatched_reads.size(); i++) {
      //cout << unmatched_reads[0].id << " " << read_id << endl;
      if (unmatched_reads[i].id == read_id) {
        fputs(buff1, unmatched_reads_file);
        fputs(buff2, unmatched_reads_file);
        fputs(buff3, unmatched_reads_file);
        fputs(buff4, unmatched_reads_file);
      break;
      }
    }
  }
  /*for (int i = 0; i < unmatched_reads.size(); i++) {
    fprintf(unmatched_reads_file, "%s HWI-BRUNOP16X_0001:3:1:%d:%d#0/1\n", unmatched_reads[i].id.c_str(), rand(), rand());
    fprintf(unmatched_reads_file, "%s\n", unmatched_reads[i].value.c_str());
    string tmp = string(unmatched_reads[i].value.size(), '#');
    fprintf(unmatched_reads_file, "+\n%s\n", tmp.c_str());
  }*/

  FILE *outGenome = fopen("tmp.fa", "w");
  fprintf(outGenome, ">chr19\n");
  for (int i = 0; i < exon_genome.size(); i+=80) {
    fprintf(outGenome,
            "%s\n",
            exon_genome.substr(i, MIN(80, exon_genome.size() - i)).c_str());
  }

  buildIndexForGenomeAlignment("tmp.fa", "tmp", alignerName);
  alignReadAgainstGenome("tmp", "unmatched.fastq", alignerName, "tmp.sam");

  return 0;  
}
