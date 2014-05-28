// Copyright 2014 Igor Jerkovic

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

struct Exon {
  int start;
  int left;
  int right;
  int coverage[1024];
  int c;

  Exon() {}
  Exon(int left_, int right_) {
    start = left_;
    left = left_;
    right = right_;
    c = 0;
    for (int i = 0; i < 1024; i++)
      coverage[i] = 0;
  }
};

vector <Exon> exon_regions;
vector < Exon > possible_exon_regions;

void addCoverage(Exon &e, int start, int end) {
  e.c++;
  for (int i = start; i <= MIN(end, 1023); i++) {
    e.coverage[i]++;
  }
}

double calcCoverage(Exon &e) {
  double cov = 0;
  for (int i = 0; i < 1024; i++) {
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
vector <string> unmatched_reads_refs;

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
      unmatched_reads_refs.push_back(cols[0]);
    } else {
      num_mapped++;
      matched_reads.push_back(Read(cols[0], cols[9], atoi(cols[3].c_str())));
    }
  }

  cout << "number of reads: " << num_reads << endl;
  cout << "number of mapped reads: " << num_mapped << " (" << (double)num_mapped/num_reads*100 << ")" << endl;
  cout << "number of unmapped reads: " << num_unmapped << " (" << (double)num_unmapped/num_reads*100 << ")" << endl;
}

void calculatePossibleExonRegions() {
  
  sort(matched_reads.begin(), matched_reads.end());

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
      unmatched_reads_refs.push_back(matched_reads[i].id);
    }
  }

  cout << unmatched_reads_refs.size() << endl;
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

  return 0;  
}
