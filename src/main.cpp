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
#include <ctime>
using namespace std;

#define MAX(A, B) (((A) > (B)) ? (A) : (B))
#define MIN(A, B) (((A) < (B)) ? (A) : (B))

int GET_DATA_FOR_HISTOGRAM = 0;

// set to true if there is no need to do the build-index phase again
int skip_build_index = 0;

int read_length = 100;

int coverage_limit = 1;

int hasResultCigFile = 0;

string result_cig_file;

string whole_genome;

int SLIDING_WINDOW_SIZE = 16;

int num_threads = 1;
  
int join_threshold = 100;
int extend_threshold = 50;

map <string, int> SEQ_MAP_SOL, SEQ_MAP_RES;

FILE *final_result;

struct Read {
  string id;
  string value;
  int map_position;
  int length;
  string chrName;
  Read(string id_, string value_, int map_position_) {
    id = id_; value = value_; map_position = map_position_, length = value_.size();
  }

  bool friend operator < (const Read &A, const Read &B) {
    return A.map_position < B.map_position;  
  }
};

#define MAX_EXON_LEN 8192*2

struct Exon {
  int start;
  int left;
  int right;
  int coverage[MAX_EXON_LEN];
  int c;
  string chrName;

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

void printExon(Exon &e) {
  cout << whole_genome.substr(e.left - 1, e.right - e.left - 1) << endl;  
}

int covm;

double calcCoverage(Exon &e) {
  double cov = 0;
  for (int i = 0; i < MAX_EXON_LEN; i++) {
    if (e.coverage[i] > covm) covm = e.coverage[i];
    cov += e.coverage[i];
  }
  return cov / (double)(e.right - e.left);
}

int newE = 0, addedE = 0;

void insertIntoExonRegions(Read &r) {
  int left = r.map_position;
  int found = 0;
  for (int i = 0; i < exon_regions.size(); i++) {
    if (left >= exon_regions[i].left && left <= exon_regions[i].right) {
      exon_regions[i].right = MAX(left + r.length, exon_regions[i].right);
      found = 1;
      addedE++;
      addCoverage(exon_regions[i],
          left - exon_regions[i].start,
          left + r.length - exon_regions[i].start);
    }
  }
  if (!found) {
    newE++;
    exon_regions.push_back(Exon(left, left + r.length));
    exon_regions[exon_regions.size() - 1].chrName = r.chrName;
    addCoverage(exon_regions[exon_regions.size() - 1], 0, r.length);
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
  } else if (aligner == "lisa") {
    string lisa_build_index = "client index ";
    lisa_build_index += genomeFileName;
    lisa_build_index += " ";
    lisa_build_index += indexPrefix;
    lisa_build_index += "\n";
    system(lisa_build_index.c_str());
  }
}

map <string, int> duplicate_reads;

void removeDuplicates(string readName, string newReadsFile) {
  
  ifstream readsFile(readName.c_str());
  string line1, line2, line3, line4;
  int n = 0, d = 0;

  FILE *withoutDupReads = fopen(newReadsFile.c_str(), "w");

  while (getline(readsFile, line1)) {
    getline(readsFile, line2);
    getline(readsFile, line3);
    getline(readsFile, line4);
    n++;
    if (duplicate_reads[line2]) {
      d++;
    } else {
      duplicate_reads[line2] = 1;
      fprintf(withoutDupReads, "%s\n", line1.c_str());  
      fprintf(withoutDupReads, "%s\n", line2.c_str());  
      fprintf(withoutDupReads, "%s\n", line3.c_str());  
      fprintf(withoutDupReads, "%s\n", line4.c_str());  
    }
  }
  cout << d << " duplicates in " << n;
}

/**
 * aligns read against the prebuilt index
 * indexPrefix - prefix name used by some aligners
 * readName - filename of the read that is going to be aligned
 * aligner - bwa, bowtie, bowtie2 or lisa
 * resultFile - filename for the result alignment
 */
void alignReadAgainstGenome(string refGenome, string indexPrefix, string readName, string aligner, string resultFile) {
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
    mapper += " -k 20 -D 15 -R 2 -N 0 -L 20 -i S,1,1.25 --gbar 4 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --score-min C,-14,0 -p 1";
    mapper += " > ";
    mapper += resultFile;
    cout << mapper << endl;
    system(mapper.c_str());  
  } else if (aligner == "lisa") {
    string lisa_mapper = "client solve ";
    lisa_mapper += refGenome;
    lisa_mapper += " ";
    lisa_mapper += indexPrefix;
    lisa_mapper += " ";
    lisa_mapper += readName;
    lisa_mapper += " ";
    lisa_mapper += resultFile;
    system(lisa_mapper.c_str());
  }
}

void saveResultCigFileMappingData() {
  ifstream samFile(result_cig_file.c_str());
  string line;

  while(getline(samFile, line)) {

    istringstream iss(line);
    string col;
    vector <string> cols;
    while(iss >> col) {
      cols.push_back(col);
    }

    SEQ_MAP_SOL[cols[0]] = atoi(cols[2].c_str());
  }
}

/*
 * checks input arguments passed from main
 */
int checkArguments(int argc, char *argv[]) {

  if (argc < 5) {
    cout << "Usage: rna-seq <aligner> <genome.fa> <index_dir> <reads.fq> <result_file>" << endl;
    cout << endl;
    cout << "<aligner>      - 'bowtie', 'bowtie2', 'bwa' or 'lisa'" << endl;
    cout << "<genome.fa>    - fasta file containing genome" << endl;
    cout << "<index_dir>    - prefix used for aligners for storing temporary files" << endl;
    cout << "<reads.fq>     - fastq file with reads to be aligned" << endl;
    cout << "<result_file>  - file in sam format containing overall alignments" << endl;
    cout << endl;
    cout << "options:" << endl;
    cout << "--cov-limit [coverage-limit-number]    - potential exons with coverage limit lower than specified are discarded" << endl;
    cout << "--skip-building-index                  - when running on same genome index again no need to wait for index building" << endl;
    cout << "--res-cig [res.cig]                    - file in cig format for validating results" << endl;
    cout << "--sliding-window-size/-sws [size]      - size of the granulation for cutting initially unmatched reads into smaller chunks - used for spliced matching" << endl;
    cout << "--join-exons-threshold/-jet [size]     - if two exons are too close join them (default: 100)" << endl;
    cout << "--extend-exons-threshold/-eet [size]   - how much to extend exons on right and left side for finding reads that were spliced (default:50)" << endl;
    cout << "--get-histogram-data                   - creates a hist.in file showing coverage data before spliced alignment" << endl;
    cout << "--num-threads [num]                    - number of threads for spliced alignment" << endl;
    cout << endl;
    return -1;
  }

  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "--skip-building-index") == 0) {
      skip_build_index = 1;
    }
    if (strcmp(argv[i], "--res-cig") == 0) {
      if ( i + 1 < argc ) {
        hasResultCigFile = 1;
        result_cig_file = argv[i + 1];
        saveResultCigFileMappingData();
      }
    }
    if (strcmp(argv[i], "--sliding-window-size") == 0 || strcmp(argv[i], "-sws") == 0) {
      SLIDING_WINDOW_SIZE = atoi(argv[i + 1]);
    }
    
    if (strcmp(argv[i], "--join-exons-threshold") == 0 || strcmp(argv[i], "-jet") == 0) {
      join_threshold = atoi(argv[i + 1]);
    }
    
    if (strcmp(argv[i], "--extend-exons-threshold") == 0 || strcmp(argv[i], "-eet") == 0) {
      extend_threshold = atoi(argv[i + 1]);
    }
    
    if (strcmp(argv[i], "--get-histogram-data") == 0) {
      GET_DATA_FOR_HISTOGRAM = 1;
    }
    
    if (strcmp(argv[i], "--num-threads") == 0) {
      num_threads = atoi(argv[i + 1]);
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

  int mapped_correctly = 0;

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
      if (hasResultCigFile) {
        if (SEQ_MAP_SOL[cols[0]] == atoi(cols[3].c_str())) {
          mapped_correctly++;
        }
      }
      matched_reads.push_back(Read(cols[0], cols[9], atoi(cols[3].c_str())));
      matched_reads[matched_reads.size() - 1].chrName = cols[2];
      fprintf(final_result, "%s %d %d %s\n", cols[2].c_str(), atoi(cols[3].c_str()) - 1, atoi(cols[3].c_str()) + 100 - 1, cols[0].c_str());
    }
  }

  cout << "number of reads: " << num_reads << endl;
  cout << "number of mapped reads: " << num_mapped << " (" << (double)num_mapped/num_reads*100 << "%)" << endl;
  cout << "number of unmapped reads: " << num_unmapped << " (" << (double)num_unmapped/num_reads*100 << "%)" << endl;
  
  if (hasResultCigFile) {
    cout << "of those " << num_mapped << " (" << (double)num_mapped/num_reads*100 << "%)" << " -> CORRECTLY MAPPED " << mapped_correctly << " (" << (double)mapped_correctly/num_mapped*100 << "%)" << endl;
  }
}

#define HIST_SIZE 160000000
int hist[HIST_SIZE] = {0};

void prepareDataForHistogramDrawing() {
  cout << "asi3" << endl;
  cout << "JURO " << matched_reads[0].map_position << " " <<  matched_reads[matched_reads.size() - 1].map_position << endl;
  int in = 0, out = 0;
  for (int i = 0; i < matched_reads.size(); i++) {
    if(matched_reads[i].map_position >= 3000000 && matched_reads[i].map_position <= 3500000) in++;
    else out++;
    for (int j = matched_reads[i].map_position; j < MIN(matched_reads[i].map_position + matched_reads[i].length, HIST_SIZE -1); j++) {
      hist[j]++;  
    }  
  }

  FILE *hist_f = fopen("hist.in", "w");
  for (int i = 0; i < HIST_SIZE; i++)
    //if (i >= 3000000 && i <= 3500000 ) {
      fprintf(hist_f, "%d\n", hist[i]);
    //}
  cout << in << " vs " << out << endl;

  cout << matched_reads[0].map_position << " " << matched_reads[matched_reads.size() - 1].map_position << endl;
}

void calculatePossibleExonRegions() {

  int total = matched_reads.size() + unmatched_reads.size();

  sort(matched_reads.begin(), matched_reads.end());

  clock_t begin = clock();
  if (GET_DATA_FOR_HISTOGRAM) {
    memset(hist, 0, sizeof(hist));
    prepareDataForHistogramDrawing();  
  }
  clock_t end = clock();

  cout << "ELAPSED in insertIntoExonRegions: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

  begin = clock();

  for (int i = 0; i < matched_reads.size(); i++) {
    // TODO REMOVE
    //if (matched_reads[i].map_position < 30000000 || matched_reads[i].map_position > 35000000)
    //  continue;
    insertIntoExonRegions(matched_reads[i]);
    // TODO: REMOVE THIS
    // if (exon_regions.size() > 30) break;
  }

  end = clock();

  cout << "ELAPSED in insertIntoExonRegions: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

  begin = clock();
  int all = 0;
  for (int i = 0; i < exon_regions.size(); i++) {
    cout << calcCoverage(exon_regions[i]) << endl;
    cout << "EXON: " << exon_regions[i].left << " " << exon_regions[i].right << endl;
    printExon(exon_regions[i]);
    if (calcCoverage(exon_regions[i]) >= coverage_limit) {
      all += exon_regions[i].c;
      possible_exon_regions.push_back(exon_regions[i]);
      //cout << "EXON: " << exon_regions[i].left << " " << exon_regions[i].right << endl;
      // TODO: REMOVE THIS
      //if (possible_exon_regions.size() > 30) break;
    }
  }
  cout << "MAX " << covm << endl;
  end = clock();

  cout << "ELAPSED in calcCov: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

  begin = clock();
  for (int i = 0; i < matched_reads.size(); i++) {
    int l = matched_reads[i].map_position;
    int r = l + matched_reads[i].length - 1;

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
  end = clock();

  cout << "ELAPSED in remainging unmatched exons: " << (double)(end - begin) / CLOCKS_PER_SEC << endl;

  cout << "AFTER INITIAL EXON COVERAGE SEARCH" << endl;
  cout << "MAPPED ON EXONS: " << total - unmatched_reads.size() << " (";
  cout <<  (total - unmatched_reads.size()) / (double)total * 100 << ")" << endl;
  cout << "UNMAPPED READS: " << unmatched_reads.size() << " (";
  cout <<  unmatched_reads.size() / (double)total * 100 << ")" << endl;
}

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

  cout << "number of possible exon regions: " << possible_exon_regions.size() << endl;
  cout << "original genome length: " << whole_genome.length() << endl;
  cout << "exon genome length: " << exon_genome.length() << " (";
  cout << exon_genome.length() / (double)whole_genome.length() * 100 << "% of original size)" << endl;
}

vector <int> exon_positions;

void addExonPositions(int start, int end) {
  for (int i = start; i <= end; i++) {
    exon_positions.push_back(i) ; 
  }  
}

void connectExonRegions2() {
  exon_genome = "";
  spliceSites.clear();

  cout << "w2 " << whole_genome.size() << endl;

  int joined = 0;

  exon_genome +=
    whole_genome.substr(
        possible_exon_regions[0].left,
        possible_exon_regions[0].right - possible_exon_regions[0].left + 1);

  addExonPositions(possible_exon_regions[0].left, possible_exon_regions[0].right);

  for (int i = 1; i < possible_exon_regions.size(); i++) {
    cout << possible_exon_regions[i].left << " " << possible_exon_regions[i].right << ", LENGTH= " << possible_exon_regions[i].right - possible_exon_regions[i].left << endl;
    /*cout << 
      whole_genome.substr(
      possible_exon_regions[i].left,
      possible_exon_regions[i].right - possible_exon_regions[i].left + 1) << endl;
      cout << endl;*/

    if ((possible_exon_regions[i].left - 1) - (possible_exon_regions[i - 1].right+1) + 1 <= join_threshold) {
      joined++;
      exon_genome += whole_genome.substr(
        possible_exon_regions[i - 1].right + 1,
        (possible_exon_regions[i].left - 1) - (possible_exon_regions[i - 1].right + 1) + 1
      );
      
      if (possible_exon_regions[i - 1].right >= 18648889 && possible_exon_regions[i - 1].right <= 18649417) {
        cout << "SE GRUPIRAT spojeni" << endl;
        cout << possible_exon_regions[i - 1].right << " " << possible_exon_regions[i].left << endl;
        cout << whole_genome.substr(
        possible_exon_regions[i - 1].right + 1,
        (possible_exon_regions[i].left - 1) - (possible_exon_regions[i - 1].right + 1) + 1
      ) << endl;
 
        cout << "SE GRUPIRAT" << endl; 
      }
      
      addExonPositions(possible_exon_regions[i - 1].right + 1, possible_exon_regions[i].left - 1);
    } else {
      exon_genome += whole_genome.substr(
        possible_exon_regions[i - 1].right + 1,
        extend_threshold
      );
      
      addExonPositions(possible_exon_regions[i - 1].right + 1, possible_exon_regions[i - 1].right + extend_threshold);

      exon_genome += whole_genome.substr(
        possible_exon_regions[i].left - 1 - extend_threshold,
        extend_threshold
      );
      
      addExonPositions(possible_exon_regions[i].left - extend_threshold, possible_exon_regions[i].left - 1);
    }

    exon_genome +=
      whole_genome.substr(
          possible_exon_regions[i].left,
          possible_exon_regions[i].right - possible_exon_regions[i].left + 1);
 
    if (possible_exon_regions[i - 1].right >= 18648889 && possible_exon_regions[i - 1].right <= 18649417) {
      cout << "SE GRUPIRAT" << endl;
      cout << whole_genome.substr(
                possible_exon_regions[i].left,
                          possible_exon_regions[i].right - possible_exon_regions[i].left + 1) << endl;
      cout << "SE GRUPIRAT" << endl; 
    }
  
    addExonPositions(possible_exon_regions[i].left, possible_exon_regions[i].right);

    spliceSites.push_back(exon_genome.size());
  }

  if (spliceSites.size()) spliceSites.pop_back();

  cout << "number of possible exon regions: " << possible_exon_regions.size() << endl;
  cout << "number of joined exons because too close: " << joined << endl;
  cout << "original genome length: " << whole_genome.length() << endl;
  cout << "exon genome length: " << exon_genome.length() << " (";
  cout << exon_genome.length() / (double)whole_genome.length() * 100 << "% of original size)" << endl;

  for (int i = 0; i < exon_genome.size(); i++) {
    exon_genome[i] = toupper(exon_genome[i]);  
  }

}

void tryConnectingExons() {

  int searchLeftSpan = 100, searchRightSpan = 100;

  vector < pair <string, string> > known_intron_signals;
  
  known_intron_signals.push_back(make_pair("GT", "AG"));
  //known_intron_signals.push_back(make_pair("GC", "AG"));
  //known_intron_signals.push_back(make_pair("AT", "AC"));

  for (int i = 0; i < 30; i++) {
    
    cout << possible_exon_regions[i].left << " " << possible_exon_regions[i].right << endl;
    cout << possible_exon_regions[i + 1].left << " " << possible_exon_regions[i + 1].right << endl;
    
    for (int j = 0; j < known_intron_signals.size(); j++) {

      int possibleStartForLeftExon = MAX(
          possible_exon_regions[i].right - searchLeftSpan,
          possible_exon_regions[i].left + searchRightSpan
        );

      int possibleEndForLeftExon = MIN(
          possible_exon_regions[i].right + searchRightSpan,
          possible_exon_regions[i + 1].left - searchLeftSpan
        );

      vector <int> signalStart, signalEnd;

      for (int x = possibleStartForLeftExon; x < possibleEndForLeftExon; x++) {
        if (whole_genome[x] == known_intron_signals[j].first[0]) {
          if (whole_genome[x + 1] == known_intron_signals[j].first[1]) {
            signalStart.push_back(x);
          }
        }
      }

      int possibleEndForRightExon = MIN(
          possible_exon_regions[i + 1].left + searchLeftSpan,
          possible_exon_regions[i + 1].right - searchRightSpan
        );

      int possibleStartForRightExon = MAX(
          possible_exon_regions[i + 1].left - searchRightSpan,
          possible_exon_regions[i].right + searchLeftSpan
        );
      for (int x = possibleStartForRightExon; x < possibleEndForRightExon; x++) {
        if (whole_genome[x] == known_intron_signals[j].second[0]) {
          if (whole_genome[x + 1] == known_intron_signals[j].second[1]) {
            signalEnd.push_back(x);
          }
        }
      }

      cout << "signal start: " << endl;
      for (int k = 0; k < signalStart.size(); k++) {
        cout << signalStart[k] << " " << whole_genome.substr(signalStart[k] - 20, 22) << endl;  
      }

      cout << "signal end: " << endl;
      for (int k = 0; k < signalEnd.size(); k++) {
        cout << signalEnd[k] << " " << whole_genome.substr(signalEnd[k] - 20, 22) << endl;  
      }

    }

  }    

}

//#define INDEX_SIZE (1LL<<32)
//vector <int> exon_index[INDEX_SIZE];

map <string, vector<int> > exon_index;


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

int junc_first, junc_last;

int printIndexValue(string hashValue) {
  cout << "(";
  vector <int> v = exon_index[hashValue];
  for (int i = 0; i < v.size(); i++ ) {
    cout << v[i] << " -> " << exon_positions[v[i]] << "," ;
    return exon_positions[v[i]];
    if (!junc_first) {
      junc_first = exon_positions[v[i]];  
    } else {
      junc_last = exon_positions[v[i]] + SLIDING_WINDOW_SIZE;  
    }
  }
  cout << ")";
  return -1;
}

int DP[1024], prev[1024];

/*
 * calculates the longest increasing subsequence
 */
void calcLIS(const vector<int> &elements) {
  
  assert(elements.size() < 1024);
  // reads maximum length can not exceed 1024

  int maxLength = 1, bestEnd = 0;
  DP[0] = 1;
  prev[0] = -1;

  for (int i = 1; i < elements.size(); i++) {
    DP[i] = 1;
    prev[i] = -1;
    
    for (int j = i - 1; j >= 0; j--) {
      if (DP[j] + 1 > DP[i] && elements[j] < elements[i]) {
        DP[i] = DP[j] + 1;
        prev[i] = j;  
      }
    }

    if (DP[i] > maxLength) {
      bestEnd = i;
      maxLength = DP[i];  
    }
    
  }

  cout << "MAX LENGTH: " << maxLength << endl;

  vector <int> lis;

  while (prev[bestEnd] != -1) {
    cout << elements[bestEnd] << " ";
    lis.push_back(elements[bestEnd]);
    bestEnd = prev[bestEnd]; 
  }
  cout << endl;

  junc_last = lis[0] + SLIDING_WINDOW_SIZE;
  junc_first = lis[lis.size() - 1];

}

/*
 * Matches the rest of reads against the exon index
 */
void unmatchedReadsAgainstIndex() {

  int couldBe = 0;
  for (int i = 0; i < unmatched_reads.size(); i++) {

    string read = unmatched_reads[i].value;
    for (int j = 0; j < read.size(); j++) {
      read[j] = toupper(read[j]);
    }

    int matched = 0;
    for (int j = 0; j < read.size() - SLIDING_WINDOW_SIZE + 1; j++) {
      matched += (exon_index[read.substr(j, SLIDING_WINDOW_SIZE)].size() > 0);
    }

    if (matched < 20) continue;

    vector <int> elements;

    junc_first = 0;
    junc_last = 0;
    
    couldBe++;
    cout << read << endl;
    cout << unmatched_reads[i].id << endl;
    for (int j = 0; j < read.size() - SLIDING_WINDOW_SIZE + 1; j++) {
      //printIndexValue(hash(read.substr(j, SLIDING_WINDOW_SIZE)));
      int el = printIndexValue(read.substr(j, SLIDING_WINDOW_SIZE));
      elements.push_back(el);
    }
    cout << endl;

    calcLIS(elements);

    string chrName = "chr1";
    fprintf(final_result, "%s %d %d %s\n", chrName.c_str(), junc_first, junc_last, unmatched_reads[i].id.c_str());
  
  }
  
  cout << "OUT OF " << unmatched_reads.size() << " - " << couldBe << " " << (double)couldBe / unmatched_reads.size() * 100 << "%" << endl;
}

/*
 * Main method of rna-seq
 */
int main(int argc, char *argv[]) {

  clock_t begin_overall = clock();

  if (checkArguments(argc, argv) == -1) return -1;

  string alignerName = argv[1];
  string refGenomeName = argv[2];
  string indexPrefix = argv[3];
  string readInputName = argv[4];
  string resultFile = argv[5];

  final_result = fopen("result.bam", "w");

  if (!skip_build_index) {
    clock_t begin = clock();
    cout << "building index for genome alignment..." << endl;
    buildIndexForGenomeAlignment(refGenomeName, indexPrefix, alignerName);
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "ELAPSED: building genome index: " <<  elapsed_secs << endl;
  }

  //removeDuplicates(readInputName, readInputName + "_");

  if (alignerName == "bwa") {
    indexPrefix = refGenomeName;  
  }

  clock_t begin = clock();
  cout << "aligning reads against genome..." << endl;
  alignReadAgainstGenome(refGenomeName, indexPrefix, readInputName, alignerName, resultFile);
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "ELAPSED: genome alignment " << elapsed_secs << endl;

  begin = clock();
  cout << "getting data from initial alignment..." << endl;
  extractReadAlignmentDataFromSamFile(resultFile);
  end = clock();
  elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "ELAPSED: extracting aln data " << elapsed_secs << endl;

  begin = clock();
  extractGenomeFromFile(refGenomeName);
  end = clock();
  elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "ELAPSED: extracting genome " << elapsed_secs << endl;
 
  begin = clock(); 
  cout << "calculating possible exon regions..." << endl;
  calculatePossibleExonRegions();
  end = clock();
  elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "ELAPSED: calculating possible exons " << elapsed_secs << endl;


  begin = clock();
  //connectExonRegions();
  cout << "connecting exon regions" << endl;
  connectExonRegions2();
  end = clock();
  elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "ELAPSED: connecting exon regions " << elapsed_secs << endl;

  //tryConnectingExons();

  //return -3; // testing

  begin = clock();
  cout << "building index for spliced alingment" << endl;
  createIndex();
  end = clock();
  elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "ELAPSED: building index for spliced aln " << elapsed_secs << endl;

  begin = clock();
  cout << "spliced alingment began" << endl;
  unmatchedReadsAgainstIndex();
  end = clock();
  elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "ELAPSED: spliced alignment " <<  elapsed_secs << endl;

  clock_t end_overall = clock();

  cout << "ELAPSED OVERALL: " << (double)(end_overall - begin_overall) / CLOCKS_PER_SEC << endl;

  return -3;

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
  fprintf(outGenome, ">chr20\n");
  for (int i = 0; i < exon_genome.size(); i+=80) {
    fprintf(outGenome,
        "%s\n",
        exon_genome.substr(i, MIN(80, exon_genome.size() - i)).c_str());
  }

  buildIndexForGenomeAlignment("tmp.fa", "tmp", alignerName);
  alignReadAgainstGenome("tmp", "tmp", "unmatched.fastq", alignerName, "tmp.sam");

  return 0;  
}
