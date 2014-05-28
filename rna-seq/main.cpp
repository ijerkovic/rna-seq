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

using std::vector;
using std::string;
using std::pair;
using std::cout;
using std::ifstream;
using std::endl;
using std::make_pair;

#define MAX_READ_LENGTH (1 << 14)

#define MAX(A, B) (((A) > (B)) ? (A) : (B))
#define MIN(A, B) (((A) < (B)) ? (A) : (B))

#define PRINT_INTRONS_GTF 0

std::vector<std::string> &split(
  const std::string &s,
  char delim,
  std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      elems.push_back(item);
    }
    return elems;
  }

std::vector<std::string> split(
  const std::string &s,
  char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
  }

struct Read {
  string id;
  string value;
  Read() {}
  Read(string id_, string value_) {
    id = id_; value = value_;
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

vector <Read> all_reads;
vector <string> all_refs;
vector <string> unmatched_reads_refs;
vector < Exon > exon_regions;
vector <pair <int, int > > gtf_exon_regions;
vector <int> start_positions;
vector < Exon > possible_exon_regions;
vector < pair<string, string> > known_intron_signals;

// TODO(ijerkovic): for now only same read_length supported of all reads
int read_length;

// TODO(ijerkovic): should be able to pass in from terminal in future
int coverage_limit;

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

int main(int argc, char *argv[]) {
  if (argc < 5) {
    cout << "Usage: rna-seq <genome.fa> <index_dir> <reads.fq> <result_file> <optional-gtf-file> --cov-limit [coverage-limit-number]" << endl;
    return -1;
  }

  freopen ("myfile.txt", "w", stdout);

  clock_t START_TIME = clock();

  coverage_limit = 20;
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "--cov-limit")) {
      if (i + 1 < argc) {
        coverage_limit = atoi(argv[i + 1]);
      }
    }
  }
  if (!coverage_limit) {
    coverage_limit = 20;  
  }
  cout << "COV_LIMIT= " << coverage_limit << endl;

  // TODO(ijerkovic): should be able to pass in as an option in future

  cout << endl;
  cout << "GENOME ALIGNMENT STARTED" << endl;

  cout << endl;
  cout << "BUILDING LISSA INDEX" << endl;
  cout << endl;
  string lisa_build_index = "client index ";
  lisa_build_index += argv[1];
  lisa_build_index += " ";
  lisa_build_index += argv[2];
  lisa_build_index += "\n";
  system(lisa_build_index.c_str());

  cout << endl;
  cout << "LISSA MAPPING READS AGAINST INDEX" << endl;
  cout << endl;
  string lisa_mapper = "client solve ";
  lisa_mapper += argv[1];
  lisa_mapper += " ";
  lisa_mapper += argv[2];
  lisa_mapper += " ";
  lisa_mapper += argv[3];
  lisa_mapper += " ";
  lisa_mapper += argv[4];
  system(lisa_mapper.c_str());

  ifstream lissa_results(argv[4]);
  string line;
  int num_all_reads = 0;
  int num_matched_reads = 0;

  while (std::getline(lissa_results, line)) {
    num_all_reads++;

    std::vector<std::string> stats = split(line, ';');

    std::vector<std::string> firstCol = split(stats[0], ',');

    int num_matches = atoi(firstCol[firstCol.size() - 1].c_str());
    
    
    //num_matched_reads += (num_matches > 0);
    
    

    if (!num_matches) {
      unmatched_reads_refs.push_back(firstCol[0]);
    } else {
      for (int i = 1; i < stats.size(); i++) {
        std::vector<std::string> tmp = split(stats[i], ',');

      
        // CONSIDER ONLY FULL MATCHES
        if(atoi(tmp[1].c_str()) >= 43 ) {
          //unmatched_reads_refs.push_back(firstCol[0]);
          start_positions.push_back(atoi(tmp[2].c_str()));
          all_refs.push_back(firstCol[0]);
          num_matched_reads++;
          break;  
        }


        // TODO(ijerkovic): for now just take the match with best score
        //start_positions.push_back(atoi(tmp[2].c_str()));
        //all_refs.push_back(firstCol[0]);
        //break;
      }
    }
  }

  cout << endl;
  cout << "GENOME ALIGNMENT FINISHED" << endl;
  cout << "READS MATCHED: " << num_matched_reads << "/" << num_all_reads
       << " (" << (double)num_matched_reads/(double)num_all_reads << ")"
       << endl;
  cout << endl;

  cout << "SPLICED ALIGNMENT STARTED" << endl;
  cout << endl;

  cout << "getting unmatched reads" << endl;
  cout << endl;

  read_length = 0;

  FILE *all_reads_file = fopen(argv[3], "r");
  char buff1[MAX_READ_LENGTH];
  char buff2[MAX_READ_LENGTH];
  while (fgets(buff1, MAX_READ_LENGTH, all_reads_file)) {
    fgets(buff2, MAX_READ_LENGTH, all_reads_file);

    string read_id = buff1;
    read_id.resize(read_id.size() - 1);

    string read_val = buff2;
    read_val.resize(read_val.size() - 1);

    if (read_val.size() != read_length && read_length) {
      cout << "ERROR: all read lengths should be same in one fastaq file!" << endl;
      return -1;
    }
    read_length = read_val.size();

    all_reads.push_back(Read(read_id, read_val));

    fgets(buff2, MAX_READ_LENGTH, all_reads_file);
    fgets(buff2, MAX_READ_LENGTH, all_reads_file);
  }

  cout << "READ LENGTH= " << read_length << endl;

  FILE *unmatched_reads_file = fopen("unmatched.fq", "w");
  int j = 0;
  for (int i = 0; i < unmatched_reads_refs.size(); i++) {
    for (; j < all_reads.size(); j++) {
      if (all_reads[j].id == unmatched_reads_refs[i]) {
        fprintf(unmatched_reads_file, "%s\n", all_reads[j].id.c_str());
        fprintf(unmatched_reads_file, "%s\n", all_reads[j].value.c_str());
        string tmp = string(all_reads[j].value.size(), 'E');
        fprintf(unmatched_reads_file, "+\n%s\n", tmp.c_str());
        break;
      }
    }
  }

  cout << "tu sam" << endl;

  // GTF file is provided - exon regions can be reconstructed
  if(argc >= 6) {

    cout << "tu sam2" << endl;
    if (strcmp(argv[5], "--cov-limit") && !atoi(argv[5])) {
      // TODO(ijerkovic): still not being considered into calculations

      cout << "extracting exon information from gtf annotation" << endl;

      ifstream gtf_annotation(argv[5]);
      while (std::getline(gtf_annotation, line)) {
        std::vector<std::string> data1 = split(line, ';');
        replace(data1[0].begin(), data1[0].end(), '\t', ' ');
        std::vector<std::string> data = split(data1[0], ' ');

        if (data[2] == "exon") {
          if (atoi(data[3].c_str()) >= 3000000 && atoi(data[4].c_str()) <= 3500000) {
            gtf_exon_regions.push_back(
              make_pair(atoi(data[3].c_str()), atoi(data[4].c_str())));
          }
    
        }
      }
    }
  }

  cout << "GTF-ANNOTATIONS_AMOUNT= " << gtf_exon_regions.size() << endl;

  std::map <string, int> MAPA;

  // PRINT OUT INTRONS
  if( PRINT_INTRONS_GTF ) {
    FILE *inGenome = fopen(argv[1], "r");
    char buff[256];
    string whole_genome = "";
    while (fgets(buff, 255, inGenome)) {
      // ignore lines that start with genome desc, they start with '>'
      if (buff[0] != '>') {
        string tmp = buff;
        tmp.resize(tmp.size() - 1);  // remove endl
        whole_genome += tmp;
      }
    }
    for(int i = 0; i < gtf_exon_regions.size(); i++) {
      cout << "INTRON " << i + 1 << endl;
      cout << whole_genome.substr(gtf_exon_regions[i].first, gtf_exon_regions[i].second - gtf_exon_regions[i].first + 1) << endl;
      string intron_junction = whole_genome.substr(gtf_exon_regions[i].first, 2) + "-" + whole_genome.substr(gtf_exon_regions[i].second - 1, 2);
      cout << intron_junction << endl;
      MAPA[intron_junction]++;
      cout << endl;  
    }
    typedef std::map<std::string, int>::iterator it_type;
    for(it_type iterator = MAPA.begin(); iterator != MAPA.end(); iterator++) {
      cout << iterator->first << " " << iterator->second << endl; 
    }
  }
  //

  cout << "reconstructing exon regions" << endl;
  cout << endl;

  sort(start_positions.begin(), start_positions.end());

  for (int i = 0; i < start_positions.size(); i++) {
    insertIntoExonRegions(start_positions[i]);
  }

  int all = 0;
  for (int i = 0; i < exon_regions.size(); i++) {
    if (calcCoverage(exon_regions[i]) > coverage_limit) {
      all += exon_regions[i].c;
      possible_exon_regions.push_back(exon_regions[i]);
    }
  }

  unmatched_reads_refs.clear();

  for (int i = 0; i < start_positions.size(); i++) {
    int l = start_positions[i];
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
      unmatched_reads_refs.push_back(all_refs[i]);
    }
  }

  cout << (num_all_reads - unmatched_reads_refs.size()) / (double)num_all_reads <<
          "% reads matched to exon regions" << endl;

  FILE *summaryFile = fopen("summary.txt", "w");
  fprintf(summaryFile,
          "%f %% reads matched to exon regions before spliced alignment\n",
          (num_all_reads - unmatched_reads_refs.size()) / (double)num_all_reads);

  // reads that dont match exon regions
  j = 0;
  for (int i = 0; i < unmatched_reads_refs.size(); i++) {
    for (; j < all_reads.size(); j++) {
      if (all_reads[j].id == unmatched_reads_refs[i]) {
        // ALL READS!!!!!!!!!! //cout << all_reads[j].value << endl;
        fprintf(unmatched_reads_file, "%s\n", all_reads[j].id.c_str());
        fprintf(unmatched_reads_file, "%s\n", all_reads[j].value.c_str());
        string tmp = string(all_reads[j].value.size(), 'E');
        fprintf(unmatched_reads_file, "+\n%s\n", tmp.c_str());
        break;
      }
    }
  }

  FILE *inGenome = fopen(argv[1], "r");
  char buff[256];
  string whole_genome = "";
  while (fgets(buff, 255, inGenome)) {
    // ignore lines that start with genome desc, they start with '>'
    if (buff[0] != '>') {
      string tmp = buff;
      tmp.resize(tmp.size() - 1);  // remove endl
      whole_genome += tmp;
    }
  }

  known_intron_signals.push_back(make_pair("GT", "AG"));
  known_intron_signals.push_back(make_pair("GC", "AG"));
  known_intron_signals.push_back(make_pair("AT", "AC"));

  // FINDING INTRONS VIA CANONIC KNOWN SIGNALS
  for (int i = 1; i < possible_exon_regions.size(); i++) {
    
    int bestDiff = -1;
    int newL, newR;
    for (int j = 0; j < known_intron_signals.size(); j++) {
      
      int found1 = -1;
      for (int x = possible_exon_regions[i - 1].right - 20; x < possible_exon_regions[i].left; x++) {
        if (whole_genome[x] == known_intron_signals[j].first[0]) {
          if (whole_genome[x + 1] == known_intron_signals[j].first[1]) {
            found1 = x;
          }
        }
        if (x > possible_exon_regions[i - 1].right + 20 && found1 != -1) {
          break;  
        }
      }
      
      int found2 = -1;
      for (int x = possible_exon_regions[i].left + 20; x > possible_exon_regions[i - 1].right; x--) {
        if (whole_genome[x] == known_intron_signals[j].second[0]) {
          if (whole_genome[x + 1] == known_intron_signals[j].second[1]) {
            found2 = x;
          }
        }
        if (x < possible_exon_regions[i - 1].left - 20 && found2 != -1) {
          break;
        }
      }

      if (found1 < found2)
      cout << "PAROVI "<< found1 << " " << found2 << endl;

      if (found1 != -1 && found2 != -1 && found1 < found2) {
        if (abs(possible_exon_regions[i - 1].right - found1) + abs(possible_exon_regions[i].left - found2) < bestDiff) {
          bestDiff = abs(possible_exon_regions[i - 1].right - found1) + abs(possible_exon_regions[i].left - found2);
          newL = found1;
          newR = found2;
        } else if (bestDiff == -1) {
          bestDiff = abs(possible_exon_regions[i - 1].right - found1) + abs(possible_exon_regions[i].left - found2);
          newL = found1;
          newR = found2;  
        }
      }
    }

    // CANONIC INTRONS FOUND
    if (bestDiff != -1) {
      cout << "POSSIBLE EXON DELTA= " << bestDiff << endl;
      cout << "NEW LEFT= " << newL << "   NEW RIGHT= " << newR << endl;
      possible_exon_regions[i - 1].right = newL;
      possible_exon_regions[i].left = newR;
    }
  }

  vector <int> spliceSites;

  string exon_genome = "";
  for (int i = 0; i < possible_exon_regions.size(); i++) {
    cout << possible_exon_regions[i].left << " " << possible_exon_regions[i].right << ", LENGTH= " << possible_exon_regions[i].right - possible_exon_regions[i].left << endl;
    cout << 
      whole_genome.substr(
        possible_exon_regions[i].left,
        possible_exon_regions[i].right - possible_exon_regions[i].left + 1) << endl;
        cout << endl;

    exon_genome +=
      whole_genome.substr(
        possible_exon_regions[i].left,
        possible_exon_regions[i].right - possible_exon_regions[i].left + 1);

    spliceSites.push_back(exon_genome.size());
  }

  cout << "POSSIBLE SPLICE JUNCTIONS FOUND: " << spliceSites.size() << endl;
  if (spliceSites.size()) spliceSites.pop_back();

  FILE *outGenome = fopen("tmp.fa", "w");
  fprintf(outGenome, ">chr19\n");
  for (int i = 0; i < exon_genome.size(); i+=80) {
    fprintf(outGenome,
            "%s\n",
            exon_genome.substr(i, MIN(80, exon_genome.size() - i)).c_str());
  }


  cout << endl;
  cout << "ALIGN AGAINST POSSIBLE EXON REGIONS" << endl;

  system("mkdir tmp");
  fclose(outGenome);

  cout << endl;
  cout << "BUILDING LISSA INDEX" << endl;
  cout << endl;
  lisa_build_index = "client index ";
  lisa_build_index += "tmp.fa";
  lisa_build_index += " ";
  lisa_build_index += "tmp";
  lisa_build_index += "\n";
  system(lisa_build_index.c_str());

  fclose(unmatched_reads_file);

  cout << endl;
  cout << "LISSA MAPPING READS AGAINST INDEX" << endl;
  cout << endl;
  lisa_mapper = "client solve ";
  lisa_mapper += "tmp.fa";
  lisa_mapper += " ";
  lisa_mapper += "tmp";
  lisa_mapper += " ";
  lisa_mapper += "unmatched.fq";
  lisa_mapper += " ";
  lisa_mapper += "tmp-res";

  system(lisa_mapper.c_str());

  ifstream lissa_results2("tmp-res");
  int still_unmatched = 0;
  int all_of_unmatched = 0;
  while (std::getline(lissa_results2, line)) {
    all_of_unmatched++;

    std::vector<std::string> stats = split(line, ';');

    std::vector<std::string> firstCol = split(stats[0], ',');

    int num_matches = atoi(firstCol[firstCol.size() - 1].c_str());
    num_matched_reads += (num_matches > 0);

    if (!num_matches) {
      still_unmatched++;
    } else {
      for (int i = 1; i < stats.size(); i++) {
        std::vector<std::string> tmp = split(stats[i], ',');

        // for now just take the match with best score
        start_positions.push_back(atoi(tmp[2].c_str()));
        all_refs.push_back(firstCol[0]);

        int isOk = 10000000;
        for (int j = 0; j < spliceSites.size(); j++) {
          if (abs(spliceSites[j]-atoi(tmp[2].c_str())) < isOk) {
            isOk = abs(spliceSites[j]-atoi(tmp[2].c_str())) ;
          }
        }
        cout << isOk << endl;

        if (isOk > read_length) {
          still_unmatched++;  
        }

        break;
      }
    }
  }
        for (int j = 0; j < spliceSites.size(); j++) {
          cout << "SPLICE " << spliceSites[j] << endl;
        }

  system("rm -rf tmp");
  //system("rm tmp*");

  cout << (double)(num_all_reads - still_unmatched) / (double)num_all_reads <<
          "% after running against possible exon regions" << endl;
  fprintf(summaryFile, "%f %% reads matched after running on discovered exon regions\n",
          (num_all_reads - still_unmatched) / (double)num_all_reads);

  clock_t END_TIME = clock();
  cout << "Running time: " << (double)(END_TIME - START_TIME) / CLOCKS_PER_SEC << endl;

  cout << "NUM EXON REGIONS= " << possible_exon_regions.size() << endl;

  return 0;
}
