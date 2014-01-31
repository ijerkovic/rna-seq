#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <iostream>
using namespace std;

#define MAX_READ_LENGTH (1 << 14)

#define MAX(A, B) (((A) > (B)) ? (A) : (B))
#define MIN(A, B) (((A) < (B)) ? (A) : (B))

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

struct Read 
{
  string id;
  string value;
  Read() {}
  Read(string id_, string value_) {
    id = id_; value = value_;  
  }
};

struct Exon
{
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
    for(int i = 0; i < 1024; i++)
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


// TODO: for now only same read_length supported of all reads
int read_length;

// TODO: should be able to pass in from terminal in future
int coverage_limit;

void addCoverage(Exon &e, int start, int end) {
  //e.c += 20;
  //cout << start << " " << end << endl;
  e.c++;
  for(int i = start; i <= MIN(end,1023); i++) {
    e.coverage[i]++;
    //cout << e.coverage[i]; 
  }  
}

double calcCoverage(Exon &e) 
{
  double cov = 0;
  for(int i = 0; i < 1024; i++) {
    //if (e.coverage[i] > 10)
    cov += e.coverage[i];
  }
  //cout << e.c << " " << cov / (double)(e.right - e.left) << endl;
  return cov / (double)(e.right - e.left);
}

int newE = 0, addedE = 0;

void insertIntoExonRegions(int left) 
{
  int found = 0;
  for(int i = 0; i < exon_regions.size(); i++) {
    if(left >= exon_regions[i].left && left <= exon_regions[i].right) {
      exon_regions[i].right = MAX(left + read_length, exon_regions[i].right);  
      found = 1;
      addedE ++;
      addCoverage(exon_regions[i], left - exon_regions[i].start, left + read_length - exon_regions[i].start);
    }  
  }
  if(!found) {
    newE++;
    exon_regions.push_back(Exon(left, left + read_length));
    addCoverage(exon_regions[exon_regions.size() - 1], 0, read_length);
  }
}

int main(int argc, char *argv[]) 
{
  // client index genome.fa index_dir
  // client solve genome.fa index_dir reads.fq result_file
  if(argc != 5 && argc != 6) {
    cout << "Usage: rna-seq <genome.fa> <index_dir> <reads.fq> <result_file> <optional-gtf-file>" << endl;
    return -1;  
  }

  clock_t START_TIME = clock();

  coverage_limit = 10; // TODO: should be able to pass in as an option in future

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

  while(std::getline(lissa_results, line)) {

    num_all_reads++;

    std::vector<std::string> stats = split(line, ';');
    
    //for(int i = 0; i < stats.size(); i++) {
    //  cout << stats[ i ] << endl;  
    //}
    
    std::vector<std::string> firstCol = split(stats[0], ',');

    //for(int i = 0; i < firstCol.size(); i++) {
    //  cout << firstCol[ i ] << endl;  
    //}

    int num_matches = atoi(firstCol[firstCol.size() - 1].c_str());
    num_matched_reads += (num_matches > 0);

    if(!num_matches) {
      unmatched_reads_refs.push_back(firstCol[0]);
    } else {
      for(int i = 1; i < stats.size(); i++) {
        std::vector<std::string> tmp = split(stats[i], ',');
        //cout << atoi(tmp[4].c_str()) << endl; continue;
        //int isEndGene = atoi(tmp[4].c_str());
        
        //if(isEndGene == 0 || 1) {

          // TODO: for now just take the match with best score
          start_positions.push_back(atoi(tmp[2].c_str()));
          all_refs.push_back(firstCol[0]);
          break;
      }  
    }
  }

  cout << endl;
  cout << "GENOME ALIGNMENT FINISHED" << endl;
  cout << "READS MATCHED: " << num_matched_reads << "/" << num_all_reads << " (" << (double)num_matched_reads/(double)num_all_reads << ")" << endl;
  cout << endl;

  cout << "SPLICED ALIGNMENT STARTED" << endl;
  cout << endl;

  cout << "getting unmatched reads" << endl;
  cout << endl;

  read_length = 0;

  FILE *all_reads_file = fopen(argv[3], "r");
  char buff1[MAX_READ_LENGTH];
  char buff2[MAX_READ_LENGTH];
  while(fgets(buff1, MAX_READ_LENGTH, all_reads_file)) {
    fgets(buff2, MAX_READ_LENGTH, all_reads_file);
    
    string read_id = buff1;
    read_id.resize(read_id.size() - 1);
    
    string read_val = buff2;
    read_val.resize(read_val.size() - 1);
  
    if(read_val.size() != read_length && read_length) {
      cout << "ERROR: all read lengths should be same in one fastaq file!" << endl;
      return -1;  
    }
    read_length = read_val.size();

    all_reads.push_back(Read(read_id, read_val));
    
    fgets(buff2, MAX_READ_LENGTH, all_reads_file);
    fgets(buff2, MAX_READ_LENGTH, all_reads_file);
  }

  FILE *unmatched_reads_file = fopen("unmatched.fq", "w");
  int j = 0;
  for(int i = 0; i < unmatched_reads_refs.size(); i++) {
    for(; j < all_reads.size(); j++) {
      if(all_reads[j].id == unmatched_reads_refs[i]) {
        fprintf(unmatched_reads_file, "%s\n", all_reads[j].id.c_str());
        fprintf(unmatched_reads_file, "%s\n", all_reads[j].value.c_str());
        string tmp = string(all_reads[j].value.size(), 'E');
        fprintf(unmatched_reads_file, "+\n%s\n", tmp.c_str());
        break;  
      }
    }
  }


  ifstream gtf_annotation(argv[5]);

  // GTF file is provided - exon regions can be reconstructed
  if(argc == 6 ) {
    
    // TODO: still not being considered into calculations

    cout << "extracting exon information from gtf annotation" << endl;
    // chr19 unknown exon  60951 61894 . - . gene_id "WASH5P"; transcript_id "NR_033266"; gene_name "WASH5P"; tss_id "TSS28072";
      
    while(std::getline(gtf_annotation, line)) {

      std::vector<std::string> data1 = split(line, ';');
      replace(data1[0].begin(), data1[0].end(), '\t', ' ');
      std::vector<std::string> data = split(data1[0], ' ');
      
      if (data[2] == "exon") {
        gtf_exon_regions.push_back(make_pair(atoi(data[3].c_str()), atoi(data[4].c_str())));
      }
    
    }

  }

  cout << "reconstructing exon regions" << endl;
  cout << endl;

  sort(start_positions.begin(), start_positions.end());
  
  for(int i = 0; i < start_positions.size(); i++) {
    insertIntoExonRegions(start_positions[i]); 
  }

  int all = 0;
  for(int i = 0; i < exon_regions.size(); i++) {
    if(calcCoverage(exon_regions[i])>coverage_limit) {
      //cout << exon_regions[i].c << endl;
      all += exon_regions[i].c;
      //cout << exon_regions[i].left << " " << exon_regions[i].right << " " << calcCoverage(exon_regions[i]) << " " << exon_regions[i].right - exon_regions[i].left << endl;  
      possible_exon_regions.push_back(exon_regions[i]);
    }
  }
  
  unmatched_reads_refs.clear(); 

  for(int i = 0; i < start_positions.size(); i++) {
    
    int l = start_positions[i];
    int r = l + read_length - 1;
    
    int found = 0;
    // this should be done in binary search
    for(int j = 0; j < possible_exon_regions.size(); j++) {
      if(l >= possible_exon_regions[ j ].left && r <= possible_exon_regions[ j ].right) {
        found = 1;
        break;  
      }
    }
    
    if(!found) {
      unmatched_reads_refs.push_back(all_refs[i]);
    }
  }
  
  cout << (num_all_reads - unmatched_reads_refs.size()) / (double)num_all_reads << "% reads matched to exon regions" << endl;

  FILE *summaryFile = fopen("summary.txt", "w");
  fprintf(summaryFile, "%f %% reads matched to exon regions before spliced alignment\n", (num_all_reads - unmatched_reads_refs.size()) / (double)num_all_reads );

  // reads that dont match exon regions
  j = 0;
  for(int i = 0; i < unmatched_reads_refs.size(); i++) {
    for(; j < all_reads.size(); j++) {
      if(all_reads[j].id == unmatched_reads_refs[i]) {
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
  while(fgets(buff, 255, inGenome)) {
    // ignore lines that start with genome desc, they start with '>'
    if(buff[0] != '>') {
      string tmp = buff;
      tmp.resize(tmp.size() - 1); // remove endl
      whole_genome += tmp;
    }
  }

  string exon_genome = "";
  for(int i = 0; i < possible_exon_regions.size(); i++) {
    exon_genome += whole_genome.substr(possible_exon_regions[i].left, possible_exon_regions[i].right - possible_exon_regions[i].left + 1);
  }

  FILE *outGenome = fopen("tmp.fa", "w");
  fprintf(outGenome, ">chr19\n");
  for(int i = 0; i < exon_genome.size(); i+=80) {
    fprintf(outGenome, "%s\n", exon_genome.substr(i, MIN(80, exon_genome.size() - i)).c_str());
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
  
  //cout << lisa_mapper << endl;

  system(lisa_mapper.c_str());
 
   
  ifstream lissa_results2("tmp-res");
  int still_unmatched = 0;
  int all_of_unmatched = 0;
  while(std::getline(lissa_results2, line)) {

    all_of_unmatched++;

    std::vector<std::string> stats = split(line, ';');
    
    //for(int i = 0; i < stats.size(); i++) {
    //  cout << stats[ i ] << endl;  
    //}
    
    std::vector<std::string> firstCol = split(stats[0], ',');

    //for(int i = 0; i < firstCol.size(); i++) {
    //  cout << firstCol[ i ] << endl;  
    //}

    int num_matches = atoi(firstCol[firstCol.size() - 1].c_str());
    num_matched_reads += (num_matches > 0);

    if(!num_matches) {
      //unmatched_reads_refs.push_back(firstCol[0]);
      still_unmatched++;
    } else {
      for(int i = 1; i < stats.size(); i++) {
        std::vector<std::string> tmp = split(stats[i], ',');
        //cout << atoi(tmp[4].c_str()) << endl; continue;
        //int isEndGene = atoi(tmp[4].c_str());
        
        //if(isEndGene == 0 || 1) {

          // for now just take the match with best score
          start_positions.push_back(atoi(tmp[2].c_str()));
          all_refs.push_back(firstCol[0]);
          break;
      }  
    }
  }

  system("rm -rf tmp");
  system("rm tmp*");

  cout << (double)(num_all_reads - still_unmatched) / (double)num_all_reads << "% after running against possible exon regions" << endl;
  fprintf(summaryFile, "%f %% reads matched after running on discovered exon regions\n", (num_all_reads - still_unmatched) / (double)num_all_reads );

  clock_t END_TIME = clock();
  cout << "Running time: " << (double)(END_TIME - START_TIME) / CLOCKS_PER_SEC << endl;
  
  return 0;  
}
