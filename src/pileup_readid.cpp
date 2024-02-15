#include "htslib/sam.h"
#include <Rcpp.h>
#include <array>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>
// [[Rcpp::plugins(cpp17)]]

const std::vector<std::string> BASES = {"A", "T", "C", "G", "del"};

typedef struct plpconf {
  const char *inname;
  samFile *infile;
  sam_hdr_t *in_samhdr;
  hts_idx_t *in_idx;
  hts_itr_t *iter;
} plpconf;

int plpconstructor(void *data, const bam1_t *b, bam_pileup_cd *cd) { return 0; }
int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd) { return 0; }

int readdata(void *data, bam1_t *b) {
  plpconf *conf = (plpconf *)data;
  if (!conf || !conf->infile) {
    return -2; // cant read data
  }

  bool skip = true;
  int ret = 0;
  do {
    ret = conf->iter ? sam_itr_next(conf->infile, conf->iter, b)
                     : sam_read1(conf->infile, conf->in_samhdr, b);
    if (b->core.tid < 0 ||
        (b->core.flag & BAM_FUNMAP)) { // exclude unmapped reads
      skip = true;
      continue;
    } else {
      skip = false;
    }
  } while (skip);

  return ret;
}

// list alleles and the read id that supports them at a given pos
std::unordered_map<std::string, std::vector<std::string>>
pileup_readid_hashmap(Rcpp::String bam_path, Rcpp::String seqname, int pos,
                      bool indel) {
  bam1_t *bamdata = NULL;
  plpconf conf = {0};
  bam_plp_t plpiter = NULL;
  const bam_pileup1_t *plp = NULL;

  conf.inname = bam_path.get_cstring();
  const char *seqname_c = seqname.get_cstring();

  if (!(bamdata = bam_init1())) {
    Rcpp::stop("Failed to initialize bamdata");
  }
  if (!(conf.infile = sam_open(conf.inname, "r"))) {
    Rcpp::stop("Failed to open file %s\n", conf.inname);
  }
  if (!(conf.in_samhdr = sam_hdr_read(conf.infile))) {
    Rcpp::stop("Failed to read header from file!\n");
  }
  if (!(plpiter = bam_plp_init(readdata, &conf))) {
    Rcpp::stop("Failed to initialize pileup data\n");
  }
  bam_plp_set_maxcnt(plpiter, INT_MAX); // caps at 2b
  if (!(conf.in_idx = sam_index_load(conf.infile, conf.inname))) {
    Rcpp::stop("Failed to load index for %s", conf.inname);
  }
  int tid = bam_name2id(conf.in_samhdr, seqname_c);
  if ((conf.iter = sam_itr_queryi(conf.in_idx, tid, pos, pos + 1)) == 0) {
    Rcpp::stop("Failed to parse region");
  }

  // set constructor destructor callbacks
  bam_plp_constructor(plpiter, plpconstructor);
  bam_plp_destructor(plpiter, plpdestructor);

  // store the alleles and the read id that supports them
  std::unordered_map<std::string, std::vector<std::string>> hashmap;
  // char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
  const std::string seq_nt16_string[] = {"=", "A", "C", "M", "G", "R",
                                         "S", "V", "T", "W", "Y", "H",
                                         "K", "D", "B", "N"};

  int actual_tid = -1, n = -1, j = 0, refpos = -1;
  while ((plp = bam_plp_auto(plpiter, &actual_tid, &refpos, &n))) {
    if (actual_tid != tid || refpos != pos) {
      continue;
    }

    // iterate all reads and print the read id and the allele
    for (j = 0; j < n; j++) {

      std::string read_id = bam_get_qname(plp[j].b);
      std::size_t id_idx = read_id.find("#");
      if (id_idx == std::string::npos) {
        Rcpp::stop("Unexpected read id format: %s", read_id);
      }

      if (plp[j].is_refskip) {
        continue;
      }

      if (!indel) { // count SNV
        if (plp[j].is_del) {
          hashmap["del"].push_back(read_id.substr(0, id_idx));
        } else {
          hashmap[seq_nt16_string[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]]
              .push_back(read_id.substr(0, id_idx));
        }
      } else { // count indel
        if (plp[j].indel == 0) {
          hashmap["."].push_back(read_id.substr(0, id_idx));
        } else if (plp[j].indel > 0) {
          // key as '+[inserted_bases]' e.g. '+A' or '+AT'
          std::string inserted_bases = "";
          for (int k = 1; k <= plp[j].indel; k++) {
            inserted_bases += seq_nt16_string[bam_seqi(bam_get_seq(plp[j].b),
                                                       plp[j].qpos + k)];
          }
          hashmap["+" + inserted_bases].push_back(read_id.substr(0, id_idx));
        } else if (plp[j].indel < 0) {
          // key as '-[deleted_number_of_bases]' e.g. '-2'
          hashmap["-" + std::to_string(-plp[j].indel)].push_back(
              read_id.substr(0, id_idx));
          // if position x has indel of -2,
          // the next 2 positions will be is_del
        }
      }
    }
  }

  if (conf.in_samhdr) {
    sam_hdr_destroy(conf.in_samhdr);
  }
  if (conf.infile) {
    sam_close(conf.infile);
  }
  if (conf.in_idx) {
    hts_idx_destroy(conf.in_idx);
  }
  if (conf.iter) {
    hts_itr_destroy(conf.iter);
  }
  if (bamdata) {
    bam_destroy1(bamdata);
  }
  if (plpiter) {
    bam_plp_destroy(plpiter);
  }

  return hashmap;
}

//' list alleles and the read id that supports them at a given pos
//' @export
// [[Rcpp::export]]
Rcpp::List pileup_readid(Rcpp::String bam_path, Rcpp::String seqname, int pos,
                         bool indel) {
  std::unordered_map<std::string, std::vector<std::string>> hashmap =
      pileup_readid_hashmap(bam_path, seqname, pos, indel);
  Rcpp::List ret = Rcpp::wrap(hashmap);
  return ret;
}

// works for both snv and indel
bool barcode_has_variant(const std::string &barcode,
                         const std::vector<std::string> &read_ids,
                         const std::size_t &umi_idx) {
  if (read_ids.empty()) {
    return false;
  }
  for (const auto &read_id : read_ids) {
    if (read_id.substr(0, umi_idx) == barcode) {
      return true;
    }
  }
  return false;
}

// from flexiplex.cpp; for UMI deduplication
unsigned int edit_distance(const std::string &s1, const std::string &s2,
                           unsigned int max_editd) {
  const std::string_view s1_view(s1);
  const std::string_view s2_view(s2);
  const std::size_t len1 = s1_view.size() + 1;
  const std::size_t len2 = s2_view.size() + 1;
  std::vector<unsigned int> dist_holder(len1 * len2);
  dist_holder[0] = 0; //[0][0]
  for (std::size_t j = 1; j < len2; ++j)
    dist_holder[j] = j; //[0][j];
  for (std::size_t i = 1; i < len1; ++i)
    dist_holder[i * len2] = 0; //[i][0];
  unsigned int best = len2;

  // loop over the distance matrix elements and calculate running distance
  for (std::size_t j = 1; j < len2; ++j) {
    bool any_below_threshold = false; // flag used for early exit
    for (std::size_t i = 1; i < len1; ++i) {
      unsigned int sub =
          (s1_view[i - 1] == s2_view[j - 1]) ? 0 : 1; // match / mismatch score

      const unsigned int &top_left = dist_holder[(i - 1) * len2 + (j - 1)];
      const unsigned int &left = dist_holder[i * len2 + (j - 1)];
      const unsigned int &top = dist_holder[(i - 1) * len2 + j];

      unsigned int min_value = std::min({top + 1, left + 1, top_left + sub});
      dist_holder[i * len2 + j] = min_value;

      if (min_value <= max_editd)
        any_below_threshold = true;
      if (j == (len2 - 1) && min_value < best) {
        best = min_value;
      }
    }
    if (!any_below_threshold) {
      return max_editd + 1;
    }
  }
  return best;
}

// compute deduped counts for each base
std::unordered_map<std::string, unsigned int> variant_umi_count(
    std::unordered_map<std::string, std::vector<std::string>> variant_readid,
    std::string barcode, const std::size_t &umi_idx) {

  const int MAX_EDITD = 2;
  std::unordered_map<std::string, unsigned int> deduped_counts;
  std::unordered_map<std::string, std::vector<std::string>> umis;

  for (const auto &variant : variant_readid) {
    for (const auto &current_read_id : variant.second) {
      if (current_read_id.substr(0, umi_idx) != barcode) {
        // not the same cell barcode
        continue;
      }
      std::string current_umi = current_read_id.substr(umi_idx + 1);
      if (umis.find(variant.first) == umis.end()) {
        // first umi
        umis[variant.first] = {current_umi};
        deduped_counts[variant.first]++;
      } else {
        // dedup subsequent umis
        std::vector<unsigned int> current_editds;
        std::transform(umis[variant.first].begin(), umis[variant.first].end(),
                       std::back_inserter(current_editds),
                       [&current_umi](const auto &umi) {
                         return edit_distance(umi, current_umi, MAX_EDITD);
                       });
        unsigned int min_editd =
            *std::min_element(current_editds.begin(), current_editds.end());
        if (min_editd <= MAX_EDITD) {
          // deduped
          continue;
        } else {
          umis[variant.first].push_back(current_umi);
          deduped_counts[variant.first]++;
        }
      }
    }
  }

  return deduped_counts;
}

// compute proportion of umi for each variant
inline std::unordered_map<std::string, double> variant_umi_proportion(
    std::unordered_map<std::string, std::vector<std::string>> variant_readid,
    std::string barcode, const std::size_t &umi_idx) {
  std::unordered_map<std::string, unsigned int> deduped_counts =
      variant_umi_count(variant_readid, barcode, umi_idx);
  std::unordered_map<std::string, double> deduped_proportion;

  // devid counts by total
  double total = 0;
  for (const auto &variant : deduped_counts) {
    total += deduped_counts[variant.first];
  }
  for (const auto &variant : deduped_counts) {
    deduped_proportion[variant.first] = deduped_counts[variant.first] / total;
  }
  return deduped_proportion;
}

// compute number of reads for each variant
std::unordered_map<std::string, unsigned int> variant_read_count(
    std::unordered_map<std::string, std::vector<std::string>> variant_readid,
    std::string barcode, const std::size_t &umi_idx) {

  std::unordered_map<std::string, unsigned int> counts;
  for (const auto &variant : variant_readid) {
    for (const auto &read_id : variant.second) {
      if (read_id.substr(0, umi_idx) == barcode) {
        counts[variant.first]++;
      }
    }
  }
  return counts;
}

//' produce a matrix of variants and cell barcodes
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix variant_proportion_matrix(Rcpp::String bam_path, Rcpp::String seqname,
                                   int pos, Rcpp::StringVector barcodes,
                                   bool indel) {
  std::unordered_map<std::string, std::vector<std::string>> hashmap =
      pileup_readid_hashmap(bam_path, seqname, pos, indel);
  if (hashmap.empty()) {
    Rcpp::stop("No reads found at %s:%i", seqname.get_cstring(),
               pos); // TODO: return empty matrix instead
  }
  const std::size_t umi_idx = hashmap.begin()->second[0].find("_");

  std::vector<std::string> variant_names;
  if (!indel) {
    variant_names = BASES;
  } else {
    variant_names = {"."};
    for (const auto &variant : hashmap) {
      if (std::find(variant_names.begin(), variant_names.end(),
                    variant.first) == variant_names.end()) {
        variant_names.push_back(variant.first);
      }
    }
  }

  Rcpp::NumericMatrix ret(variant_names.size(), barcodes.size());
  Rcpp::colnames(ret) = barcodes;
  Rcpp::rownames(ret) = Rcpp::wrap(variant_names);

  // 0 if no read from that cell barcode support the variant
  // 1 otherwise
  for (unsigned int i = 0; i < barcodes.size(); i++) {
    Rcpp::String current_barcode = barcodes[i];

    int current_variants_sum = 0;
    for (unsigned int j = 0; j < variant_names.size(); j++) {
      if (hashmap.find(variant_names[j]) == hashmap.end()) {
        // no such variant in any cell barcode
        continue;
      } else if (barcode_has_variant(current_barcode, hashmap[variant_names[j]],
                                     umi_idx)) {
        ret(j, i) = 1;
        current_variants_sum++;
      }
    }

    if (current_variants_sum > 1) {
      std::unordered_map<std::string, double> deduped_counts =
          variant_umi_proportion(hashmap, current_barcode, umi_idx);
      for (auto j = 0; j < variant_names.size(); j++) {
        if (deduped_counts.find(variant_names[j]) != deduped_counts.end()) {
          ret(j, i) = deduped_counts[variant_names[j]];
        } else {
          ret(j, i) = 0;
        }
      }
    }
  }

  return ret;
}
