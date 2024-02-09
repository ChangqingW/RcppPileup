#include "htslib/sam.h"
#include <Rcpp.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
// [[Rcpp::plugins(cpp17)]]

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

//' list alleles and the read id that supports them at a given pos
//'
//' @export
// [[Rcpp::export]]
int pileup_readid(Rcpp::String bam_path, Rcpp::String seqname, int pos) {

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

  int actual_tid = -1, n = -1, j = 0, refpos = -1;
  while ((plp = bam_plp_auto(plpiter, &actual_tid, &refpos, &n))) {
    if (actual_tid != tid || refpos != pos) {
      continue;
    }

    // iterate all reads and print the read id and the allele
    for (j = 0; j < n; j++) {
      if (plp[j].is_del) {
        Rcpp::Rcout << "-"
                    << "\t" << bam_get_qname(plp[j].b) << "\n";
        continue;
      }
      Rcpp::Rcout << seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]
                  << "\t" << bam_get_qname(plp[j].b) << "\n";
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

  return 0;
}