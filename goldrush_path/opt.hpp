#ifndef OPT_HPP
#define OPT_HPP

#include <cstdlib>
#include <getopt.h>
#include <string>

namespace opt {

extern size_t assigned_max;
extern size_t unassigned_min;
extern size_t tile_length;
extern size_t block_size;
extern uint64_t hash_universe;
extern uint64_t genome_size;
extern size_t kmer_size;
extern uint32_t phred_min;
extern uint32_t phred_delta;
extern size_t weight;
extern size_t min_length;
extern size_t hash_num;
extern double occupancy;
extern double ratio;
extern size_t jobs;
extern size_t max_paths;
extern size_t threshold;
extern std::string prefix_file;
extern std::string input;
extern std::string seed_preset;
extern std::string filter_file;
extern int help;
extern int ntcard;
extern int silver_path;

}

static const struct option longopts[] = {
  { "silver_path", no_argument, &opt::silver_path, 1 },
  { "help", no_argument, &opt::help, 1 },
  { "ntcard", no_argument, &opt::ntcard, 1 },
  { nullptr, 0, nullptr, 0 }
};

static constexpr bool verbose = true;

void
process_options(int argc, char** argv);

#endif
