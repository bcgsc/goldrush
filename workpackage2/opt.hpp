#ifndef OPT_HPP
#define OPT_HPP

#include <cstdlib>
#include <getopt.h>
#include <string>

namespace opt {

extern size_t assigned_max;
extern size_t unassigned_min;
extern size_t tile_length;
extern size_t genome_size;
extern size_t target_size;
extern size_t kmer_size;
extern size_t weight;
extern size_t min_length;
extern size_t hash_num;
extern double occupancy;
extern double ratio;
extern size_t levels;
extern size_t jobs;
extern size_t max_paths;
extern size_t threshold;
extern std::string prefix_file;
extern std::string input;
extern std::string seed_preset;
extern int help;
extern int ntcard;
extern int second_pass;
extern int temp_mode;
extern int new_temp_mode;

}

static const struct option longopts[] = {
  { "new_temp_mode", no_argument, &opt::new_temp_mode, 1 },
  { "temp_mode", no_argument, &opt::temp_mode, 1 },
  { "second_pass", no_argument, &opt::second_pass, 1 },
  { "help", no_argument, &opt::help, 1 },
  { "ntcard", no_argument, &opt::ntcard, 1 },
  { nullptr, 0, nullptr, 0 }
};

static constexpr bool verbose = false;

void
process_options(int argc, char** argv);

#endif