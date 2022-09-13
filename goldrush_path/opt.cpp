#include "opt.hpp"
#include <cmath>
#include <iostream>

namespace opt {

size_t assigned_max = 1;
size_t unassigned_min = 5;
size_t tile_length = 1000;
uint64_t hash_universe = 0;
uint64_t genome_size = 0;
size_t kmer_size = 0;
size_t weight = 0;
size_t min_length = 20000;
size_t hash_num = 3;
double occupancy = 0.1;
double ratio = 0.9;
size_t jobs = 48;
size_t block_size = 10;
size_t max_paths = 1;
size_t threshold = 10;
uint32_t phred_min = 15;
uint32_t phred_delta = 5;
std::string prefix_file = "goldrush_out";
std::string input = "";
std::string seed_preset = "";
int help = 0;
int ntcard = 0;
int silver_path = 0;
std::string filter_file = "";

}

void
print_usage(const std::string& progname)
{
  std::cout
    << "Usage:  " << progname
    << "  -k K -w W -i INPUT -g G [-p prefix] [-P PHRED_AVG] [-o O] [-t T] [-f "
       "F] [-h H] [-u U] [-m M] [-H HASH_UNIVERSE] [-s S] [-x X] [-M MAX_PATHS]"
       "[-a A] [-j J] [-b B] [-d D] [--silver_path] [--ntcard] [--help] \n\n"
       "  -i INPUT                find golden paths from INPUT [required]\n"
       "  -g G                    estimated genome size [required]\n"
       "  -b B                    during insertion, B number of consecutive tiles to be inserted with the same ID [10]\n"
       "  -d D                    remove reads with greater or equal then D phred average "
       "between first half and second half of the read [5]\n"
       "  -f F                    don't use reads from F. Expects one read per line\n"
       "  -o O                    use O as occupancy [0.1]\n"
       "  -h H                    use h as number of spaced seed patterns [1]\n"
       "  -H HASH_UNIVERSE        determine MiBF size based on HASH_UNIVERSE "
       "[Calculated based on W and h]\n"
       "  -t T                    tile length [1000]\n"
       "  -k K                    span of spaced seed [required]\n"
       "  -w W                    weight of spaced seed [required]\n"
       "  -m M                    use reads longer than M [20000]\n"
       "  -u U                    U minimum unassigned tiles for read to be unassigned "
       "[5]\n"
       "  -a A                    A maximum assigned tiles for read to be unassigned [1]\n"
       "  -p prefix               write output to files with prefix "
       "[goldrush_out]\n"
       "  -P PHRED_AVG            minimum average phred score for each read [15]\n"
       "  -j J                    number of threads [48]\n"
       "  -s S                    use S seed preset. Must be consistent with k and w [n/a, "
       "generate one randomly based on k and w]\n"
       "  -x X                    require X hits for a tile to be assigned [10]\n"
       "  -M MAX_PATHS            output MAX_PATHS [5, used with --silver_path]\n"
       "  --ntcard                use ntcard to estimate genome size [false, assume max "
       "entries]\n"
       "  --silver_path           generate silver path(s) instead of golden path. "
       "Silver paths terminate when the number of bases recruited equals or "
       "exceeds T * r\n"
       "  --help                  display this help and exit\n";
}

void
process_options(int argc, char** argv)
{
  int optindex = 0;
  int c;
  char* end = nullptr;
  while ((c = getopt_long(argc,
                          argv,
                          "a:b:d:f:g:h:i:j:k:m:M:o:r:s:t:u:w:x:p:P:H:",
                          longopts,
                          &optindex)) != -1) {
    switch (c) {
      case 0:
        break;
      case 'a': {
        opt::assigned_max = strtoul(optarg, &end, 10);
        break;
      }
      case 'b':
        opt::block_size = strtoul(optarg, &end, 10);
        break;
      case 'd':
        opt::phred_delta = strtoul(optarg, &end, 10);
        break;
      case 'f':
        opt::filter_file = optarg;
        break;
      case 'H':
        opt::hash_universe = strtoull(optarg, &end, 10);
        break;
      case 'h':
        opt::hash_num = strtoul(optarg, &end, 10);
        break;
      case 'i':
        opt::input = optarg;
        break;
      case 'j':
        opt::jobs = strtoul(optarg, &end, 10);
        break;
      case 'k':
        opt::kmer_size = strtoul(optarg, &end, 10);
        break;
      case 'm':
        opt::min_length = strtoul(optarg, &end, 10);
        break;
      case 'M':
        opt::max_paths = strtoul(optarg, &end, 10);
        break;
      case 'o':
        opt::occupancy = strtod(optarg, &end);
        break;
      case 'r':
        opt::ratio = strtod(optarg, &end);
        break;
      case 'p':
        opt::prefix_file = optarg;
        break;
      case 'P':
        opt::phred_min = strtoul(optarg, &end, 10);
        break;
      case 's':
        opt::seed_preset = optarg;
        break;
      case 't':
        opt::tile_length = strtoul(optarg, &end, 10);
        break;
      case 'g':
        opt::genome_size = (uint64_t)strtod(optarg, &end);
        break;
      case 'u': {
        opt::unassigned_min = strtoul(optarg, &end, 10);
        break;
      }
      case 'w': {
        opt::weight = strtoul(optarg, &end, 10);
        break;
      }
      case 'x': {
        opt::threshold = strtoul(optarg, &end, 10);
        break;
      }
      default:
        exit(EXIT_FAILURE);
    }
  }

  if (opt::help) {
    print_usage("goldrush_path");
    exit(0);
  }

  if (!opt::kmer_size) {
    std::cerr << "span of spaced seed cannot be 0" << std::endl;
    print_usage("goldrush_path");
    exit(1);
  }

  if (!opt::weight) {
    std::cerr << "weight of spaced seed cannot be 0" << std::endl;
    print_usage("goldrush_path");
    exit(1);
  }

  if (opt::genome_size == 0) {
    std::cerr << "genome size cannot be 0" << std::endl;
    print_usage("goldrush_path");
    exit(1);
  }

  if (!opt::seed_preset.empty()) {
    if (opt::kmer_size != opt::seed_preset.size()) {
      std::cerr << "seed preset must be the same size of k" << std::endl;
      print_usage("goldrush_path");
      exit(1);
    }
    uint8_t num_1s_in_seed = 0;
    for (char c : opt::seed_preset) {
      if (c == '1') {
        ++num_1s_in_seed;
      }
    }
    if (opt::weight != num_1s_in_seed) {
      std::cerr << "seed preset must have the same weight as w" << std::endl;
      print_usage("goldrush_path");
      exit(1);
    }
  }
}
