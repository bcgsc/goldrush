#include "opt.hpp"
#include <cmath>
#include <iostream>

namespace opt {

size_t assigned_max = 5;
size_t unassigned_min = 5;
size_t tile_length = 1000;
uint64_t hash_universe = 0;
uint64_t target_size = 0;
size_t kmer_size = 0;
size_t weight = 0;
size_t min_length = 5000;
size_t hash_num = 1;
double occupancy = 0.1;
double ratio = 0.1;
size_t jobs = 1;
size_t block_size = 10;
size_t max_paths = 1;
size_t threshold = 10;
uint32_t phred_min = 10;
uint32_t phred_delta = 5;
std::string prefix_file = "";
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
    << "  -k K -w W -i INPUT [-p prefix] [-o O] [-t T] [-h H] [-u U] [-m M]  "
       "[-a A] [-j J]\n\n"
       "  -i INPUT    find golden paths from INPUT [required]\n"
       "  -o O        use O as occupancy[0.1]\n"
       "  -h H        use h as number of spaced seed patterns [1]\n"
       "  -t T        use T as tile length [1000]\n"
       "  -k K        use K as span of spaced seed [required]\n"
       "  -w W        use W as weight of spaced seed [required]\n"
       "  -m M        use reads longer than M [5000]\n"
       "  -u U        U minimum unassigned tiles for read to be unassigned "
       "[5]\n"
       "  -a A        A maximum assigned tiles for read to be unassigned [5]\n"
       "  -p prefix   write output to files with prefix, e.g.\n"
       "  -P phred    minimum averge phred score for each read\n"
       "prefix_golden_path_0.fa [workpackage2]\n"
       "  -j J        use J number of threads [1]\n"
       "  -x X        require X hits for a tile to be assigned [10]\n"
       "  --ntcard    use ntcard to estimate genome size [false, assume max "
       "entries]\n"
       "  --silver_path    generate silver path instead of golden path. Silver paths terminate when the number bases recruited equals or exceeds T * r"
       "  --help      display this help and exit\n";
}

void
process_options(int argc, char** argv)
{
  int optindex = 0;
  int c;
  char* end = nullptr;
  while ((c = getopt_long(argc,
                          argv,
                          "a:b:d:f:g:h:i:j:k:m:M:o:r:s:t:u:w:x:p:P:T:",
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
      case 'g':
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
      case 'T':
        opt::target_size = (uint64_t)strtod(optarg, &end);
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
    exit(0);
  }

  if (!opt::weight) {
    std::cerr << "weight of spaced seed cannot be 0" << std::endl;
    print_usage("goldrush_path");
    exit(0);
  }
}
