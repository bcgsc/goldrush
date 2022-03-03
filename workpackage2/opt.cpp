#include "opt.hpp"

#include <iostream>

namespace opt {

size_t assigned_max = 5;
size_t unassigned_min = 5;
size_t tile_length = 1000;
size_t genome_size = 0;
size_t target_size = 0;
size_t kmer_size = 0;
size_t weight = 0;
size_t min_length = 5000;
size_t hash_num = 1;
double occupancy = 0.1;
double ratio = 0.1;
size_t levels = 1;
size_t jobs = 1;
size_t block_size = 10;
size_t max_paths = 1;
size_t threshold = 10;
std::string prefix_file = "workpackage2";
std::string input = "";
std::string seed_preset = "";
int help = 0;
int ntcard = 0;
int second_pass = 0;
int temp_mode = 0;
int new_temp_mode = 0;
std::string filter_file ="";

}

void
print_usage(const std::string& progname)
{
  std::cout
    << "Usage:  " << progname
    << "  -k K -w W -i INPUT [-p prefix] [-o O] [-t T] [-h H] [-u U] [-m M]  "
       "[-a A] [-l L] [-j J]\n\n"
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
       "  -l L        output L golden paths [1]\n"
       "  -p prefix   write output to files with prefix, e.g. "
       "prefix_golden_path_0.fa [workpackage2]\n"
       "  -j J        use J number of threads [1]\n"
       "  -x X        require X hits for a tile to be assigned [10]\n"
       "  --ntcard    use ntcard to estimate genome size [false, assume max "
       "entries]\n"
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
                          "a:b:f:g:h:i:j:k:l:m:M:o:r:s:t:u:w:x:p:T:",
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
      case 'f':
          opt::filter_file = optarg;
          break;
      case 'g':
        opt::genome_size = strtoull(optarg, &end, 10);
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
      case 'l':
        opt::levels = strtoul(optarg, &end, 10);
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
      case 's':
        opt::seed_preset = optarg;
        break;
      case 't':
        opt::tile_length = strtoul(optarg, &end, 10);
        break;
      case 'T':
        opt::target_size = strtoul(optarg, &end, 10);
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
    print_usage("workpackage2");
    exit(0);
  }

  if (!opt::kmer_size) {
    std::cerr << "span of spaced seed cannot be 0" << std::endl;
    print_usage("workpackage2");
    exit(0);
  }

  if (!opt::weight) {
    std::cerr << "weight of spaced seed cannot be 0" << std::endl;
    print_usage("workpackage2");
    exit(0);
  }
}