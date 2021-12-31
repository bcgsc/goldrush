#include "btllib/seq_reader.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/util.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream> 
#include <getopt.h>

#include <unordered_set>

#if _OPENMP
#include <omp.h>
#endif

#include "multiLensfrHashIterator.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>
#include <sstream>
#include <getopt.h>
#include <sys/stat.h>
#include "MIBFClassifier.hpp"
#include "config.h"
#include "Common/Options.h"
#include "bloom_filter.hpp"
#include "MIBloomFilter.hpp"
#include "MIBFConstructSupport.hpp"
//#include "btl_bloomfilter/vendor/stHashIterator.hpp"
//#include "Common/sntHashIterator.hpp"

//#include "btl_bloomfilter/BloomFilter.hpp"

#include "Common/Options.h"

#include "tsl/robin_map.h"
#include "tsl/robin_set.h"
#include "ntcard.hpp"

#include <tuple>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>
#include <google/dense_hash_set>
#include <sdsl/int_vector.hpp>
//#include "Common/sntHashIterator.hpp"
#include "MIBFConstructSupport.hpp"


#include <zlib.h>
#include <stdio.h>

namespace opt {
    size_t assigned_max = 5;
    size_t unassigned_min = 5;
    size_t tile_length = 1000;
    size_t genome_size = 0;
    size_t kmer_size = 0;
    size_t weight = 0;
    size_t min_length = 5000;
    size_t hash_num = 1;
    double occupancy = 0.1;
    size_t levels = 1;
    size_t jobs = 1;
    size_t threshold = 10;
    std::string prefix_file = "workpackage2";
    std::string input = "";
    std::string seed_preset = "";
    int help = 0;
    int ntcard = 0;
    int second_pass = 0;

}

bool sort_by_sec(const pair<size_t, size_t> &a,
              const pair<size_t, size_t> &b)
{
    return (a.second > b.second);
}

static void
printUsage(const std::string& progname)
{
	std::cout << "Usage:  " << progname
	          << "  -k K -w W -i INPUT [-p prefix] [-o O] [-t T] [-h H] [-u U] [-m M]  [-a A] [-l L] [-j J]\n\n"
                 "  -i INPUT    find golden paths from INPUT [required]\n"
                 "  -o O        use O as occupancy[0.1]\n"
                 "  -h H        use h as number of spaced seed patterns [1]\n"
                 "  -t T        use T as tile length [1000]\n"
	             "  -k K        use K as span of spaced seed [required]\n"
	             "  -w W        use W as weight of spaced seed [required]\n"
                 "  -m M        use reads longer than M [5000]\n"
                 "  -u U        U minimum unassigned tiles for read to be unassigned [5]\n"
	             "  -a A        A maximum assigned tiles for read to be unassigned [5]\n"
	             "  -l L        output L golden paths [1]\n"
	             "  -p prefix   write output to files with prefix, e.g. prefix_golden_path_0.fa [workpackage2]\n"
	             "  -j J        use J number of threads [1]\n"
                 "  -x X        require X hits for a tile to be assigned [10]\n"
                 "  --ntcard    use ntcard to estimate genome size [false, assume max entries]\n"
	             "  --help      display this help and exit\n";
}


uint64_t calc_ntcard_genome_size(const std::vector<std::string>& seed_string_vec) {

    vector<string> inFiles;
    inFiles.push_back(opt::input);
    uint64_t genome_size = 0;

    /* indices 0 and 1 are reserved for F0 and 1.
       Indicex 2 to 10001 store kmer multiplicities up to 10000 */
    static const size_t DEFAULT_NTCARD_HIST_COV_MAX = 10000;
    std::cerr << "Calculating expected entries" << std::endl;
    const auto histArray = getHist(inFiles, opt::kmer_size, opt::jobs, DEFAULT_NTCARD_HIST_COV_MAX, seed_string_vec);
    for (size_t i = 0; i < seed_string_vec.size(); ++i) {
        std::cerr << "Expected entries for seed pattern " << seed_string_vec[i] << " : " << histArray[i][1] << std::endl;
       
        genome_size += (histArray[i][1]);
    }
    std::cerr << "Total expected entries for seed patterns: " << genome_size << std::endl;
    return genome_size;


}

const std::vector<std::string> make_seed_pattern () {

    std::vector<std::string> seed_string_vec;
    std::string left_seed_str;
    std::string right_seed_str;

    if (opt::seed_preset == "") {
        srand(time(NULL));
        // seed generation
        std::cerr << "Designing base symmetrical spaced seed" << "\n"
                << "Using:" << "\n"
                << "span: " << opt::kmer_size << "\n"
                << "weight: " << opt::weight << std::endl;

        std::vector<unsigned> left_seed_vec (opt::kmer_size / 2, 0);
        left_seed_vec[0] = 1; // left most val in seed must be a 1
        size_t weight_count = 0;

        while (weight_count != opt::weight / 2 ){
            for (size_t i = 1; i < opt::kmer_size / 2; ++i) {
                left_seed_vec[i] = rand() % 2;            
            }
            weight_count = std::count (left_seed_vec.begin(), left_seed_vec.end(), 1);
        }
        
        std::stringstream temp_ss;
        for(const auto val : left_seed_vec) {
            temp_ss << val;
        } 

        left_seed_str = temp_ss.str();
        right_seed_str = std::string(left_seed_str.rbegin(), left_seed_str.rend());      

    } else {
        std::cerr << "Using preset spaced seed" << "\n"
                << "with:" << "\n"
                << "span: " << opt::seed_preset.size() << "\n"
                << "weight: " << std::count (opt::seed_preset.begin(), opt::seed_preset.end(), '1') << std::endl;
       left_seed_str = opt::seed_preset.substr(0, opt::seed_preset.size() / 2);
       right_seed_str = opt::seed_preset.substr(opt::seed_preset.size() / 2, opt::seed_preset.size() / 2);        
    }

    for (size_t i = 0; i < opt::hash_num; ++i) {
        seed_string_vec.push_back(left_seed_str + std::string(i , '0') + right_seed_str);
    }

    return  seed_string_vec;
}


size_t calc_num_assigned_tiles (const std::unique_ptr<MIBloomFilter<uint32_t>>& miBF, const std::vector<std::vector<uint64_t>> hashed_values, std::vector<uint32_t>& tiles_assigned_id_vec, std::vector<bool>& tiles_assigned_bool_vec) {

    
    size_t num_assigned_tiles = 0;
    size_t num_tiles = hashed_values.size();
    //std::vector<uint32_t> tiles_assigned_id_vec (num_tiles, 0);


#if _OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < num_tiles; ++i) { // for each tile except first and last tiles. We consider them erroneous and not used in tile assignment.

        tsl::robin_map<uint32_t, std::pair<uint32_t, uint32_t>> id_counts; // store counts of ids

        // Reusable vector for ranks 
        vector<uint64_t> m_rank_pos(miBF->getHashNum());
        // Reusable vector for IDs 
        vector<uint32_t> m_data(miBF->getHashNum());
        
        const auto& hashed_values_flat_array = hashed_values[i];
        std::vector<uint64_t> hashes(opt::hash_num, 0);
        for (size_t curr_frame = 0; curr_frame < (hashed_values_flat_array.size() / opt::hash_num); ++curr_frame) {
            
            for (size_t curr_hash = 0; curr_hash < opt::hash_num; ++curr_hash) {
                hashes[curr_hash] = hashed_values_flat_array[curr_frame * opt::hash_num + curr_hash];
            }

            tsl::robin_set<uint32_t> unique_ids;
            if (miBF->atRank(hashes, m_rank_pos)) {                      // if its a hit
                m_data = miBF->getData(m_rank_pos);                 // m_data has ID's
                for(unsigned m = 0; m < miBF->getHashNum(); m++){   // iterate over ID's
                    if(m_data[m] > miBF->s_mask){                     // if ID is saturated
                        uint32_t new_id = m_data[m] & miBF->s_antiMask;
                        if (new_id == 0) {
                            continue;
                        }                                
                        unique_ids.insert(new_id);
                    }
                    else{
                        if (m_data[m] == 0) {
                            continue;
                        }
                        unique_ids.insert(m_data[m]);
                    }
                }
            }

            for (const auto& unique_id : unique_ids) { // tabulate all unique ids to count table
                if (id_counts.find(unique_id) != id_counts.end()) {
                    ++id_counts[unique_id].second;
                } else {
                    id_counts[unique_id] = std::make_pair(unique_id, 1);
                }
            }

        } 

        uint32_t curr_id = 0;
        uint32_t curr_id_count = 0;

        for (auto id_counts_it = id_counts.begin(); id_counts_it != id_counts.end(); ++id_counts_it) { //find id with highest count in a tile
            if (id_counts_it->second.second > curr_id_count) {
                curr_id = id_counts_it->first;
                curr_id_count = id_counts_it->second.second;
            }
        }

        if (curr_id_count > opt::threshold) {
#if _OPENMP
#pragma omp atomic
#endif
            num_assigned_tiles += 1;
            tiles_assigned_bool_vec[i] = true;
        }

        tiles_assigned_id_vec[i] = curr_id;
        
    }

    //print which id each tile is assigned to
    for (const auto& tiles_assigned_id : tiles_assigned_id_vec) {
        std::cerr << tiles_assigned_id << "\t";


    }
    std::cerr << std::endl;

    for (const auto& tiles_assigned_bool : tiles_assigned_bool_vec) {
        std::cerr << tiles_assigned_bool << "\t";


    }
    std::cerr << std::endl;
    num_assigned_tiles = 0;
    for (const auto& is_tile_assigned : tiles_assigned_bool_vec) {
        if (is_tile_assigned) {
            ++num_assigned_tiles;
        }
    }


    for (size_t i = 1; i < num_tiles -1; ++i) {
        if (tiles_assigned_bool_vec[i] == false && tiles_assigned_id_vec[i] > 100) {
            if ((tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i - 1] && tiles_assigned_bool_vec[i - 1] == true) || (tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i + 1] && tiles_assigned_bool_vec[i + 1] == true)) {
                tiles_assigned_bool_vec[i] = true;
            } else if ((tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i - 1] + 1 && tiles_assigned_bool_vec[i - 1] == true) || (tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i + 1] + 1 && tiles_assigned_bool_vec[i + 1] == true)) {
                tiles_assigned_bool_vec[i] = true;
            } else if ((tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i - 1] - 1 && tiles_assigned_bool_vec[i - 1] == true) || (tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i + 1] - 1 && tiles_assigned_bool_vec[i + 1] == true)) {
                tiles_assigned_bool_vec[i] = true;
            } else if (tiles_assigned_id_vec[i -1] == tiles_assigned_id_vec[i + 1]  && tiles_assigned_bool_vec[i -1] == true &&  tiles_assigned_bool_vec[i + 1] == true) {
                tiles_assigned_bool_vec[i] = true;
                tiles_assigned_id_vec[i] = tiles_assigned_id_vec[i -1];
            }
        }
    }

    for (size_t i = num_tiles - 2; i >= 1; --i) {
        if (tiles_assigned_bool_vec[i] == false  && tiles_assigned_id_vec[i] > 100) {
            if ((tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i - 1] && tiles_assigned_bool_vec[i - 1] == true) || (tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i + 1] && tiles_assigned_bool_vec[i + 1] == true)) {
                tiles_assigned_bool_vec[i] = true;
            } else if ((tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i - 1] + 1 && tiles_assigned_bool_vec[i - 1] == true) || (tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i + 1] + 1 && tiles_assigned_bool_vec[i + 1] == true)) {
                tiles_assigned_bool_vec[i] = true;
            } else if ((tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i - 1] - 1 && tiles_assigned_bool_vec[i - 1] == true) || (tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i + 1] - 1 && tiles_assigned_bool_vec[i + 1] == true)) {
                tiles_assigned_bool_vec[i] = true;
            } else if (tiles_assigned_id_vec[i -1] == tiles_assigned_id_vec[i + 1]  && tiles_assigned_bool_vec[i -1] == true &&  tiles_assigned_bool_vec[i + 1] == true) {
                tiles_assigned_bool_vec[i] = true;
                tiles_assigned_id_vec[i] = tiles_assigned_id_vec[i -1];
            }
        }
    }

    //print which id each tile is assigned to
    for (const auto& tiles_assigned_id : tiles_assigned_id_vec) {
        std::cerr << tiles_assigned_id << "\t";


    }
    std::cerr << std::endl;

    for (const auto& tiles_assigned_bool : tiles_assigned_bool_vec) {
        std::cerr << tiles_assigned_bool << "\t";


    }
    std::cerr << std::endl;

    size_t start_idx = 0;
    size_t end_idx = 0;
    //size_t curr_stretch = 0;
    std::vector<std::pair<size_t, size_t>> coord_vec;
    for (size_t i = 1; i < num_tiles - 1; ++i){
        if (tiles_assigned_bool_vec[i] == false && tiles_assigned_bool_vec[i - 1] == true) {
            start_idx = i;
        } else if (tiles_assigned_bool_vec[i] == true && tiles_assigned_bool_vec[i - 1] == false) {
            end_idx = i - 1;
            coord_vec.push_back(std::make_pair(start_idx, end_idx));
        }
    }
    for (const auto& coords : coord_vec) {
        const auto& left  = tiles_assigned_id_vec[coords.first - 1];
        const auto& right  = tiles_assigned_id_vec[coords.second + 1];
        if (left == right|| left == right + 1 || left == right - 1) {
            for (auto i = coords.first; i <= coords.second; ++i) {
                tiles_assigned_bool_vec[i] = true;
                tiles_assigned_id_vec[i] = left;
            }
        }

    }

    //print which id each tile is assigned to
    for (const auto& tiles_assigned_id : tiles_assigned_id_vec) {
        std::cerr << tiles_assigned_id << "\t";


    }
    std::cerr << std::endl;

    for (const auto& tiles_assigned_bool : tiles_assigned_bool_vec) {
        std::cerr << tiles_assigned_bool << "\t";


    }
    std::cerr << std::endl;


    for (size_t i = 2; i < num_tiles -2; ++i) {
        if (tiles_assigned_bool_vec[i] == true && tiles_assigned_id_vec[i] > 100) {
            if (tiles_assigned_bool_vec[i - 1] == false && tiles_assigned_bool_vec[i + 1] == false) {
                tiles_assigned_bool_vec[i] = false;
            }
        }
    }

    for (size_t i = num_tiles - 3; i >= 2; --i) {
        if (tiles_assigned_bool_vec[i] == true && tiles_assigned_id_vec[i] > 100) {
            if (tiles_assigned_bool_vec[i - 1] == false && tiles_assigned_bool_vec[i + 1] == false) {
                tiles_assigned_bool_vec[i] = false;
            }
        }
    }


    //print which id each tile is assigned to
    for (const auto& tiles_assigned_id : tiles_assigned_id_vec) {
        std::cerr << tiles_assigned_id << "\t";


    }
    std::cerr << std::endl;

    for (const auto& tiles_assigned_bool : tiles_assigned_bool_vec) {
        std::cerr << tiles_assigned_bool << "\t";


    }
    std::cerr << std::endl;

    std::cerr << std::endl;
    num_assigned_tiles = 0;
    for (const auto& is_tile_assigned : tiles_assigned_bool_vec) {
        if (is_tile_assigned) {
            ++num_assigned_tiles;
        }
    }
    return num_assigned_tiles;
}


int main(int argc, char** argv) {

	static const struct option longopts[] = { { "second_pass", no_argument, &opt::second_pass, 1 },
                                              { "help", no_argument, &opt::help, 1 },
                                              { "ntcard", no_argument, &opt::ntcard, 1 },
		                                      { nullptr, 0, nullptr, 0 } };

    int optindex = 0;
    int c;
    char* end = nullptr;
    while ((c = getopt_long(argc, argv, "a:g:h:i:j:k:l:m:o:s:t:u:w:x:p:", longopts, &optindex)) != -1) {
        switch (c) {
        case 0:
            break;
        case 'a': {
            opt::assigned_max = strtoul(optarg, &end, 10);
            break;
        }
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
            opt::min_length =  strtoul(optarg, &end, 10);
            break;
        case 'o':
            opt::occupancy = strtod(optarg, &end);
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
        printUsage("workpackage2");
        exit(0);
    }


    if (!opt::kmer_size) {
        std::cerr << "span of spaced seed cannot be 0" << std::endl;
        printUsage("workpackage2");
        exit(0);
    }

    if (!opt::weight) {
        std::cerr << "weight of spaced seed cannot be 0" << std::endl;
        printUsage("workpackage2");
        exit(0);
    }

    if (opt::second_pass) {
        std::cerr << "second_pass" << std::endl;
    }
#if _OPENMP
	omp_set_num_threads(opt::jobs);
#endif

    size_t min_seq_len = opt::min_length;


    //srand (1); // for testing, change to srand(time(NULL)) for actual code
    const std::vector<std::string> seed_string_vec = make_seed_pattern();


    
	if (opt::genome_size == 0) {
        if (opt::ntcard) {
            opt::genome_size = calc_ntcard_genome_size(seed_string_vec);
        } else {
            static const uint8_t BASES = 4;
            static const float HASH_UNIVERSE_COEFFICIENT = 0.5;
            opt::genome_size = pow(BASES, opt::weight) * HASH_UNIVERSE_COEFFICIENT * opt::hash_num;
        }
        
	}



    
    uint32_t id = 1;

    std::cerr << "Calculating " << opt::levels << " golden path(s)" << "\n"
              << "Using:" << "\n"
              << "tile length: " << opt::tile_length << "\n" 
              << "seed pattens: " << opt::hash_num << "\n"
              << "threshold: " << opt::threshold << "\n"
              << "base seed pattern: " << seed_string_vec[0] << "\n"
              << "minimum unassigned tiles: " << opt::unassigned_min << "\n"
              << "maximum assigned tiles: " << opt::assigned_max << "\n"
              << "genome size: " << opt::genome_size << "\n"
              << "occupancy: " << opt::occupancy << "\n"
              << "jobs: " << opt::jobs << std::endl;

            
    
    std::vector<std::ofstream> golden_path_vec;
    for (size_t level = 0; level < opt::levels; ++level){
        golden_path_vec.emplace_back(std::ofstream(opt::prefix_file + "_golden_path_" + std::to_string(level) + ".fa"));
    }

    double sTime = omp_get_wtime();
    std::cerr << "allocating bit vector" << std::endl;

    size_t filter_size = MIBloomFilter<uint32_t>::calcOptimalSize(opt::genome_size,	1, opt::occupancy);
    MIBFConstructSupport<uint32_t, multiLensfrHashIterator> miBFCS(opt::genome_size, opt::kmer_size,
            seed_string_vec.size(), opt::occupancy, filter_size, seed_string_vec);
    
    std::cerr << "finished allocating bit vector" << std::endl;
    std::cerr << "in "<< setprecision(4) << fixed
	          << omp_get_wtime() - sTime << "\n";


    
    std::cerr << "opening: " << opt::input << std::endl;

    uint32_t ids_inserted = 0;
    uint64_t bases = 0;


    std::cerr << "inserting bit vector" << std::endl;

    sTime = omp_get_wtime();
    {
        btllib::SeqReader reader(opt::input, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
        for (const auto record : reader) {
            if (record.seq.size() < min_seq_len) {
                continue;
            }
            multiLensfrHashIterator itr(record.seq, seed_string_vec); 
            miBFCS.insertBV(itr);
        }
    }

    std::cerr << "finished inserting bit vector" << std::endl;
    std::cerr << "in "<< setprecision(4) << fixed
	          << omp_get_wtime() - sTime << "\n";



    // setting up MIBF
    miBFCS.setup();
    std::vector<std::unique_ptr<MIBloomFilter<uint32_t>>> mibf_vec;
    for (size_t i = 0; i < opt::levels; ++i) {
        mibf_vec.emplace_back(std::move(std::unique_ptr<MIBloomFilter<uint32_t>>(miBFCS.getEmptyMIBF())));
    }


    std::cerr << "assigning tiles" << std::endl;
    sTime = omp_get_wtime();

    btllib::SeqReader reader(opt::input, btllib::SeqReader::Flag::LONG_MODE);
    for (const auto record : reader) {
        if (record.seq.size() < min_seq_len) {
                std::cerr << "too short" << std::endl;
                std::cerr << "skipping: " << record.id << std::endl;
                //wood_path <<  record.id << '\n' << record.seq <<  std::endl; skipping wood path output to reduce time
                ++id;
                continue;

        }
        if (id % 1000 == 0) {
            std::cerr << "processed " << id << " reads" <<std::endl;
        }
        size_t len = record.seq.size();
        size_t num_tiles = len / opt::tile_length;
        
        std::cerr << "name: " << record.id << std::endl;
        std::cerr << "num tiles: " << num_tiles - 2 << std::endl;

        
        bool assigned = true;

        //precompute hash to avoid rehashing
        std::vector<std::vector<uint64_t>> hashed_values(num_tiles, std::vector<uint64_t>()); 
#if _OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < num_tiles; ++i) {
            std::string tile_seq = record.seq.substr(i * opt::tile_length, opt::tile_length + opt::kmer_size - 1);
            multiLensfrHashIterator itr(tile_seq, seed_string_vec); 
            while(itr != itr.end()){
                for (size_t curr_hash = 0; curr_hash < seed_string_vec.size(); ++curr_hash) {
                    hashed_values[i].push_back((*itr)[curr_hash]);
                }
                ++itr;
            }
        }

        for (size_t level = 0; level < opt::levels; ++level) {
            auto& miBF = mibf_vec[level];
            std::cerr << "current level : " << level << std::endl;
            size_t num_tiles = hashed_values.size();
            std::vector<uint32_t> tiles_assigned_id_vec (num_tiles, 0);
            std::vector<bool> tiles_assigned_bool_vec (num_tiles, false);

            const size_t num_assigned_tiles = calc_num_assigned_tiles( miBF, hashed_values, tiles_assigned_id_vec, tiles_assigned_bool_vec);
            std::cerr << "num assigned tiles: " << num_assigned_tiles << std::endl;
            size_t num_unassigned_tiles = num_tiles - num_assigned_tiles;
            std::cerr << "num unassigned tiles: " << num_unassigned_tiles << std::endl;
            
            // assignment logic
            if (num_unassigned_tiles >= opt::unassigned_min && num_assigned_tiles <= opt::assigned_max) {
                assigned = false;
            }
            if (opt::second_pass) {
                if (num_assigned_tiles == num_tiles) {
                    assigned = true;
                } else {
                    assigned = false;
                }
            }

            if (!assigned) {
                std::cerr << "unassigned" << std::endl;
                ++ids_inserted;

#if _OPENMP
#pragma omp parallel for
#endif
                for (size_t i = 0; i < num_tiles; ++i) {
                    uint32_t curr_ids_inserted = ids_inserted + uint32_t( (i + 1) * opt::tile_length / 10000);
                    const auto& hashed_values_flat_array = hashed_values[i];
                    miBFCS.insertMIBF(*miBF, hashed_values_flat_array, curr_ids_inserted);//, non_singletons_bf_vec);
                    //miBFCS.insertSaturation(*miBF, Hhashes, ids_inserted); // don't care about saturation atm so skipping for speed
                        //}
                }
                ids_inserted = ids_inserted + uint32_t(record.seq.size() / 10000);
                //output read to golden path
                golden_path_vec[level] << ">" << record.id << '\n' << record.seq << std::endl;
                break; //breaks the level loop
                
            } else {
                if (num_assigned_tiles == num_tiles || opt::second_pass == true) {
                    std::cerr << "complete assignment" << std::endl;
                    continue;
                }
                /*if (num_tiles < 15) {
                    continue;
                }*/
                /*if (num_unassigned_tiles < 2) {
                    continue;
                }*/
                /*std::unordered_map<size_t, size_t> left_flank;
                for (int i = 1; i < 11; ++i ){
                    if (left_flank.contains(left_flank[tiles_assigned_id_vec[i]])) {
                        ++left_flank[tiles_assigned_id_vec[i]];
                    } else {
                        left_flank[tiles_assigned_id_vec[i]] == 1;
                    }
                }

                //size_t curr_max_id = 0;
                //size_t curr_max_hits = 0;
                std::vector<std::pair<size_t, size_t>> left_flank_vec;
                for ( const auto &[id, hits]: left_flank ) {
                    left_flank.push_back(std::make_pair(ids, hits));
                }                
                sort(vect.begin(), vect.end(), sort_by_sec);
                */
                size_t start_idx = 0;
                size_t end_idx = 0;
                ssize_t longest_start_idx = 0;
                ssize_t longest_end_idx = 0;
                size_t curr_stretch = 0;
                size_t longest_stretch = 0;
                std::cerr << "checkpoint 1" <<std::endl;
                for (size_t i = 1; i < num_tiles - 1; ++i){
                    if (tiles_assigned_bool_vec[i] == false && tiles_assigned_bool_vec[i - 1] == true) {
                        start_idx = i;
                        curr_stretch = 1;
                    } else if ((tiles_assigned_bool_vec[i] == false && tiles_assigned_bool_vec[i] == tiles_assigned_bool_vec[i-1]) && (i + 1 != num_tiles - 1)) {
                        ++curr_stretch;
                    } else if (tiles_assigned_bool_vec[i] == true && tiles_assigned_bool_vec[i] != tiles_assigned_bool_vec[i-1]) {
                        end_idx = i - 1;
                        if (longest_stretch < curr_stretch) {
                            longest_stretch = curr_stretch;
                            longest_start_idx = (ssize_t)start_idx;
                            longest_end_idx = (ssize_t)end_idx;
                        }                       
                    } else if (i + 1 == num_tiles - 1 && end_idx < start_idx) {
                        end_idx = i;
                        ++curr_stretch;
                        if (longest_stretch < curr_stretch) {
                            longest_stretch = curr_stretch;
                            longest_start_idx = (ssize_t)start_idx;
                            longest_end_idx = (ssize_t)end_idx;
                        }
                    }
                }

                bool good_flank = false;
                size_t trim_start_idx = longest_start_idx - 1;
                size_t trim_end_idx = longest_end_idx + 1;

                if (num_tiles < 15) {
                    bool good_right_flank = false;
                    bool good_left_flank = false;
                    std::unordered_map<size_t, size_t> left_flank;
                    for (ssize_t i = longest_start_idx - 1; i >= 1; --i ){
                        if (left_flank.contains(tiles_assigned_id_vec[i])) {
                            ++left_flank[tiles_assigned_id_vec[i]];
                        } else {
                            left_flank[tiles_assigned_id_vec[i]] = 1;
                        }
                    }
                    std::cerr << "checkpoint 2.1" <<std::endl;
                    std::vector<std::pair<size_t, size_t>> left_flank_vec;
                    for (  const auto &myPair: left_flank ) {
                        left_flank_vec.push_back(std::make_pair(myPair.first, myPair.second));
                    }
                    std::cerr << "checkpoint 2.2" <<std::endl;                
                    sort(left_flank_vec.begin(), left_flank_vec.end(), sort_by_sec);
                    std::cerr << "checkpoint 2.3" <<std::endl;
                    std::cerr << left_flank_vec.size() <<std::endl;
                    if (left_flank_vec.size() != 0) {
                        if (left_flank_vec[0].second >= 2) {
                            trim_start_idx = longest_start_idx - 2;
                            good_left_flank = true;
                        } else if ( left_flank_vec.size() >= 2 && (left_flank_vec[0].second + left_flank_vec[1].second> 3 && (left_flank_vec[0].first - left_flank_vec[1].first == 1 || left_flank_vec[1].first - left_flank_vec[0].first == 1))) {
                            trim_start_idx = longest_start_idx - 2;
                            good_left_flank = true;       

                        }
                    }
                    std::unordered_map<size_t, size_t> right_flank;
                    for (ssize_t i = longest_end_idx + 1; i < (ssize_t)num_tiles - 1; ++i ){
                        if (right_flank.contains(tiles_assigned_id_vec[i])) {
                            ++right_flank[tiles_assigned_id_vec[i]];
                        } else {
                            right_flank[tiles_assigned_id_vec[i]] = 1;
                        }
                    }
                    std::vector<std::pair<size_t, size_t>> right_flank_vec;
                    for ( const auto &myPair :  right_flank ) {
                        right_flank_vec.push_back(std::make_pair(myPair.first, myPair.second));
                    }                
                    sort(right_flank_vec.begin(), right_flank_vec.end(), sort_by_sec);
                    if (right_flank_vec.size() != 0) {
                        if (right_flank_vec[0].second >= 2) {
                            trim_end_idx = longest_end_idx + 2;
                            good_right_flank = true;
                        } else if ( right_flank_vec.size() >= 2 && (right_flank_vec[0].second + right_flank_vec[1].second> 3 && (right_flank_vec[0].first - right_flank_vec[1].first == 1 || right_flank_vec[1].first - right_flank_vec[0].first == 1))) {
                            trim_end_idx = longest_end_idx + 2;
                            good_right_flank = true;    

                        }
                    }
                    if (good_left_flank && good_right_flank) {
                        good_flank = true;
                    }
                } else {
                    std::cerr << "checkpoint 2" <<std::endl;
                    //bool valid = true;

                    if (longest_start_idx - 5 >= 1 ) {
                        std::unordered_map<size_t, size_t> left_flank;
                        for (ssize_t i = longest_start_idx - 5; i < longest_start_idx; ++i ){
                            if (left_flank.contains(tiles_assigned_id_vec[i])) {
                                ++left_flank[tiles_assigned_id_vec[i]];
                            } else {
                                left_flank[tiles_assigned_id_vec[i]] = 1;
                            }
                        }
                        std::cerr << "checkpoint 2.1" <<std::endl;
                        std::vector<std::pair<size_t, size_t>> left_flank_vec;
                        for (  const auto &myPair: left_flank ) {
                            left_flank_vec.push_back(std::make_pair(myPair.first, myPair.second));
                        }
                        std::cerr << "checkpoint 2.2" <<std::endl;                
                        sort(left_flank_vec.begin(), left_flank_vec.end(), sort_by_sec);
                        std::cerr << "checkpoint 2.3" <<std::endl;
                        std::cerr << left_flank_vec.size() <<std::endl;
                        if (left_flank_vec[0].second >= 2) {
                            trim_start_idx = longest_start_idx - 2;
                            good_flank = true;
                        } else if (left_flank_vec[0].second + left_flank_vec[1].second> 3 && (left_flank_vec[0].first - left_flank_vec[1].first == 1 || left_flank_vec[1].first - left_flank_vec[0].first == 1)) {
                            trim_start_idx = longest_start_idx - 2;
                            good_flank = true;       

                        }  
                        std::cerr << "checkpoint 2.4" <<std::endl;
                    } else {
                        trim_start_idx = 0;
                    }
                    std::cerr << "checkpoint 3" <<std::endl;

                    if (longest_end_idx + 5 < (ssize_t)num_tiles - 1 ) {
                        std::unordered_map<size_t, size_t> right_flank;
                        for (ssize_t i = longest_end_idx + 5; i > longest_end_idx; --i ){
                            if (right_flank.contains(tiles_assigned_id_vec[i])) {
                                ++right_flank[tiles_assigned_id_vec[i]];
                            } else {
                                right_flank[tiles_assigned_id_vec[i]] = 1;
                            }
                        }
                        std::vector<std::pair<size_t, size_t>> right_flank_vec;
                        for ( const auto &myPair :  right_flank ) {
                            right_flank_vec.push_back(std::make_pair(myPair.first, myPair.second));
                        }                
                        sort(right_flank_vec.begin(), right_flank_vec.end(), sort_by_sec);
                        if (right_flank_vec[0].second >= 2) {
                            trim_end_idx = longest_end_idx + 2;
                            good_flank = true;
                        } else if (right_flank_vec[0].second + right_flank_vec[1].second> 3 && (right_flank_vec[0].first - right_flank_vec[1].first == 1 || right_flank_vec[1].first - right_flank_vec[0].first == 1)) {
                        trim_end_idx = longest_end_idx + 2;
                        good_flank = true;    

                        }

                    } else {
                        trim_end_idx = (ssize_t)num_tiles - 1;
                    }
                }
                std::cerr << "checkpoint 4" <<std::endl;
                if (good_flank) {
                    assigned = false;
                    std::cerr << "trimmed" << std::endl;
                    ++ids_inserted;

#if _OPENMP
#pragma omp parallel for
#endif
                    for (size_t i = trim_start_idx; i <= trim_end_idx; ++i) {
                        uint32_t curr_ids_inserted = ids_inserted + uint32_t( (i - trim_start_idx + 1) * opt::tile_length / 10000);
                        const auto& hashed_values_flat_array = hashed_values[i];
                        miBFCS.insertMIBF(*miBF, hashed_values_flat_array, curr_ids_inserted);//, non_singletons_bf_vec);
                        //miBFCS.insertSaturation(*miBF, Hhashes, ids_inserted); // don't care about saturation atm so skipping for speed
                            //}
                    }
                    ids_inserted = ids_inserted + uint32_t((trim_end_idx - trim_start_idx) * 1000 / 10000);
                    //output read to golden path
                    if (trim_end_idx == num_tiles -1) {
                        golden_path_vec[level] << ">trimmed" << record.id << '\n' << record.seq.substr(trim_start_idx * 1000,std::string::npos) << std::endl;
                    } else {
                        golden_path_vec[level] << ">trimmed" << record.id << '\n' << record.seq.substr(trim_start_idx * 1000, (trim_end_idx - trim_start_idx + 1) * 1000) << std::endl;
                    }
                    

                    break; //breaks the level loop
                    
                }

            }

        }
        if (assigned) {
            std::cerr << "assigned" << std::endl;
            //output read to wood path
            //wood_path <<  record.id << '\n' << record.seq <<  std::endl; skipping wood path output to reduce time
            
        }
        //tracking number of bases inserted at an id
        std::cerr << "inserted: " << id << " ";
        if (assigned) {
            std::cerr << bases << std::endl;
        } else {
            bases += record.seq.size();
            std::cerr << bases << std::endl;
        }
        
        ++id;
    }
    std::cerr << "assgined" << std::endl;
        	std::cerr << "in "<< setprecision(4) << fixed
	          << omp_get_wtime() - sTime << "\n";
}