#include "opt.hpp"
#include "read_hashing.hpp"
#include "spaced_seeds.hpp"

#include "Common/Options.h"
#include "MIBFConstructSupport.hpp"
#include "MIBloomFilter.hpp"
#include "multiLensfrHashIterator.hpp"
#include "calc_phred_average.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/util.hpp"
#include "ntcard.hpp"
#include <sdsl/int_vector.hpp>

#if _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <tuple>
#include <set>
#include <vector>


void log_tile_states(std::vector<uint32_t>& tiles_assigned_id_vec, std::vector<uint8_t>& tiles_assigned_bool_vec) {
    //print which id each tile is assigned to
    for (const auto& tiles_assigned_id : tiles_assigned_id_vec) {
        std::cerr << tiles_assigned_id << "\t";
    }

    std::cerr << std::endl;
    //print whether the id to each tile is assigned
    for (const auto& tiles_assigned_bool : tiles_assigned_bool_vec) {
        std::cerr << +tiles_assigned_bool << "\t";
    }
    std::cerr << std::endl;
}


bool
sort_by_sec(const pair<size_t, size_t>& a, const pair<size_t, size_t>& b)
{
  return (a.second > b.second);
}

void
fill_bit_vector(const std::string& input_file,
                MIBFConstructSupport<uint32_t, multiLensfrHashIterator>& miBFCS,
                const size_t min_seq_len,
                const std::vector<std::string>& spaced_seeds,
                std::unordered_set<std::string>& filter_out_reads)
{
  std::cerr << "inserting bit vector" << std::endl;

  auto sTime = omp_get_wtime();

  btllib::SeqReader reader(input_file, btllib::SeqReader::Flag::LONG_MODE);
    if (reader.get_format() != btllib::SeqReader::Format::FASTQ) {
    std::cerr << "Gold Path requires fastq format" << std::endl;
    exit(1);
    }
#pragma omp parallel
  for (const auto record : reader) {
    if (record.seq.size() < min_seq_len) {
      continue;
    }
    const auto phred_stat = calc_phred_average(record.qual);
    if ( phred_stat.first < opt::phred_min || phred_stat.second >= 5) {
#pragma omp critical
{
        filter_out_reads.insert(record.id);
}
        continue;
    }
    multiLensfrHashIterator itr(record.seq, spaced_seeds);
    miBFCS.insertBV(itr);
  }

  std::cerr << "finished inserting bit vector" << std::endl;
  std::cerr << "in " << setprecision(4) << fixed << omp_get_wtime() - sTime
            << "\n";
}

size_t
calc_num_assigned_tiles(const MIBloomFilter<uint32_t>& miBF,
                        const std::vector<std::vector<uint64_t>> hashed_values,
                        std::vector<uint32_t>& tiles_assigned_id_vec,
                        std::vector<uint8_t>& tiles_assigned_bool_vec)
{

    size_t num_assigned_tiles = 0;
    size_t num_tiles = hashed_values.size();
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> tiles_assigned_all_id_vec(hashed_values.size());


#if _OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < num_tiles; ++i) {

        std::map<uint32_t, std::pair<uint32_t, uint32_t>> id_counts; // store counts of ids

        // Reusable vector for ranks 
        vector<uint64_t> m_rank_pos(miBF.getHashNum());
        // Reusable vector for IDs 
        vector<uint32_t> m_data(miBF.getHashNum());
        
        const auto& hashed_values_flat_array = hashed_values[i];
        std::vector<uint64_t> hashes(opt::hash_num, 0);
        for (size_t curr_frame = 0; curr_frame < (hashed_values_flat_array.size() / opt::hash_num); ++curr_frame) {
            
            for (size_t curr_hash = 0; curr_hash < opt::hash_num; ++curr_hash) {
                hashes[curr_hash] = hashed_values_flat_array[curr_frame * opt::hash_num + curr_hash];
            }

            std::set<uint32_t> unique_ids;
            if (miBF.atRank(hashes, m_rank_pos)) {                      // if its a hit
                m_data = miBF.getData(m_rank_pos);                 // m_data has ID's
                for(unsigned m = 0; m < miBF.getHashNum(); m++){   // iterate over ID's
                    if(m_data[m] > miBF.s_mask){                     // if ID is saturated
                        uint32_t new_id = m_data[m] & miBF.s_antiMask;
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
        std::vector<std::pair<uint32_t, uint32_t>> id_counts_vec;
        for (auto id_counts_it = id_counts.begin(); id_counts_it != id_counts.end(); ++id_counts_it) { //find id with highest count in a tile
            if (id_counts_it->second.second > curr_id_count) {
                curr_id = id_counts_it->first;
                curr_id_count = id_counts_it->second.second;
            }
            if (id_counts_it->second.second > 2) {
                id_counts_vec.emplace_back(std::make_pair(id_counts_it->first, id_counts_it->second.second));
            }
        }

        sort(id_counts_vec.begin(), id_counts_vec.end(), sort_by_sec);

        tiles_assigned_id_vec[i] = curr_id;
        tiles_assigned_all_id_vec[i] = id_counts_vec;

        
    }

    for (size_t i = 0; i < num_tiles; ++i) {
        if (!tiles_assigned_all_id_vec[i].empty()){
            if (tiles_assigned_all_id_vec[i][0].second > opt::threshold){
                tiles_assigned_bool_vec[i] = 1;
            }
        }
    }
if (num_tiles >= 3) {
if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}
    num_assigned_tiles = 0;
    for (const auto& is_tile_assigned : tiles_assigned_bool_vec) {
        if (is_tile_assigned) {
            ++num_assigned_tiles;
        }
    }

    for (size_t i = 1; i < num_tiles; ++i) {
        uint32_t curr_id = tiles_assigned_id_vec[i];
        uint32_t prev_id = tiles_assigned_id_vec[i - 1];
        if (curr_id != prev_id) {
            for (auto& pair : tiles_assigned_all_id_vec[i]) {
                if (pair.first == prev_id) {
                    tiles_assigned_id_vec[i] = prev_id;
                    if (pair.second > opt::threshold) {
                        tiles_assigned_bool_vec[i] = 1;
                    } else {
                        tiles_assigned_bool_vec[i] = 0;
                    }
                }
            }
        }
    }

if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}

    for (ssize_t i = (ssize_t)num_tiles - 2; i >= 0; --i) {
        uint32_t curr_id = tiles_assigned_id_vec[i];
        uint32_t prev_id = tiles_assigned_id_vec[i + 1];
        if (curr_id != prev_id) {
            for (auto& pair : tiles_assigned_all_id_vec[i]) {
                if (pair.first == prev_id) {
                    tiles_assigned_id_vec[i] = prev_id;
                    if (pair.second > opt::threshold) {
                        tiles_assigned_bool_vec[i] = 1;
                    } else {
                        tiles_assigned_bool_vec[i] = 0;
                    }
                }
            }
        }
    }

if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}



    for (size_t i = 1; i < num_tiles -1; ++i) {
        auto& curr_assign = tiles_assigned_bool_vec[i];
        auto& curr_id = tiles_assigned_id_vec[i];
        const auto& prev_assign = tiles_assigned_bool_vec[i - 1];
        const auto& prev_id = tiles_assigned_id_vec[i - 1];
        const auto& next_assign = tiles_assigned_bool_vec[i + 1];
        const auto& next_id = tiles_assigned_id_vec[i + 1];
        if (!curr_assign) {
            if ((curr_id == prev_id && prev_assign) || (curr_id == next_id && next_assign)) {
                curr_assign = 1;
            } else if ((curr_id == prev_id + 1 && prev_assign) || (curr_id == next_id + 1 && next_assign)) {
                curr_assign = 1;
            } else if ((curr_id == prev_id - 1 && prev_assign ) || (curr_id == next_id - 1 && next_assign )) {
                curr_assign = 1;
            } else if (prev_id == next_id  && prev_assign  &&  next_assign ) {
                tiles_assigned_bool_vec[i] = prev_assign;
                curr_id = prev_id;
            }
        }
    }

    for (size_t i = num_tiles - 2; i >= 1; --i) {
        auto& curr_assign = tiles_assigned_bool_vec[i];
        auto& curr_id = tiles_assigned_id_vec[i];
        const auto& prev_assign = tiles_assigned_bool_vec[i - 1];
        const auto& prev_id = tiles_assigned_id_vec[i - 1];
        const auto& next_assign = tiles_assigned_bool_vec[i + 1];
        const auto& next_id = tiles_assigned_id_vec[i + 1];
        if (!curr_assign) {
            if ((curr_id == prev_id && prev_assign ) || (curr_id == next_id && next_assign )) {
                curr_assign = 1;
            } else if ((curr_id == prev_id + 1 && prev_assign ) || (curr_id == next_id + 1 && next_assign )) {
                curr_assign = 1;
            } else if ((curr_id == prev_id - 1 && prev_assign ) || (curr_id == next_id - 1 && next_assign )) {
                curr_assign = 1;
            } else if (prev_id == next_id  && prev_assign  &&  next_assign ) {
                tiles_assigned_bool_vec[i] = prev_assign;
                curr_id = prev_id;
            }
        }
    }

if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}
    size_t start_idx = 0;
    size_t end_idx = 0;
    //size_t curr_stretch = 0;
    std::vector<std::pair<size_t, size_t>> coord_vec;
    for (size_t i = 1; i < num_tiles - 1; ++i){
        const auto& curr_assign = tiles_assigned_bool_vec[i];
        const auto& prev_assign = tiles_assigned_bool_vec[i - 1];
        if (!curr_assign && prev_assign ) {
            start_idx = i;
        } else if (curr_assign  && !prev_assign) {
            end_idx = i - 1;
            coord_vec.push_back(std::make_pair(start_idx, end_idx));
        }
    }
    for (const auto& coords : coord_vec) {
        if (coords.first== 0 || coords.second == num_tiles - 1 ) {
            continue;
        } 

        const auto& left  = tiles_assigned_id_vec[coords.first - 1];
        const auto& right  = tiles_assigned_id_vec[coords.second + 1];
        if (left == right|| left == right + 1 || left == right - 1) {
            for (auto i = coords.first; i <= coords.second; ++i) {
                tiles_assigned_bool_vec[i] = 1;
                tiles_assigned_id_vec[i] = left;
            }
        }

    }

if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}
    if (num_tiles >= 3) {
        for (size_t i = 2; i < num_tiles -2; ++i) {
            const auto& curr_assign = tiles_assigned_bool_vec[i];
            const auto& prev_assign = tiles_assigned_bool_vec[i - 1];
            const auto& next_assign = tiles_assigned_bool_vec[i + 1];
            if (curr_assign ) {
                if (!prev_assign && !next_assign) {
                    tiles_assigned_bool_vec[i] = 0;
                }
            }
        }

        for (size_t i = num_tiles - 3; i >= 2; --i) {
            const auto& curr_assign = tiles_assigned_bool_vec[i];
            const auto& prev_assign = tiles_assigned_bool_vec[i - 1];
            const auto& next_assign = tiles_assigned_bool_vec[i + 1];
            if (curr_assign ) {
                if (!prev_assign && !next_assign) {
                    tiles_assigned_bool_vec[i] = 0;
                }
            }
        }
    }


if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}

    std::map<uint32_t, std::vector<uint32_t>> id_to_idx;
    for (size_t i = 0; i < num_tiles; ++i) {
        uint32_t curr_id = tiles_assigned_id_vec[i];
        uint32_t curr_id_bool = tiles_assigned_bool_vec[i];
        if (curr_id_bool) {
            id_to_idx[curr_id].emplace_back(i);
        }
    }


    for (auto& id_idx_vec_pair : id_to_idx) {
        auto& idx_vec = id_idx_vec_pair.second;
        sort(idx_vec.begin(), idx_vec.end());
        
        for (size_t i = 1; i < idx_vec.size(); ++i){
            uint32_t prev_idx = idx_vec[i - 1];
            uint32_t curr_idx = idx_vec[i];
            if (curr_idx > prev_idx + 1) {
                uint32_t prev_id = tiles_assigned_id_vec[prev_idx];
                for (size_t j = prev_idx + 1; j <= curr_idx; ++j) {
                     tiles_assigned_id_vec[j] = prev_id;
                }
            }
        }
    }
if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}

    size_t last_id = tiles_assigned_id_vec[num_tiles - 1];
    size_t second_last_id = tiles_assigned_id_vec[num_tiles - 2];
    size_t start_id = tiles_assigned_id_vec[0];
    size_t second_start_id = tiles_assigned_id_vec[1];
    if (last_id == second_last_id || last_id == second_last_id + 1 || last_id == second_last_id - 1) {
        tiles_assigned_bool_vec[num_tiles - 1] = 1;
    }
    if (start_id == second_start_id || start_id == second_start_id + 1 || start_id == second_start_id - 1) {
        tiles_assigned_bool_vec[0] = 1;
    }


    for (size_t i = 1; i < num_tiles - 1; ++i) {
        auto& curr_assign = tiles_assigned_bool_vec[i];
        const auto& curr_id = tiles_assigned_id_vec[i];
        const auto& prev_id = tiles_assigned_id_vec[i - 1];
        const auto& next_id = tiles_assigned_id_vec[i + 1];
        if (curr_id != next_id && curr_id != next_id  - 1 && curr_id != next_id  + 1 && curr_id != prev_id && curr_id != prev_id  - 1 && curr_id != prev_id  + 1) {
            curr_assign = 0;
        }
        
    }


if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}

    start_idx = 0;
    end_idx = 0;

    coord_vec.clear();
    for (size_t i = 1; i < num_tiles - 1; ++i){
        const auto& curr_assign = tiles_assigned_bool_vec[i];
        const auto& prev_assign = tiles_assigned_bool_vec[i - 1];
        if (curr_assign  && !prev_assign) {
            start_idx = i;
        } else if (!curr_assign && prev_assign  ) {
            end_idx = i - 1;
            coord_vec.push_back(std::make_pair(start_idx, end_idx));
        }
    }

    for (const auto& coords : coord_vec) {
        if (coords.second - coords.first + 1 <= 5 ) {
            for (auto i = coords.first; i <= coords.second; ++i) {
                tiles_assigned_bool_vec[i] = 0;
            }
        }

    }


if (verbose) {
    log_tile_states(tiles_assigned_id_vec, tiles_assigned_bool_vec);
}
}
    num_assigned_tiles = 0;
    for (const auto& is_tile_assigned : tiles_assigned_bool_vec) {
        if (is_tile_assigned) {
            ++num_assigned_tiles;
        }
    }
    return num_assigned_tiles;
}

inline void
process_read(const btllib::SeqReader::Record& record,
             const std::vector<std::vector<uint64_t>>& hashed_values,
             std::vector<std::ofstream>& golden_path_vec,
             std::vector<std::unique_ptr<MIBloomFilter<uint32_t>>>& mibf_vec,
             MIBFConstructSupport<uint32_t, multiLensfrHashIterator>& miBFCS,
             uint64_t& inserted_bases,
             uint64_t& target_bases,
             uint64_t& curr_path,
             uint32_t& id,
             uint32_t& ids_inserted,
             const size_t min_seq_len,
             const std::unordered_set<std::string>& filter_out_reads)
{
  if (record.seq.size() < min_seq_len) {
    if (verbose) {
      std::cerr << "too short" << std::endl;
      std::cerr << "skipping: " << record.id << std::endl;
    }

    ++id;
    return;
    }
    if (!filter_out_reads.empty()) {
        if (filter_out_reads.find(record.id) != filter_out_reads.end()) {
    if (verbose) {
      std::cerr << "hairpin or quality too low" << std::endl;
      std::cerr << "skipping: " << record.id << std::endl;
    }

    ++id;
    return;
        }
    }
  if (id % 10000 == 0) {
    std::cerr << "processed " << id << " reads" << std::endl;
  }
  size_t len = record.seq.size();
  size_t num_tiles = len / opt::tile_length;

  if (verbose) {
    std::cerr << "name: " << record.id << std::endl;
    std::cerr << "num tiles: " << num_tiles << std::endl;
 }

  bool assigned = true;

  auto& miBF = mibf_vec[0];


  std::vector<uint32_t> tiles_assigned_id_vec(num_tiles, 0);
  std::vector<uint8_t> tiles_assigned_bool_vec(num_tiles, 0);
  const size_t num_assigned_tiles = calc_num_assigned_tiles(
  *miBF, hashed_values, tiles_assigned_id_vec, tiles_assigned_bool_vec);
  if (verbose) {
    std::cerr << "num assigned tiles: " << num_assigned_tiles << std::endl;
  }
  const size_t num_unassigned_tiles = num_tiles - num_assigned_tiles;
  if (verbose) {
    std::cerr << "num unassigned tiles: " << num_unassigned_tiles << std::endl;
   }

   // assignment logic
    if (num_unassigned_tiles >= opt::unassigned_min &&
        num_assigned_tiles <= opt::assigned_max) {
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
if (verbose) {
                std::cerr << "unassigned" << std::endl;
}
                ++ids_inserted;
                size_t block_start = 0;
                while (block_start < num_tiles) {
                    size_t block_end = std::min(block_start + opt::block_size, num_tiles);
                    uint32_t curr_ids_inserted = ids_inserted + uint32_t( (block_start ) / opt::block_size);
                    std::vector<uint64_t> hashed_values_new_array;
                    hashed_values_new_array.reserve(hashed_values[block_start].size() * opt::block_size);
                    for (size_t i = block_start; i < block_end; ++i) {
                        hashed_values_new_array.insert(hashed_values_new_array.end(), hashed_values[i].begin(), hashed_values[i].end());
                    }
                    miBFCS.insertMIBF(*miBF, hashed_values_new_array, curr_ids_inserted);
                    block_start = block_start + opt::block_size;
                }
                ids_inserted = ids_inserted + uint32_t(record.seq.size() / (opt::tile_length * opt::block_size));
                //output read to golden path
                const auto phred_stats = calc_phred_average(record.qual);
                golden_path_vec[0] << "@" << record.id << '_' << phred_stats.first << '_' << phred_stats.second << '\n' << record.seq << '\n' << '+' << '\n' << record.qual << std::endl;
                
                inserted_bases += record.seq.size();
                if (opt::temp_mode || opt::new_temp_mode) {
                    if (target_bases < inserted_bases) {
                        ++curr_path;
                        if (opt::max_paths < curr_path) {
                            exit(0);
                        }
                        inserted_bases = 0;
                        mibf_vec.pop_back();
                        mibf_vec.emplace_back(std::move(std::unique_ptr<MIBloomFilter<uint32_t>>(miBFCS.getEmptyMIBF())));
                        golden_path_vec.pop_back();
                        golden_path_vec.emplace_back(std::ofstream(opt::prefix_file + "_golden_path_" + std::to_string(curr_path) + ".fq"));
                        ids_inserted = 0;
                    }
                }
                
            } else {
                if (opt::temp_mode) {
                    return;
                }
                if (num_assigned_tiles == num_tiles || opt::second_pass == true) {
if (verbose) {
                    std::cerr << "complete assignment" << std::endl;
}
                    return;
                }
                size_t start_idx = 0;
                size_t end_idx = 0;
                ssize_t longest_start_idx = 0;
                ssize_t longest_end_idx = 0;
                size_t curr_stretch = 0;
                size_t longest_stretch = 0;
                for (size_t i = 1; i < num_tiles - 1; ++i){
                    if (!tiles_assigned_bool_vec[i] && tiles_assigned_bool_vec[i - 1] ) {
                        start_idx = i;
                        curr_stretch = 1;
                    } else if ((!tiles_assigned_bool_vec[i] && tiles_assigned_bool_vec[i] == tiles_assigned_bool_vec[i-1]) && (i + 1 != num_tiles - 1)) {
                        ++curr_stretch;
                    } else if (tiles_assigned_bool_vec[i]  && tiles_assigned_bool_vec[i] != tiles_assigned_bool_vec[i-1]) {
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
                size_t trim_start_idx = 0;
                if (longest_start_idx != 0) {
                    trim_start_idx = longest_start_idx - 1;
                } else {
                    trim_start_idx = longest_start_idx;
                }
                size_t trim_end_idx = longest_end_idx + 1;

                if (num_tiles < 15) {
                    bool good_right_flank = false;
                    bool good_left_flank = false;
                    std::map<size_t, size_t> left_flank;
                    for (ssize_t i = longest_start_idx - 1; i >= 0; --i ){
                        if (left_flank.find(tiles_assigned_id_vec[i]) != left_flank.end()) {
                            ++left_flank[tiles_assigned_id_vec[i]];
                        } else {
                            left_flank[tiles_assigned_id_vec[i]] = 1;
                        }
                    }
                    std::vector<std::pair<size_t, size_t>> left_flank_vec;
                    for (  const auto &myPair: left_flank ) {
                        left_flank_vec.push_back(std::make_pair(myPair.first, myPair.second));
                    }              
                    sort(left_flank_vec.begin(), left_flank_vec.end(), sort_by_sec);
                    if (left_flank_vec.size() != 0) {
                        if (left_flank_vec[0].second >= 2) {
                            if (longest_start_idx != 0) {
                                trim_start_idx = longest_start_idx - 1;
                            } else {
                                trim_start_idx = longest_start_idx;
                            }
                            good_left_flank = true;
                        } else if ( left_flank_vec.size() >= 2 && (left_flank_vec[0].second + left_flank_vec[1].second> 3 && (left_flank_vec[0].first - left_flank_vec[1].first == 1 || left_flank_vec[1].first - left_flank_vec[0].first == 1))) {
                            if (longest_start_idx != 0) {
                                trim_start_idx = longest_start_idx - 1;
                            } else {
                                trim_start_idx = longest_start_idx;
                            }
                            good_left_flank = true;       

                        } else if (trim_start_idx == 0) {
                            good_left_flank = true; 
                        }
                    }
                    if (trim_start_idx == 0) {
                            good_left_flank = true; 
                    }
                    std::map<size_t, size_t> right_flank;
                    for (ssize_t i = longest_end_idx + 1; i < (ssize_t)num_tiles; ++i ){
                        if (right_flank.find(tiles_assigned_id_vec[i]) != right_flank.end()) {
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
                            trim_end_idx = longest_end_idx + 1;
                            good_right_flank = true;
                        } else if ( right_flank_vec.size() >= 2 && (right_flank_vec[0].second + right_flank_vec[1].second> 3 && (right_flank_vec[0].first - right_flank_vec[1].first == 1 || right_flank_vec[1].first - right_flank_vec[0].first == 1))) {
                            trim_end_idx = longest_end_idx + 1;
                            good_right_flank = true;    

                        } else if (trim_end_idx == num_tiles - 1) {
                            good_right_flank = true; 
                        }
                    }
                    if (trim_end_idx == num_tiles - 1) {
                            good_right_flank = true; 
                    }

                    if (good_left_flank && good_right_flank) {
                        good_flank = true;
                    }


                } else {


                    if (longest_start_idx - 5 >= 1 ) {
                        std::map<size_t, size_t> left_flank;
                        for (ssize_t i = longest_start_idx - 5; i < longest_start_idx; ++i ){
                            if (left_flank.find(tiles_assigned_id_vec[i]) != left_flank.end()) {
                                ++left_flank[tiles_assigned_id_vec[i]];
                            } else {
                                left_flank[tiles_assigned_id_vec[i]] = 1;
                            }
                        }

                        std::vector<std::pair<size_t, size_t>> left_flank_vec;
                        for (  const auto &myPair: left_flank ) {
                            left_flank_vec.push_back(std::make_pair(myPair.first, myPair.second));
                        }
               
                        sort(left_flank_vec.begin(), left_flank_vec.end(), sort_by_sec);


                        if (left_flank_vec[0].second >= 2) {
                            if (longest_start_idx != 0) {
                                trim_start_idx = longest_start_idx - 1;
                            } else {
                                trim_start_idx = longest_start_idx;
                            }
                            good_flank = true;
                        } else if (left_flank_vec[0].second + left_flank_vec[1].second> 3 && (left_flank_vec[0].first - left_flank_vec[1].first == 1 || left_flank_vec[1].first - left_flank_vec[0].first == 1)) {
                            if (longest_start_idx != 0) {
                                trim_start_idx = longest_start_idx - 1;
                            } else {
                                trim_start_idx = longest_start_idx;
                            }
                            good_flank = true;       

                        }  

                    } else {

                        good_flank = true;
                        trim_start_idx = 0;
                    }


                    if (longest_end_idx + 5 < (ssize_t)num_tiles - 1 ) {
                        std::map<size_t, size_t> right_flank;
                        for (ssize_t i = longest_end_idx + 5; i > longest_end_idx; --i ){
                            if (right_flank.find(tiles_assigned_id_vec[i]) != right_flank.end()) {
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
                            trim_end_idx = longest_end_idx + 1;
                            good_flank = true;
                        } else if (right_flank_vec[0].second + right_flank_vec[1].second> 3 && (right_flank_vec[0].first - right_flank_vec[1].first == 1 || right_flank_vec[1].first - right_flank_vec[0].first == 1)) {
                        trim_end_idx = longest_end_idx + 1;
                        good_flank = true;    

                        }

                    } else {
                        good_flank = true;
                        trim_end_idx = (ssize_t)num_tiles - 1;
                    }
                }
                if (good_flank) {
                    assigned = false;
if (verbose) {
                    std::cerr << "trimmed" << std::endl;
}
                    ++ids_inserted;
                    size_t block_start = trim_start_idx;
                    while (block_start <= trim_end_idx) {
                        size_t block_end = std::min(block_start + opt::block_size - 1, trim_end_idx);
                        uint32_t curr_ids_inserted = ids_inserted + uint32_t( (block_start - trim_start_idx + 1)  / opt::block_size);
                        std::vector<uint64_t> hashed_values_new_array;
                        hashed_values_new_array.reserve(hashed_values[block_start].size() * opt::block_size);
                        for (size_t i = block_start; i <= block_end; ++i) {
                            hashed_values_new_array.insert(hashed_values_new_array.end(), hashed_values[i].begin(), hashed_values[i].end());
                        }
                        miBFCS.insertMIBF(*miBF, hashed_values_new_array, curr_ids_inserted);
                        block_start = block_start + opt::block_size;
                    }
                    ids_inserted = ids_inserted + uint32_t((trim_end_idx - trim_start_idx)  / opt::block_size);
                    //output read to golden path
                    if (trim_end_idx == num_tiles - 1) {
                        std::string new_seq = record.seq.substr(trim_start_idx * opt::tile_length, std::string::npos);
                        std::string new_qual = record.qual.substr(trim_start_idx * opt::tile_length, std::string::npos);
                        const auto phred_stats = calc_phred_average(new_qual);
                        golden_path_vec[0] << "@trimmed" << record.id << '_' << phred_stats.first << '_' << phred_stats.second << '\n' << new_seq << '\n' << '+' << '\n' << new_qual << std::endl;
                        inserted_bases += new_seq.size();
                    } else {
                        std::string new_seq = record.seq.substr(trim_start_idx * opt::tile_length, (trim_end_idx - trim_start_idx + 1) * opt::tile_length);
                        std::string new_qual = record.qual.substr(trim_start_idx * opt::tile_length, (trim_end_idx - trim_start_idx + 1) * opt::tile_length);
                        const auto phred_stats = calc_phred_average(new_qual);
                        golden_path_vec[0] << "@trimmed" << record.id << '_' << phred_stats.first << '_' << phred_stats.second << '\n' << new_seq << '\n' << '+' << '\n' << new_qual << std::endl;
                        inserted_bases += new_seq.size();
                    }
                    
                    if (opt::temp_mode || opt::new_temp_mode) {
                        if (target_bases < inserted_bases) {
                            ++curr_path;
                            if (opt::max_paths < curr_path) {
                                exit(0);
                            }
                            inserted_bases = 0;
                            mibf_vec.pop_back();
                            mibf_vec.emplace_back(std::move(std::unique_ptr<MIBloomFilter<uint32_t>>(miBFCS.getEmptyMIBF())));
                            golden_path_vec.pop_back();
                            golden_path_vec.emplace_back(std::ofstream(opt::prefix_file + "_golden_path_" + std::to_string(curr_path) + ".fq"));
                            ids_inserted = 0;
                        }
                    }
                    
                    
                }

            }


  
  if (assigned) {
    if (verbose) {
      std::cerr << "assigned" << std::endl;
    }
    // output read to wood path
   
  }
  ++id;
}

int
main(int argc, char** argv)
{
  process_options(argc, argv);

#if _OPENMP
  omp_set_num_threads(opt::jobs);
#endif

  if (opt::second_pass) {
    std::cerr << "second_pass" << std::endl;
  }

  // srand (1); // for testing, change to srand(time(NULL)) for actual code
  const auto seed_string_vec = make_seed_pattern(
    opt::seed_preset, opt::kmer_size, opt::weight, opt::hash_num);

  if (opt::genome_size == 0) {
    if (opt::ntcard) {
      opt::genome_size = calc_ntcard_genome_size(
        opt::input, opt::kmer_size, seed_string_vec, opt::jobs);
    } else {
      static const uint8_t BASES = 4;
      static const float HASH_UNIVERSE_COEFFICIENT = 0.5;
      opt::genome_size =
        pow(BASES, opt::weight) * HASH_UNIVERSE_COEFFICIENT * opt::hash_num;
    }
  }

  std::cerr << "Calculating " << opt::levels << " golden path(s)"
            << "\n"
            << "Using:"
            << "\n"
            << "tile length: " << opt::tile_length << "\n"
            << "block size: " << opt::block_size << "\n"
            << "seed pattens: " << opt::hash_num << "\n"
            << "threshold: " << opt::threshold << "\n"
            << "base seed pattern: " << seed_string_vec[0] << "\n"
            << "minimum unassigned tiles: " << opt::unassigned_min << "\n"
            << "maximum assigned tiles: " << opt::assigned_max << "\n"
            << "genome size: " << opt::genome_size << "\n"
            << "minimum average phred quality score: " << opt::phred_min << "\n"
            << "occupancy: " << opt::occupancy << "\n"
            << "jobs: " << opt::jobs << std::endl;

  std::unordered_set<std::string> filter_out_reads;
  if (opt::filter_file != "") {
    std::cerr << "Using only reads not found in: " << opt::filter_file << std::endl;
    std::string read_name = "";
    std::ifstream infileStream(opt::filter_file);
    while (infileStream >> read_name) {
      filter_out_reads.insert(read_name);
    }
  }       
    

  std::vector<std::ofstream> golden_path_vec;

  golden_path_vec.emplace_back(std::ofstream(
      opt::prefix_file + "_golden_path_1.fq"));

  double sTime = omp_get_wtime();
  std::cerr << "allocating bit vector" << std::endl;

  size_t filter_size = MIBloomFilter<uint32_t>::calcOptimalSize(
    opt::genome_size, 1, opt::occupancy);
  MIBFConstructSupport<uint32_t, multiLensfrHashIterator> miBFCS(
    opt::genome_size,
    opt::kmer_size,
    seed_string_vec.size(),
    opt::occupancy,
    filter_size,
    seed_string_vec);

  std::cerr << "finished allocating bit vector" << std::endl;
  std::cerr << "in " << setprecision(4) << fixed << omp_get_wtime() - sTime
            << "\n";

  std::cerr << "opening: " << opt::input << std::endl;

  fill_bit_vector(opt::input, miBFCS, opt::min_length, seed_string_vec, filter_out_reads);

  // setting up MIBF
  miBFCS.setup();
  std::vector<std::unique_ptr<MIBloomFilter<uint32_t>>> mibf_vec;
    mibf_vec.emplace_back(miBFCS.getEmptyMIBF());

  std::cerr << "assigning tiles" << std::endl;
  sTime = omp_get_wtime();

  btllib::OrderQueueMPMC<ReadTileHashes> precomputed_hash_queue(
    btllib::SeqReader::LONG_MODE_BUFFER_SIZE,
    btllib::SeqReader::LONG_MODE_BLOCK_SIZE);
  start_read_hashing(opt::input,
                     opt::tile_length,
                     opt::min_length,
                     opt::kmer_size,
                     seed_string_vec,
                     precomputed_hash_queue,
                     6,
                     filter_out_reads);

  uint64_t inserted_bases = 0;
  uint64_t target_bases = opt::ratio * opt::target_size;
  uint64_t curr_path = 1;
  uint32_t id = 1;
  uint32_t ids_inserted = 0;
  //std::unordered_map<uint32_t, uint8_t> id_to_num_tiles_inserted;
  while (true) {
    decltype(precomputed_hash_queue)::Block block(
      btllib::SeqReader::LONG_MODE_BLOCK_SIZE);

    precomputed_hash_queue.read(block);

    if (block.count == 0) {
      break;
    }
    for (size_t idx = 0; idx < block.count; idx++) {
      const auto& read_hashes = block.data[idx];
      const auto& record = read_hashes.read;
      const auto& hashed_values = read_hashes.tile_hashes;
      process_read(record,
                   hashed_values,
                   golden_path_vec,
                   mibf_vec,
                   miBFCS,
                   inserted_bases,
                   target_bases,
                   curr_path,
                   id,
                   ids_inserted,
                   opt::min_length,
                   filter_out_reads);
    }
  }

  std::cerr << "assigned" << std::endl;
  std::cerr << "in " << setprecision(4) << fixed << omp_get_wtime() - sTime
            << "\n";
}
