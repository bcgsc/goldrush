#include "opt.hpp"
#include "read_hashing.hpp"
#include "spaced_seeds.hpp"

#include "Common/Options.h"
#include "MIBFConstructSupport.hpp"
#include "MIBloomFilter.hpp"
#include "multiLensfrHashIterator.hpp"

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
#include <unordered_set>
#include <vector>

bool
sort_by_sec(const pair<size_t, size_t>& a, const pair<size_t, size_t>& b)
{
  return (a.second > b.second);
}

void fill_bit_vector(const std::string& input_file, MIBFConstructSupport<uint32_t, multiLensfrHashIterator> &miBFCS, const size_t min_seq_len, const std::vector<std::string>& spaced_seeds) {
  std::cerr << "inserting bit vector" << std::endl;

  auto sTime = omp_get_wtime();

  btllib::SeqReader reader(input_file, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
  for (const auto record : reader) {
    if (record.seq.size() < min_seq_len) {
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
                        std::vector<bool>& tiles_assigned_bool_vec)
{

  size_t num_assigned_tiles = 0;
  size_t num_tiles = hashed_values.size();

  static std::unique_ptr<uint32_t[]> tiles_assigned_id_array;
  if (verbose) {
    tiles_assigned_id_array = std::make_unique<uint32_t[]>(num_tiles);
    std::memset(tiles_assigned_id_array.get(), 0, num_tiles * sizeof(uint32_t));
  }

#if _OPENMP
#pragma omp parallel for reduction(+ : num_assigned_tiles)
#endif
  for (size_t i = 1; i < num_tiles - 1;
       ++i) { // for each tile except first and last tiles. We consider them
              // erroneous and not used in tile assignment.
    std::unordered_map<uint32_t, uint32_t> id_counts;

    // Reusable vector for ranks
    vector<uint64_t> m_rank_pos(miBF.getHashNum());
    // Reusable vector for IDs
    vector<uint32_t> m_data(miBF.getHashNum());

    const auto& hashed_values_flat_array = hashed_values[i];
    for (size_t curr_frame = 0;
         curr_frame < (hashed_values_flat_array.size() / opt::hash_num);
         ++curr_frame) {

      std::unordered_set<uint32_t> unique_ids;
      if (miBF.atRank(hashed_values_flat_array.data() +
                        curr_frame * opt::hash_num,
                      m_rank_pos)) {                       // if its a hit
        m_data = miBF.getData(m_rank_pos);                 // m_data has ID's
        for (unsigned m = 0; m < miBF.getHashNum(); m++) { // iterate over ID's
          if (m_data[m] > miBF.s_mask) {                   // if ID is saturated
            const uint32_t new_id = m_data[m] & miBF.s_antiMask;
            if (new_id == 0) {
              continue;
            }
            unique_ids.insert(new_id);
          } else {
            if (m_data[m] == 0) {
              continue;
            }
            unique_ids.insert(m_data[m]);
          }
        }
      }

      for (const auto& unique_id :
           unique_ids) { // tabulate all unique ids to count table
        auto id_counts_it = id_counts.find(unique_id);
        if (id_counts_it != id_counts.end()) {
          ++(id_counts_it->second);
        } else {
          id_counts[unique_id] = 1;
        }
      }
    }

    uint32_t curr_id = 0;
    uint32_t curr_id_count = 0;

    for (auto id_counts_it = id_counts.begin(); id_counts_it != id_counts.end();
         ++id_counts_it) { // find id with highest count in a tile
      if (id_counts_it->second > curr_id_count) {
        curr_id = id_counts_it->first;
        curr_id_count = id_counts_it->second;
      }
    }

    if (curr_id_count > opt::threshold) {
      num_assigned_tiles++;
      tiles_assigned_bool_vec[i] = true;
    }

    if (verbose) {
      tiles_assigned_id_array[i] = curr_id;
    }
  }

  // print which id each tile is assigned to
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

  for (size_t i = 1; i < num_tiles - 1; ++i) {
    if (tiles_assigned_bool_vec[i] == false && tiles_assigned_id_vec[i] > 100) {
      if ((tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i - 1] &&
           tiles_assigned_bool_vec[i - 1] == true) ||
          (tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i + 1] &&
           tiles_assigned_bool_vec[i + 1] == true)) {
        tiles_assigned_bool_vec[i] = true;
      } else if ((tiles_assigned_id_vec[i] ==
                    tiles_assigned_id_vec[i - 1] + 1 &&
                  tiles_assigned_bool_vec[i - 1] == true) ||
                 (tiles_assigned_id_vec[i] ==
                    tiles_assigned_id_vec[i + 1] + 1 &&
                  tiles_assigned_bool_vec[i + 1] == true)) {
        tiles_assigned_bool_vec[i] = true;
      } else if ((tiles_assigned_id_vec[i] ==
                    tiles_assigned_id_vec[i - 1] - 1 &&
                  tiles_assigned_bool_vec[i - 1] == true) ||
                 (tiles_assigned_id_vec[i] ==
                    tiles_assigned_id_vec[i + 1] - 1 &&
                  tiles_assigned_bool_vec[i + 1] == true)) {
        tiles_assigned_bool_vec[i] = true;
      } else if (tiles_assigned_id_vec[i - 1] == tiles_assigned_id_vec[i + 1] &&
                 tiles_assigned_bool_vec[i - 1] == true &&
                 tiles_assigned_bool_vec[i + 1] == true) {
        tiles_assigned_bool_vec[i] = true;
        tiles_assigned_id_vec[i] = tiles_assigned_id_vec[i - 1];
      }
    }
  }

  for (size_t i = num_tiles - 2; i >= 1; --i) {
    if (tiles_assigned_bool_vec[i] == false && tiles_assigned_id_vec[i] > 100) {
      if ((tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i - 1] &&
           tiles_assigned_bool_vec[i - 1] == true) ||
          (tiles_assigned_id_vec[i] == tiles_assigned_id_vec[i + 1] &&
           tiles_assigned_bool_vec[i + 1] == true)) {
        tiles_assigned_bool_vec[i] = true;
      } else if ((tiles_assigned_id_vec[i] ==
                    tiles_assigned_id_vec[i - 1] + 1 &&
                  tiles_assigned_bool_vec[i - 1] == true) ||
                 (tiles_assigned_id_vec[i] ==
                    tiles_assigned_id_vec[i + 1] + 1 &&
                  tiles_assigned_bool_vec[i + 1] == true)) {
        tiles_assigned_bool_vec[i] = true;
      } else if ((tiles_assigned_id_vec[i] ==
                    tiles_assigned_id_vec[i - 1] - 1 &&
                  tiles_assigned_bool_vec[i - 1] == true) ||
                 (tiles_assigned_id_vec[i] ==
                    tiles_assigned_id_vec[i + 1] - 1 &&
                  tiles_assigned_bool_vec[i + 1] == true)) {
        tiles_assigned_bool_vec[i] = true;
      } else if (tiles_assigned_id_vec[i - 1] == tiles_assigned_id_vec[i + 1] &&
                 tiles_assigned_bool_vec[i - 1] == true &&
                 tiles_assigned_bool_vec[i + 1] == true) {
        tiles_assigned_bool_vec[i] = true;
        tiles_assigned_id_vec[i] = tiles_assigned_id_vec[i - 1];
      }
    }
  }

  // print which id each tile is assigned to
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
  // size_t curr_stretch = 0;
  std::vector<std::pair<size_t, size_t>> coord_vec;
  for (size_t i = 1; i < num_tiles - 1; ++i) {
    if (tiles_assigned_bool_vec[i] == false &&
        tiles_assigned_bool_vec[i - 1] == true) {
      start_idx = i;
    } else if (tiles_assigned_bool_vec[i] == true &&
               tiles_assigned_bool_vec[i - 1] == false) {
      end_idx = i - 1;
      coord_vec.push_back(std::make_pair(start_idx, end_idx));
    }
  }
  for (const auto& coords : coord_vec) {
    const auto& left = tiles_assigned_id_vec[coords.first - 1];
    const auto& right = tiles_assigned_id_vec[coords.second + 1];
    if (left == right || left == right + 1 || left == right - 1) {
      for (auto i = coords.first; i <= coords.second; ++i) {
        tiles_assigned_bool_vec[i] = true;
        tiles_assigned_id_vec[i] = left;
      }
    }
  }

  // print which id each tile is assigned to
  for (const auto& tiles_assigned_id : tiles_assigned_id_vec) {
    std::cerr << tiles_assigned_id << "\t";
  }
  std::cerr << std::endl;

  for (const auto& tiles_assigned_bool : tiles_assigned_bool_vec) {
    std::cerr << tiles_assigned_bool << "\t";
  }
  std::cerr << std::endl;

  for (size_t i = 2; i < num_tiles - 2; ++i) {
    if (tiles_assigned_bool_vec[i] == true && tiles_assigned_id_vec[i] > 100) {
      if (tiles_assigned_bool_vec[i - 1] == false &&
          tiles_assigned_bool_vec[i + 1] == false) {
        tiles_assigned_bool_vec[i] = false;
      }
    }
  }

  for (size_t i = num_tiles - 3; i >= 2; --i) {
    if (tiles_assigned_bool_vec[i] == true && tiles_assigned_id_vec[i] > 100) {
      if (tiles_assigned_bool_vec[i - 1] == false &&
          tiles_assigned_bool_vec[i + 1] == false) {
        tiles_assigned_bool_vec[i] = false;
      }
    }
  }

  // print which id each tile is assigned to
  for (const auto& tiles_assigned_id : tiles_assigned_id_vec) {
    std::cerr << tiles_assigned_id << "\t";
  }
  std::cerr << std::endl;

  for (const auto& tiles_assigned_bool : tiles_assigned_bool_vec) {
    std::cerr << tiles_assigned_bool << "\t";
  }
  std::cerr << std::endl;

  for (size_t i = 1; i < num_tiles - 1; ++i) {
    /*if (i == 0) {
        if (tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i + 1] &&
    tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i + 1]  - 1 &&
    tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i + 1]  + 1) {
            tiles_assigned_bool_vec[i] = false;
        }
        continue;
    } else if (i == num_tiles - 1) {
        if (tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i - 1] &&
    tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i - 1]  - 1 &&
    tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i - 1]  + 1) {
            tiles_assigned_bool_vec[i] = false;
        }
        continue;
    }*/
    if (tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i + 1] &&
        tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i + 1] - 1 &&
        tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i + 1] + 1 &&
        tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i - 1] &&
        tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i - 1] - 1 &&
        tiles_assigned_id_vec[i] != tiles_assigned_id_vec[i - 1] + 1) {
      tiles_assigned_bool_vec[i] = false;
    }
  }

  // print which id each tile is assigned to
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
  return num_assigned_tiles;
}

void process_read(const btllib::SeqReader::Record& record, const std::vector<std::vector<uint64_t>>& hashed_values, std::vector<std::ofstream>& golden_path_vec, std::vector<std::unique_ptr<MIBloomFilter<uint32_t>>>& mibf_vec, MIBFConstructSupport<uint32_t, multiLensfrHashIterator>& miBFCS, uint64_t& inserted_bases, uint64_t& target_bases, uint64_t& curr_path, uint32_t& id, uint32_t& ids_inserted, const size_t min_seq_len) {
      if (record.seq.size() < min_seq_len) {
        if (verbose) {
          std::cerr << "too short" << std::endl;
          std::cerr << "skipping: " << record.id << std::endl;
        }
        // wood_path <<  record.id << '\n' << record.seq <<  std::endl; skipping
        // wood path output to reduce time
        ++id;
        return;
      }
      if (id % 10000 == 0) {
        std::cerr << "processed " << id << " reads" << std::endl;
      }
      size_t len = record.seq.size();
      size_t num_tiles = len / opt::tile_length;

      if (verbose) {
        std::cerr << "name: " << record.id << std::endl;
        std::cerr << "num tiles: " << num_tiles - 2 << std::endl;
      }

      bool assigned = true;

      for (size_t level = 0; level < opt::levels; ++level) {
        auto& miBF = mibf_vec[level];
        if (verbose) {
          std::cerr << "current level : " << level << std::endl;
        }

        std::vector<uint32_t> tiles_assigned_id_vec(num_tiles, 0);
        std::vector<bool> tiles_assigned_bool_vec(num_tiles, false);
        const size_t num_assigned_tiles = calc_num_assigned_tiles(
          *miBF, hashed_values, tiles_assigned_id_vec, tiles_assigned_bool_vec);
        if (verbose) {
          std::cerr << "num assigned tiles: " << num_assigned_tiles
                    << std::endl;
        }
        const size_t num_unassigned_tiles = num_tiles - 2 - num_assigned_tiles;
        if (verbose) {
          std::cerr << "num unassigned tiles: " << num_unassigned_tiles
                    << std::endl;
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
          std::cerr << "unassigned" << std::endl;
          ++ids_inserted;

#if _OPENMP
#pragma omp parallel for
#endif
          for (size_t i = 0; i < num_tiles; ++i) {
            uint32_t curr_ids_inserted =
              ids_inserted + uint32_t((i + 1) * opt::tile_length / 10000);
            const auto& hashed_values_flat_array = hashed_values[i];
            miBFCS.insertMIBF(*miBF,
                              hashed_values_flat_array,
                              curr_ids_inserted); //, non_singletons_bf_vec);
            // miBFCS.insertSaturation(*miBF, Hhashes, ids_inserted); // don't
            // care about saturation atm so skipping for speed
            // }
          }

          ids_inserted = ids_inserted + uint32_t(record.seq.size() / 10000);
          // output read to golden path

          golden_path_vec[level] << ">" << record.id << '\n'
                                 << record.seq << std::endl;

          inserted_bases += record.seq.size();
          if (opt::temp_mode || opt::new_temp_mode) {
            if (target_bases < inserted_bases) {
              ++curr_path;
              if (opt::max_paths < curr_path) {
                exit(0);
              }
              inserted_bases = 0;
              mibf_vec.pop_back();
              mibf_vec.emplace_back(std::unique_ptr<MIBloomFilter<uint32_t>>(
                miBFCS.getEmptyMIBF()));
              golden_path_vec.pop_back();
              golden_path_vec.emplace_back(
                std::ofstream(opt::prefix_file + "_golden_path_" +
                              std::to_string(curr_path) + ".fa"));
              ids_inserted = 0;
            }
          }
          break; // breaks the level loop

        } else {
          if (opt::temp_mode) {
            continue;
          }
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
          std::cerr << "checkpoint 1" << std::endl;
          for (size_t i = 1; i < num_tiles - 1; ++i) {
            if (tiles_assigned_bool_vec[i] == false &&
                tiles_assigned_bool_vec[i - 1] == true) {
              start_idx = i;
              curr_stretch = 1;
            } else if ((tiles_assigned_bool_vec[i] == false &&
                        tiles_assigned_bool_vec[i] ==
                          tiles_assigned_bool_vec[i - 1]) &&
                       (i + 1 != num_tiles - 1)) {
              ++curr_stretch;
            } else if (tiles_assigned_bool_vec[i] == true &&
                       tiles_assigned_bool_vec[i] !=
                         tiles_assigned_bool_vec[i - 1]) {
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
          std::cerr << "trim_start_idx: " << trim_start_idx << std::endl;
          std::cerr << "trim_end_idx: " << trim_end_idx << std::endl;

          if (num_tiles < 15) {
            bool good_right_flank = false;
            bool good_left_flank = false;
            std::unordered_map<size_t, size_t> left_flank;
            for (ssize_t i = longest_start_idx - 1; i >= 0; --i) {
              if (left_flank.find(tiles_assigned_id_vec[i]) !=
                  left_flank.end()) {
                ++left_flank[tiles_assigned_id_vec[i]];
              } else {
                left_flank[tiles_assigned_id_vec[i]] = 1;
              }
            }
            std::cerr << "checkpoint 2.1" << std::endl;
            std::vector<std::pair<size_t, size_t>> left_flank_vec;
            for (const auto& myPair : left_flank) {
              left_flank_vec.push_back(
                std::make_pair(myPair.first, myPair.second));
            }
            std::cerr << "checkpoint 2.2" << std::endl;
            sort(left_flank_vec.begin(), left_flank_vec.end(), sort_by_sec);
            std::cerr << "checkpoint 2.3" << std::endl;
            std::cerr << left_flank_vec.size() << std::endl;
            if (left_flank_vec.size() != 0) {
              if (left_flank_vec[0].second >= 2) {
                trim_start_idx = longest_start_idx - 2;
                good_left_flank = true;
              } else if (left_flank_vec.size() >= 2 &&
                         (left_flank_vec[0].second + left_flank_vec[1].second >
                            3 &&
                          (left_flank_vec[0].first - left_flank_vec[1].first ==
                             1 ||
                           left_flank_vec[1].first - left_flank_vec[0].first ==
                             1))) {
                trim_start_idx = longest_start_idx - 2;
                good_left_flank = true;

              } else if (trim_start_idx == 0) {
                good_left_flank = true;
              }
            }
            if (good_left_flank) {
              std::cerr << "good left flank: true" << std::endl;
            } else {
              std::cerr << "good left flank: false" << std::endl;
            }
            std::unordered_map<size_t, size_t> right_flank;
            for (ssize_t i = longest_end_idx + 1; i < (ssize_t)num_tiles; ++i) {
              if (right_flank.find(tiles_assigned_id_vec[i]) !=
                  right_flank.end()) {
                ++right_flank[tiles_assigned_id_vec[i]];
              } else {
                right_flank[tiles_assigned_id_vec[i]] = 1;
              }
            }
            std::vector<std::pair<size_t, size_t>> right_flank_vec;
            for (const auto& myPair : right_flank) {
              right_flank_vec.push_back(
                std::make_pair(myPair.first, myPair.second));
            }
            sort(right_flank_vec.begin(), right_flank_vec.end(), sort_by_sec);
            if (right_flank_vec.size() != 0) {
              if (right_flank_vec[0].second >= 2) {
                trim_end_idx = longest_end_idx + 2;
                good_right_flank = true;
              } else if (right_flank_vec.size() >= 2 &&
                         (right_flank_vec[0].second +
                              right_flank_vec[1].second >
                            3 &&
                          (right_flank_vec[0].first -
                               right_flank_vec[1].first ==
                             1 ||
                           right_flank_vec[1].first -
                               right_flank_vec[0].first ==
                             1))) {
                trim_end_idx = longest_end_idx + 2;
                good_right_flank = true;

              } else if (trim_end_idx == num_tiles - 1) {
                good_right_flank = true;
              }
            }
            if (good_right_flank) {
              std::cerr << "good right flank: true" << std::endl;
            } else {
              std::cerr << "good right flank: false" << std::endl;
            }
            if (good_left_flank && good_right_flank) {
              good_flank = true;
            }

          } else {
            std::cerr << "checkpoint 2" << std::endl;
            // bool valid = true;

            if (longest_start_idx - 5 >= 1) {
              std::unordered_map<size_t, size_t> left_flank;
              for (ssize_t i = longest_start_idx - 5; i < longest_start_idx;
                   ++i) {
                if (left_flank.find(tiles_assigned_id_vec[i]) !=
                    left_flank.end()) {
                  ++left_flank[tiles_assigned_id_vec[i]];
                } else {
                  left_flank[tiles_assigned_id_vec[i]] = 1;
                }
              }
              std::cerr << "checkpoint 2.1" << std::endl;
              std::vector<std::pair<size_t, size_t>> left_flank_vec;
              for (const auto& myPair : left_flank) {
                left_flank_vec.push_back(
                  std::make_pair(myPair.first, myPair.second));
              }
              std::cerr << "checkpoint 2.2" << std::endl;
              sort(left_flank_vec.begin(), left_flank_vec.end(), sort_by_sec);
              std::cerr << "checkpoint 2.3" << std::endl;
              std::cerr << left_flank_vec.size() << std::endl;
              if (left_flank_vec[0].second >= 2) {
                trim_start_idx = longest_start_idx - 2;
                good_flank = true;
              } else if (left_flank_vec[0].second + left_flank_vec[1].second >
                           3 &&
                         (left_flank_vec[0].first - left_flank_vec[1].first ==
                            1 ||
                          left_flank_vec[1].first - left_flank_vec[0].first ==
                            1)) {
                trim_start_idx = longest_start_idx - 2;
                good_flank = true;
              }
              std::cerr << "checkpoint 2.4" << std::endl;
            } else {
              trim_start_idx = 0;
            }
            std::cerr << "checkpoint 3" << std::endl;

            if (longest_end_idx + 5 < (ssize_t)num_tiles - 1) {
              std::unordered_map<size_t, size_t> right_flank;
              for (ssize_t i = longest_end_idx + 5; i > longest_end_idx; --i) {
                if (right_flank.find(tiles_assigned_id_vec[i]) !=
                    right_flank.end()) {
                  ++right_flank[tiles_assigned_id_vec[i]];
                } else {
                  right_flank[tiles_assigned_id_vec[i]] = 1;
                }
              }
              std::vector<std::pair<size_t, size_t>> right_flank_vec;
              for (const auto& myPair : right_flank) {
                right_flank_vec.push_back(
                  std::make_pair(myPair.first, myPair.second));
              }
              sort(right_flank_vec.begin(), right_flank_vec.end(), sort_by_sec);
              if (right_flank_vec[0].second >= 2) {
                trim_end_idx = longest_end_idx + 2;
                good_flank = true;
              } else if (right_flank_vec[0].second + right_flank_vec[1].second >
                           3 &&
                         (right_flank_vec[0].first - right_flank_vec[1].first ==
                            1 ||
                          right_flank_vec[1].first - right_flank_vec[0].first ==
                            1)) {
                trim_end_idx = longest_end_idx + 2;
                good_flank = true;
              }

            } else {
              trim_end_idx = (ssize_t)num_tiles - 1;
            }
          }
          std::cerr << "checkpoint 4" << std::endl;
          if (good_flank) {
            assigned = false;
            std::cerr << "trimmed" << std::endl;
            ++ids_inserted;

#if _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = trim_start_idx; i <= trim_end_idx; ++i) {
              uint32_t curr_ids_inserted =
                ids_inserted +
                uint32_t((i - trim_start_idx + 1) * opt::tile_length / 10000);
              const auto& hashed_values_flat_array = hashed_values[i];
              miBFCS.insertMIBF(*miBF,
                                hashed_values_flat_array,
                                curr_ids_inserted); //, non_singletons_bf_vec);
              // miBFCS.insertSaturation(*miBF, Hhashes, ids_inserted); // don't
              // care about saturation atm so skipping for speed
              // }
            }
            ids_inserted =
              ids_inserted +
              uint32_t((trim_end_idx - trim_start_idx) * 1000 / 10000);
            // output read to golden path
            if (trim_end_idx == num_tiles - 1) {
              std::string new_seq =
                record.seq.substr(trim_start_idx * 1000, std::string::npos);
              golden_path_vec[level] << ">trimmed" << record.id << '\n'
                                     << new_seq << std::endl;
              inserted_bases += new_seq.size();
            } else {
              std::string new_seq =
                record.seq.substr(trim_start_idx * 1000,
                                  (trim_end_idx - trim_start_idx + 1) * 1000);
              golden_path_vec[level] << ">trimmed" << record.id << '\n'
                                     << new_seq << std::endl;
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
                mibf_vec.emplace_back(std::unique_ptr<MIBloomFilter<uint32_t>>(
                  miBFCS.getEmptyMIBF()));
                golden_path_vec.pop_back();
                golden_path_vec.emplace_back(
                  std::ofstream(opt::prefix_file + "_golden_path_" +
                                std::to_string(curr_path) + ".fa"));
                ids_inserted = 0;
              }
            }

            break; // breaks the level loop
          }
        }
      }
      if (assigned) {
        std::cerr << "assigned" << std::endl;
        // output read to wood path
        // wood_path <<  record.id << '\n' << record.seq <<  std::endl; skipping
        // wood path output to reduce time

        ++id;
      }
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
            << "seed pattens: " << opt::hash_num << "\n"
            << "threshold: " << opt::threshold << "\n"
            << "base seed pattern: " << seed_string_vec[0] << "\n"
            << "minimum unassigned tiles: " << opt::unassigned_min << "\n"
            << "maximum assigned tiles: " << opt::assigned_max << "\n"
            << "genome size: " << opt::genome_size << "\n"
            << "occupancy: " << opt::occupancy << "\n"
            << "jobs: " << opt::jobs << std::endl;

  std::vector<std::ofstream> golden_path_vec;
  for (size_t level = 0; level < opt::levels; ++level) {
    golden_path_vec.emplace_back(std::ofstream(
      opt::prefix_file + "_golden_path_" + std::to_string(level) + ".fa"));
  }

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

  fill_bit_vector(opt::input, miBFCS, opt::min_length, seed_string_vec);

  // setting up MIBF
  miBFCS.setup();
  std::vector<std::unique_ptr<MIBloomFilter<uint32_t>>> mibf_vec;
  for (size_t i = 0; i < opt::levels; ++i) {
    mibf_vec.emplace_back(
      std::unique_ptr<MIBloomFilter<uint32_t>>(miBFCS.getEmptyMIBF()));
  }

  std::cerr << "assigning tiles" << std::endl;
  sTime = omp_get_wtime();

  btllib::OrderQueueMPMC<ReadTileHashes> precomputed_hash_queue(
    btllib::SeqReader::LONG_MODE_BUFFER_SIZE,
    btllib::SeqReader::LONG_MODE_BLOCK_SIZE);
  start_read_hashing(opt::input,
                     opt::min_length,
                     opt::tile_length,
                     opt::kmer_size,
                     seed_string_vec,
                     precomputed_hash_queue,
                     6);

  uint64_t inserted_bases = 0;
  uint64_t target_bases = opt::ratio * opt::target_size;
  uint64_t curr_path = 1;
  uint32_t id = 1;
  uint32_t ids_inserted = 0;
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

      process_read(record, hashed_values, golden_path_vec, mibf_vec, miBFCS, inserted_bases, target_bases, curr_path, id, ids_inserted, opt::min_length);
    }
  }

  std::cerr << "assigned" << std::endl;
  std::cerr << "in " << setprecision(4) << fixed << omp_get_wtime() - sTime
            << "\n";
}
