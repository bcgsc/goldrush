#ifndef READ_HASHING_HPP
#define READ_HASHING_HPP

#include "btllib/order_queue.hpp"
#include "btllib/seq_reader.hpp"

#include <string>
#include <unordered_set>
#include <vector>

struct ReadTileHashes
{
  btllib::SeqReader::Record read;
  std::vector<std::vector<uint64_t>> tile_hashes;
};

void
start_read_hashing(
  const std::string& reads_filepath,
  size_t tile_size,
  size_t min_seq_len,
  size_t k,
  const std::vector<std::string>& spaced_seeds,
  btllib::OrderQueueMPMC<ReadTileHashes>& read_tile_hashes_queue,
  const unsigned worker_num,
  const std::unordered_set<std::string>& filter_out_reads);
#endif