#include "multiLensfrHashIterator.hpp"

#include "read_hashing.hpp"

#include <future>

size_t
read_hashing(btllib::SeqReader& reader,
             const size_t tile_size,
             const size_t min_seq_len,
             const size_t k,
             const std::vector<std::string>& spaced_seeds,
             btllib::OrderQueueMPMC<ReadTileHashes>& read_tile_hashes_queue)
{
  btllib::OrderQueueMPMC<btllib::SeqReader::Record>::Block record_block(
    reader.get_block_size());
  btllib::OrderQueueMPMC<ReadTileHashes>::Block block(reader.get_block_size());

  size_t last_block_num = 0;
  while (true) {
    record_block = reader.read_block();
    if (record_block.count == 0) {
      break;
    }

    for (size_t idx = 0; idx < record_block.count; idx++) {
      auto record = record_block.data[idx];
      const size_t len = record.seq.size();
      const size_t num_tiles = len / tile_size;

      std::vector<std::vector<uint64_t>> hashed_values(num_tiles,
                                                       std::vector<uint64_t>());
      if (record.seq.size() >= min_seq_len) {
        for (size_t i = 0; i < num_tiles; ++i) {
          std::string tile_seq =
            record.seq.substr(i * tile_size, tile_size + k - 1);
          multiLensfrHashIterator itr(tile_seq, spaced_seeds);
          while (itr != itr.end()) {
            for (size_t curr_hash = 0; curr_hash < spaced_seeds.size();
                 ++curr_hash) {
              hashed_values[i].push_back((*itr)[curr_hash]);
            }
            ++itr;
          }
        }
      }

      block.data[block.count++] = { std::move(record),
                                    std::move(hashed_values) };
      if (block.count == block.data.size()) {
        block.num = record_block.num;
        last_block_num = block.num;
        read_tile_hashes_queue.write(block);
        block.count = 0;
        block.num = 0;
      }
    }
  }
  if (block.count > 0) {
    block.num = last_block_num;
    read_tile_hashes_queue.write(block);
    block.count = 0;
    block.num = 0;
  }
  return last_block_num;
}

void
start_read_hashing(
  const std::string& reads_filepath,
  const size_t tile_size,
  const size_t min_seq_len,
  const size_t k,
  const std::vector<std::string>& spaced_seeds,
  btllib::OrderQueueMPMC<ReadTileHashes>& read_tile_hashes_queue,
  const unsigned worker_num)
{
  (new std::thread([&]() {
    btllib::SeqReader reader(reads_filepath,
                             btllib::SeqReader::Flag::LONG_MODE);
    std::vector<std::future<size_t>> last_block_nums_futures;
    std::vector<size_t> last_block_nums;

    for (unsigned i = 0; i < worker_num; i++) {
      last_block_nums_futures.push_back(
        std::async(read_hashing,
                   std::ref(reader),
                   tile_size,
                   min_seq_len,
                   k,
                   std::cref(spaced_seeds),
                   std::ref(read_tile_hashes_queue)));
    }
    for (auto& f : last_block_nums_futures) {
      last_block_nums.push_back(f.get());
    }

    btllib::OrderQueueMPMC<ReadTileHashes>::Block dummy(
      reader.get_block_size());
    dummy.num =
      *(std::max_element(last_block_nums.begin(), last_block_nums.end())) + 1;
    dummy.count = 0;
    read_tile_hashes_queue.write(dummy);
  }))
    ->detach();
}