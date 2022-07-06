/*
 * MIBFConstructSupport.hpp
 *
 * Purpose: To provide support for easier filter construction
 *
 * IDs are still managed by the class calling this object
 *
 * In order to use object one must have a iterator with same interface
 * as ntHashIterator
 *
 * Assumes saturation bit is being used
 * TODO: add functionality of strand bit
 *
 *  Created on: Mar 2028
 *      Author: Justin Chu
 */

#ifndef MIBFCONSTRUCTSUPPORT_HPP_
#define MIBFCONSTRUCTSUPPORT_HPP_

#include "MIBloomFilter.hpp"
#include <google/dense_hash_map>
#include <google/dense_hash_set>
#include <google/sparse_hash_map>
#include <iostream>
#include <sdsl/bit_vector_il.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <tuple>
#include <vector>

#if _OPENMP
#include <omp.h>
#endif

// T = T ID type, H = rolling hash itr
template<typename T, class H>
class MIBFConstructSupport
{
public:
  static const unsigned BLOCKSIZE = 512;
  /*
   * numhashfunctions may also mean number of spaced seeds
   */
  MIBFConstructSupport(size_t expectedEntries,
                       unsigned k,
                       unsigned numHashFunction,
                       double occupancy,
                       const vector<string>& spacedSeeds = vector<string>(0))
    : m_isBVMade(false)
    , m_isMIBFMade(false)
    , m_expectedEntries(expectedEntries)
    , m_k(k)
    , m_h(numHashFunction)
    , m_occupancy(occupancy)
    , m_spacedSeeds(spacedSeeds)
    , m_counts(vector<T>())
  {
    m_filterSize =
      MIBloomFilter<T>::calcOptimalSize(m_expectedEntries, m_h, m_occupancy);
    std::cerr << "m_filterSize: " << m_filterSize << std::endl;
    m_bv = sdsl::bit_vector(m_filterSize);
  }

  MIBFConstructSupport(size_t expectedEntries,
                       unsigned k,
                       unsigned numHashFunction,
                       double occupancy,
                       size_t filterSize,
                       const vector<string>& spacedSeeds = vector<string>(0))
    : m_isBVMade(false)
    , m_isMIBFMade(false)
    , m_expectedEntries(expectedEntries)
    , m_k(k)
    , m_h(numHashFunction)
    , m_occupancy(occupancy)
    , m_spacedSeeds(spacedSeeds)
    , m_counts(vector<T>())
    , m_filterSize(filterSize)
  {
    std::cerr << "m_filterSize: " << m_filterSize << std::endl;
    m_bv = sdsl::bit_vector(m_filterSize);
  }

  ~MIBFConstructSupport()
  {
    // assert(m_isBVMade & m_isMIBFMade);
  }

  /*
   * Returns count of collisions (counts unique k-mers)
   */
  inline size_t insertBVColli(H& itr)
  {
    // assert(!m_isBVMade);
    size_t count = 0;
    /* init rolling hash state and compute hash values for first k-mer */
    for (; itr != itr.end(); ++itr) {
      unsigned colliCount = 0;
      for (unsigned i = 0; i < m_h; ++i) {
        uint64_t pos = (*itr)[i] % m_bv.size();
        uint64_t* dataIndex = m_bv.data() + (pos >> 6);
        uint64_t bitMaskValue = (uint64_t)1 << (pos & 0x3F);
        colliCount +=
          __sync_fetch_and_or(dataIndex, bitMaskValue) >> (pos & 0x3F) & 1;
      }
      if (colliCount == m_h) {
        ++count;
      }
    }
    return count;
  }



  void insertBV(H& itr)
  {

    /* init rolling hash state and compute hash values for first k-mer */
    for (; itr != itr.end(); ++itr) {
      for (unsigned i = 0; i < m_h; ++i) {
        uint64_t pos = (*itr)[i] % m_bv.size();
        uint64_t* dataIndex = m_bv.data() + (pos >> 6);
        uint64_t bitMaskValue = (uint64_t)1 << (pos & 0x3F);
        (void)(__sync_fetch_and_or(dataIndex, bitMaskValue) >> (pos & 0x3F) &
               1);
      }
    }
  }

  void insertBV(const std::vector<uint64_t>& hash_vec)
  {

    /* init rolling hash state and compute hash values for first k-mer */
    for (const auto& hash : hash_vec) {
      uint64_t pos = hash % m_bv.size();
      uint64_t* dataIndex = m_bv.data() + (pos >> 6);
      uint64_t bitMaskValue = (uint64_t)1 << (pos & 0x3F);
      (void)(__sync_fetch_and_or(dataIndex, bitMaskValue) >> (pos & 0x3F) & 1);
    }
  }

  /*
   * Set up MIBF data structures for reuse
   */

  void setup()
  {
    m_bv_il = sdsl::bit_vector_il<BLOCKSIZE>(m_bv);
    m_bv = sdsl::bit_vector();
    m_rankSupport = sdsl::rank_support_il<1>(&m_bv_il);
  }

  /*
   * Generate empty miBF, can currently only be called once per object
   */
  MIBloomFilter<T>* getEmptyMIBF()
  {
    MIBloomFilter<T>* miBF =
      new MIBloomFilter<T>(m_h, m_k, m_bv_il, m_rankSupport, m_spacedSeeds);
    m_counts = vector<T>(miBF->getPop(), 0);
    return miBF;
  }

  void reset_counts()
  {
    memset(&m_counts[0], 0, m_counts.size() * sizeof m_counts[0]);
  }

  /*
   * Uses single value Reservoir sampling
   * pair<ID,ID> first ID stores the currentID, and ID stores the current
   * observation count If the second ID exceeds max possible count, the ID is a
   * critical ID Critical IDs are needed for partial hits and always replace
   * existing IDs
   *
   * Once saturation is set, insertions are prevented
   */
  void insertMIBF(MIBloomFilter<T>& miBF, H& itr, T id)
  {
    // assert(m_isBVMade & !m_isMIBFMade);
    // get positions
    hashSet values;
    values.set_empty_key(miBF.size());
    while (itr != itr.end()) {
      for (unsigned i = 0; i < m_h; ++i) {
        values.insert((*itr)[i]);
      }
      ++itr;
    }
    for (hashSet::iterator itr = values.begin(); itr != values.end(); itr++) {
      uint64_t randomSeed = *itr ^ id;
      uint64_t rank = miBF.getRankPos(*itr);
      T count = __sync_add_and_fetch(&m_counts[rank], 1);
      T randomNum = std::hash<T>{}(randomSeed) % count;
      if (randomNum == count - 1) {
        miBF.setData(rank, id);
      }
    }
  }

  void insertMIBF(MIBloomFilter<T>& miBF,
                  const std::vector<uint64_t>& hash_vec,
                  T id)
  {
    // assert(m_isBVMade & !m_isMIBFMade);
    // get positions
#if _OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < hash_vec.size(); ++i) {
      const auto& hash = hash_vec[i];
      uint64_t randomSeed = hash ^ id;
      uint64_t rank = miBF.getRankPos(hash);
      T count = __sync_add_and_fetch(&m_counts[rank], 1);
      T randomNum = std::hash<T>{}(randomSeed) % count;
      // std::cerr << "id: " << id << " randomNum: " << randomNum << "
      // randomSeed: " << randomSeed <<std::endl;
      if (randomNum == count - 1) {
        miBF.setData(rank, id);
      }
    }
  }

  void insertMIBF(MIBloomFilter<T>& miBF,
                  const std::vector<std::vector<uint64_t>>& hash_vec,
                  size_t start,
                  size_t end,
                  T id)
  {
    // assert(m_isBVMade & !m_isMIBFMade);
    // get positions
    size_t vec_size = hash_vec[0].size();
    size_t num_elements = 0;
    for (size_t i = start; i < end; ++i) {
      num_elements += hash_vec[i].size();
    }

#if _OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < num_elements; ++i) {
      size_t vec_num = i / vec_size;
      size_t hash_loc = i % vec_size;
      const auto& hash = hash_vec[start + vec_num][hash_loc];
      uint64_t randomSeed = hash ^ id;
      uint64_t rank = miBF.getRankPos(hash);
      T count = __sync_add_and_fetch(&m_counts[rank], 1);
      T randomNum = std::hash<T>{}(randomSeed) % count;
      if (randomNum == count - 1) {
        miBF.setData(rank, id);
      }
    }
  }


  void insertSaturation(MIBloomFilter<T>& miBF, H& itr, T id)
  {

    typedef google::dense_hash_set<uint64_t> SatSet;
    SatSet satVal;
    satVal.set_empty_key(miBF.size());
    setSatIfMissing(miBF, id, itr);
  }

  /*
   * Returns number of bits in top level bit vector
   */
  size_t getFilterSize() const { return m_filterSize; }

private:
  typedef google::dense_hash_set<uint64_t> hashSet;

  bool m_isBVMade;
  bool m_isMIBFMade;
  size_t m_expectedEntries;
  unsigned m_k, m_h;
  double m_occupancy;
  const vector<string>& m_spacedSeeds;
  sdsl::bit_vector m_bv;
  sdsl::bit_vector_il<BLOCKSIZE> m_bv_il;
  sdsl::rank_support_il<1> m_rankSupport;
  vector<T> m_counts;
  size_t m_filterSize;

  /*
   * Attempts of mutate values to prevent saturation
   * If unable to save values it will saturate values
   * Small chance that mutation may erase entries
   */
  inline void setSatIfMissing(MIBloomFilter<T>& miBF, T id, H& itr)
  {
    while (itr != itr.end()) {
      // for each set of hash values, check for saturation
      vector<uint64_t> rankPos = miBF.getRankPos(*itr);
      vector<T> results = miBF.getData(rankPos);
      vector<T> replacementIDs(m_h);
      bool valueFound = false;
      vector<T> seenSet(m_h);
      for (unsigned i = 0; i < m_h; ++i) {
        T currentResult = results[i] & MIBloomFilter<T>::s_antiMask;
        if (currentResult == id) {
          valueFound = true;
          break;
        }
        if (find(seenSet.begin(), seenSet.end(), currentResult) ==
            seenSet.end()) {
          seenSet.push_back(currentResult);
        } else {
          replacementIDs.push_back(currentResult);
        }
      }
      if (!valueFound) {
        uint64_t replacementPos = m_counts.size();
        T minCount = numeric_limits<T>::min();
        for (unsigned i = 0; i < m_h; ++i) {
          T currentResult = results[i] & MIBloomFilter<T>::s_antiMask;
          if (find(replacementIDs.begin(),
                   replacementIDs.end(),
                   currentResult) != replacementIDs.end()) {
            if (minCount < m_counts[rankPos[i]]) {
              minCount = m_counts[rankPos[i]];
              replacementPos = rankPos[i];
            }
          }
        }
        // mutate if possible
        if (replacementPos != m_counts.size()) {
          miBF.setData(replacementPos, id);
#pragma omp atomic update
          ++m_counts[replacementPos];
        } else {
          miBF.saturate(*itr);
        }
      }
      ++itr;
    }
  }
};
#endif /* MIBFCONSTRUCTSUPPORT_HPP_ */
