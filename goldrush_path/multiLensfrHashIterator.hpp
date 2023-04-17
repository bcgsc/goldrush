/*
 * HashFunction.hpp
 *
 * 	Created to get around templating issues revolving ntHashIterator and
 * sfrHashIterator
 *
 *
 *  Created on: Aug. 25, 2020
 *      Author: cjustin
 */

#ifndef MULTI_LENGTH_STHASHITERATOR_HPP_
#define MULTI_LENGTH_STHASHITERATOR_HPP_
#include <btllib/nthash.hpp>
#include <assert.h>
#include <vector>

class multiLensfrHashIterator
{
public:
  /*
   * Default constructor.
   */
  multiLensfrHashIterator()
    : m_hVec(NULL)
    , m_pos(std::numeric_limits<std::size_t>::max())
  {}

  multiLensfrHashIterator(const std::string& seq,
                          const std::vector<std::string>& seedVec)
    : m_itrVec(std::vector<btllib::SeedNtHash>())
    , m_seed(seedVec.size())
    , m_hash(m_seed)
    , m_hVec(new uint64_t[m_hash])
    , m_pos(0)
  {
    for (size_t i = 0; i < m_seed; ++i) {
      std::vector<std::string> seed{ seedVec[i] };
      m_itrVec.push_back(btllib::SeedNtHash(seq, seed, 1, seedVec[i].size()));
      m_itrVec[i].roll();
      m_hVec[i] = m_itrVec[i].hashes()[0];
    }
  }

  /** get pointer to hash values for current k-mer */
  const uint64_t* operator*() const { return m_hVec; }

  /** pre-increment operator */
  multiLensfrHashIterator& operator++()
  {
    bool update = false;
    ++m_pos;
    for (size_t i = 0; i < m_seed; ++i) {
      if (m_itrVec[i].roll()) {
        update = true;
      } else {
        continue;
      }

      m_hVec[i] = m_itrVec[i].hashes()[0];
    }

    if (!update) {
      m_pos = std::numeric_limits<std::size_t>::max();
    }

    return *this;
  }
  /** iterator pointing to one past last element */
  static const multiLensfrHashIterator end()
  {
    return multiLensfrHashIterator();
  }

  /** destructor */
  ~multiLensfrHashIterator()
  {
    if (m_hVec != NULL) {
      delete[] m_hVec;
    }
  }
  /** test equality with another iterator */
  bool operator==(const multiLensfrHashIterator& it) const
  {
    return m_pos == it.m_pos;
  }

  /** test inequality with another iterator */
  bool operator!=(const multiLensfrHashIterator& it) const
  {
    return !(*this == it);
  }

private:
  std::vector<btllib::SeedNtHash> m_itrVec;
  unsigned m_seed;
  unsigned m_hash;
  uint64_t* m_hVec;
  size_t m_pos;
};

#endif /* MULTI_LENGTH_STHASHITERATOR_HPP_ */
