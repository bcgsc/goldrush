/*
 *
 * ntcard.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef NTCARD_H_
#define NTCARD_H_

#include "btllib/seq_reader.hpp"
#include "multiLensfrHashIterator.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace nts {
unsigned nThrd = 1;
unsigned kmLen = 64;
size_t rBuck;
unsigned rBits = 27;
unsigned sBits = 11;
unsigned sMask;
unsigned covMax = 10000;
size_t nSamp = 2;
bool samH = true;
unsigned num_seeds = 0;
std::vector<std::string> seed_vec;
} // namespace nts

size_t
getInf(const char* inFile)
{
  std::ifstream in(inFile, std::ifstream::ate | std::ifstream::binary);
  return in.tellg();
}

unsigned
getftype(std::ifstream& in, std::string& samSeq)
{
  std::string hseq;
  getline(in, hseq);
  if (hseq[0] == '>') {
    return 1;
  }
  if (hseq[0] == '@') {
    if ((hseq[1] == 'H' && hseq[2] == 'D') ||
        (hseq[1] == 'S' && hseq[2] == 'Q') ||
        (hseq[1] == 'R' && hseq[2] == 'G') ||
        (hseq[1] == 'P' && hseq[2] == 'G') ||
        (hseq[1] == 'C' && hseq[2] == 'O')) {
      return 2;
    } else
      return 0;
  }
  std::istringstream alnSec(hseq);
  std::string s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11;
  alnSec >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9 >> s10 >> s11;
  if ((s2.find_first_not_of("0123456789") == std::string::npos) &&
      (s5.find_first_not_of("0123456789") == std::string::npos)) {
    nts::samH = false;
    samSeq = hseq;
    return 2;
  }
  return 3;
}

inline void
ntComp(const uint64_t hVal, std::vector<uint16_t>& t_Counter)
{
  uint64_t indBit = nts::nSamp;
  if (hVal >> (63 - nts::sBits) == 1)
    indBit = 0;
  if (hVal >> (64 - nts::sBits) == nts::sMask)
    indBit = 1;
  if (indBit < nts::nSamp) {
    size_t shVal = hVal & (nts::rBuck - 1);
#pragma omp atomic
    ++t_Counter[indBit * nts::rBuck + shVal];
  }
}

inline void
stRead(const string& seq,
       std::vector<std::vector<uint16_t>>& t_Counters,
       std::vector<size_t>& totKmer_vec)
{
  multiLensfrHashIterator itr(seq, nts::seed_vec);

  while (itr != itr.end()) {
    for (size_t i = 0; i < nts::num_seeds; ++i) {
      auto& t_Counter = t_Counters[i];
      auto& totKmer = totKmer_vec[i];
      ntComp((*itr)[i], t_Counter);
      ++totKmer;
    }
    ++itr;
  }
}

void
compEst(const std::vector<uint16_t>& t_Counter, double& F0Mean, double fMean[])
{
  unsigned** p = new unsigned*[nts::nSamp];
  for (size_t i = 0; i < nts::nSamp; i++) {
    p[i] = new unsigned[65536];
  }
  for (size_t i = 0; i < nts::nSamp; i++)
    for (size_t j = 0; j < 65536; j++)
      p[i][j] = 0;

  for (size_t i = 0; i < nts::nSamp; i++)
    for (size_t j = 0; j < nts::rBuck; j++)
      ++p[i][t_Counter[i * nts::rBuck + j]];

  double pMean[65536];
  for (size_t i = 0; i < 65536; i++)
    pMean[i] = 0.0;
  for (size_t i = 0; i < 65536; i++) {
    for (size_t j = 0; j < nts::nSamp; j++)
      pMean[i] += p[j][i];
    pMean[i] /= 1.0 * nts::nSamp;
  }

  F0Mean = (ssize_t)((nts::rBits * log(2) - log(pMean[0])) * 1.0 *
                     ((size_t)1 << (nts::sBits + nts::rBits)));
  for (size_t i = 0; i < 65536; i++)
    fMean[i] = 0;
  fMean[1] =
    -1.0 * pMean[1] / (pMean[0] * (log(pMean[0]) - nts::rBits * log(2)));
  for (size_t i = 2; i < 65536; i++) {
    double sum = 0.0;
    for (size_t j = 1; j < i; j++)
      sum += j * pMean[i - j] * fMean[j];
    fMean[i] =
      -1.0 * pMean[i] / (pMean[0] * (log(pMean[0]) - nts::rBits * log(2))) -
      sum / (i * pMean[0]);
  }
  for (size_t i = 1; i < 65536; i++)
    fMean[i] = abs((ssize_t)(fMean[i] * F0Mean));
}

std::vector<std::vector<size_t>>
getHist(const vector<string>& inFiles,
        const unsigned kLen,
        const unsigned nThr,
        const unsigned covMax,
        std::vector<std::string> seed_vec = std::vector<std::string>(0))
{
#ifdef _OPENMP
  double sTime = omp_get_wtime();
#endif

  nts::nThrd = nThr;
  nts::kmLen = kLen;
  nts::num_seeds = seed_vec.size();
  nts::seed_vec = seed_vec;
  nts::covMax = covMax;
  size_t num_histArray = 1;
  if (nts::num_seeds) {
    num_histArray = nts::num_seeds;
  }
  std::vector<std::vector<size_t>> histArray_vec(
    num_histArray, std::vector<size_t>(covMax + 2, 0));

  size_t totalSize = 0;
  for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i)
    totalSize += getInf(inFiles[file_i].c_str());
  if (totalSize < 50000000000)
    nts::sBits = 7;

  std::vector<size_t> totalKmers_vec(num_histArray, 0);

  nts::rBuck = ((size_t)1) << nts::rBits;
  nts::sMask = (((size_t)1) << (nts::sBits - 1)) - 1;
  std::vector<std::vector<uint16_t>> t_Counters(
    num_histArray, std::vector<uint16_t>(nts::nSamp * nts::rBuck, 0));

#ifdef _OPENMP
  omp_set_num_threads(nts::nThrd);
#endif

  //#pragma omp parallel for schedule(dynamic)
  for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
    std::vector<size_t> totKmer_vec(num_histArray, 0);
    ;
    btllib::SeqReader reader(inFiles[file_i],
                             btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
    for (const auto record : reader) {

      stRead(record.seq, t_Counters, totKmer_vec);

      /*for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
              size_t totKmer = 0;
              std::ifstream in(inFiles[file_i].c_str());
              std::string samSeq;
              unsigned ftype = getftype(in, samSeq);
              if (ftype == 0)
                      getEfq(in, t_Counter, totKmer);
              else if (ftype == 1)
                      getEfa(in, t_Counter, totKmer);
              else if (ftype == 2)
                      getEsm(in, samSeq, t_Counter, totKmer);
              else {
                      std::cerr << "Error in reading file: " << inFiles[file_i]
         << std::endl; exit(EXIT_FAILURE);
              }
              in.close();*/
    }
    for (size_t i = 0; i < num_histArray; ++i) {
#pragma omp atomic
      totalKmers_vec[i] += totKmer_vec[i];
    }
  }
  for (size_t curr_histArray = 0; curr_histArray < num_histArray;
       ++curr_histArray) {
    auto& histArray = histArray_vec[curr_histArray];
    double F0Mean = 0.0;
    double fMean[65536];
    auto& t_Counter = t_Counters[curr_histArray];
    compEst(t_Counter, F0Mean, fMean);
    histArray[0] = totalKmers_vec[curr_histArray];
    histArray[1] = (size_t)F0Mean;
    for (size_t i = 2; i <= nts::covMax + 1; i++)
      histArray[i] = (size_t)fMean[i - 1];
  }
#ifdef _OPENMP
  std::cerr << "Reapeat profile estimated using ntCard in (sec): "
            << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
#endif
  return histArray_vec;
}

uint64_t
calc_ntcard_genome_size(const std::string& input_path,
                        unsigned k,
                        const std::vector<std::string>& seed_string_vec,
                        unsigned threads)
{

  vector<string> inFiles;
  inFiles.push_back(input_path);
  uint64_t genome_size = 0;

  /* indices 0 and 1 are reserved for F0 and 1.
     Indicex 2 to 10001 store kmer multiplicities up to 10000 */
  static const size_t DEFAULT_NTCARD_HIST_COV_MAX = 10000;
  std::cerr << "Calculating expected entries" << std::endl;
  const auto histArray =
    getHist(inFiles, k, threads, DEFAULT_NTCARD_HIST_COV_MAX, seed_string_vec);
  for (size_t i = 0; i < seed_string_vec.size(); ++i) {
    std::cerr << "Expected entries for seed pattern " << seed_string_vec[i]
              << " : " << histArray[i][1] << std::endl;

    genome_size += (histArray[i][1]);
  }
  std::cerr << "Total expected entries for seed patterns: " << genome_size
            << std::endl;
  return genome_size;
}

#endif /* NTCARD_H_ */
