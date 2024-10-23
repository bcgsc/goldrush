#include "calc_phred_average.hpp"

#include <cmath>
#include <cstdlib> // for abs
#include <string>
#include <utility>

std::pair<uint32_t, uint32_t>
calc_phred_average(const std::string qual)
{
  double phred_sum = 0.0;
  double first_avg = 0.0;
  double second_avg = 0.0;
  size_t qual_size = qual.size();

  for (size_t i = 0; i < qual_size; ++i) {
    // Convert ASCII character to Phred score
    int phred_score = (int)qual.at(i) - 33;

    // Delog the Phred score: 10^(-Q/10)
    double delog_phred = pow(10.0, -phred_score / 10.0);

    phred_sum += delog_phred;

    // Store sum for the first half
    if (i == qual_size / 2 - 1) {
      first_avg = phred_sum;
    }
  }

  // Compute the second half average
  second_avg = phred_sum - first_avg;

  // Calculate averages by dividing by half of the length
  second_avg = second_avg / (qual_size * 0.5);
  first_avg = first_avg / (qual_size * 0.5);

  // Return pair: total average Phred score and absolute difference between
  // first and second halves
  return std::make_pair((uint32_t)(-10 * log10(phred_sum / qual_size)),
                        (uint32_t)abs((int32_t)(-10 * log10(first_avg)) -
                                      (int32_t)(-10 * log10(second_avg))));
}