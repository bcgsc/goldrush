#ifndef CALC_PHRED_AVERAGE_H_
#define CALC_PHRED_AVERAGE_H_
#include <cstdint>
#include <string>

std::pair<uint32_t, uint32_t>
calc_phred_average(const std::string qual);

double
sum_phred(const std::string& qual);
#endif
