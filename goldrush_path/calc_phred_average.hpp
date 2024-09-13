#ifndef CALC_PHRED_AVERAGE_H_
#define CALC_PHRED_AVERAGE_H_
#include <string>
#include <cstdint>

std::pair<uint32_t, uint32_t>
calc_phred_average(const std::string qual);
#endif
