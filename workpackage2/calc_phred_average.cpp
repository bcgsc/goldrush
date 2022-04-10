#include "calc_phred_average.hpp"
#include <string>


uint32_t calc_phred_average(const std::string qual) {
    uint32_t phred_sum = 0;
    for (char c : qual){ 
        phred_sum += ((uint32_t)c - 33);
    }
    return (phred_sum / qual.size());
}
