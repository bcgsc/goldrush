#include "calc_phred_average.hpp"
#include <string>


std::pair<uint32_t, uint32_t> calc_phred_average(const std::string qual) {
    uint32_t phred_sum = 0;
    /*for (char c : qual){ 
        phred_sum += (uint32_t)c ;
    }*/
    uint32_t first_avg = 0;
    uint32_t second_avg = 0;
    for (size_t i = 0; i < qual.size(); ++i) {
        phred_sum += (uint32_t)qual.at(i);
        if (i == qual.size() / 2) {
            first_avg= phred_sum;
        }
    }
    second_avg = phred_sum - first_avg;
    second_avg = second_avg / (qual.size() * 0.5);
    first_avg = first_avg / (qual.size() * 0.5);
    return std::make_pair((phred_sum / qual.size()) - 33, (uint32_t)abs((int32_t)first_avg - (int32_t)second_avg));
}
