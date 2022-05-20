#ifndef SPACED_SEEDS_HPP
#define SPACED_SEEDS_HPP

#include <string>
#include <vector>

std::vector<std::string>
make_seed_pattern(const std::string& seed_preset,
                  unsigned k,
                  unsigned weight,
                  unsigned h);

#endif