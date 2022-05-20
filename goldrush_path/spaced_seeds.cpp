#include "spaced_seeds.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

std::vector<std::string>
make_seed_pattern(const std::string& seed_preset,
                  unsigned k,
                  unsigned weight,
                  unsigned h)
{

  std::vector<std::string> seed_string_vec;
  std::string left_seed_str;
  std::string right_seed_str;

  if (seed_preset == "") {
    srand(time(NULL));
    // seed generation
    std::cerr << "Designing base symmetrical spaced seed"
              << "\n"
              << "Using:"
              << "\n"
              << "span: " << k << "\n"
              << "weight: " << weight << std::endl;

    std::vector<unsigned> left_seed_vec(k / 2, 0);
    left_seed_vec[0] = 1; // left most val in seed must be a 1
    size_t weight_count = 0;

    while (weight_count != weight / 2) {
      for (size_t i = 1; i < k / 2; ++i) {
        left_seed_vec[i] = rand() % 2;
      }
      weight_count = std::count(left_seed_vec.begin(), left_seed_vec.end(), 1);
    }

    std::stringstream temp_ss;
    for (const auto val : left_seed_vec) {
      temp_ss << val;
    }

    left_seed_str = temp_ss.str();
    right_seed_str = std::string(left_seed_str.rbegin(), left_seed_str.rend());

  } else {
    std::cerr << "Using preset spaced seed"
              << "\n"
              << "with:"
              << "\n"
              << "span: " << seed_preset.size() << "\n"
              << "weight: "
              << std::count(seed_preset.begin(), seed_preset.end(), '1')
              << std::endl;
    left_seed_str = seed_preset.substr(0, seed_preset.size() / 2);
    right_seed_str =
      seed_preset.substr(seed_preset.size() / 2, seed_preset.size() / 2);
  }

  for (size_t i = 0; i < h; ++i) {
    seed_string_vec.push_back(left_seed_str + std::string(i, '0') +
                              right_seed_str);
  }

  return seed_string_vec;
}