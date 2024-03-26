// Copyright <2024> [OMObuan]

#include <rename_type.h>

#include <cstddef>
#include <random>

class Chromosome {
 public:
    Chromosome() = default;

    Chromosome(f64, std::byte);

    void variateChormosome();

    std::byte getValue() const;
    f64 getDominance() const;

 private:
    std::byte m_value{};
    f64 m_dominance{};

    static std::default_random_engine engine;
};
