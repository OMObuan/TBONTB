// Copyright <2024> [OMObuan]

#include <TBONTB/rename_type.h>

#include <cstddef>
#include <random>

class Chromosome {
 public:
    friend std::ofstream &operator<<(std::ofstream &, Chromosome const &);

    friend std::ifstream &operator>>(std::ifstream &, Chromosome &);

    Chromosome() = default;

    explicit Chromosome(bool);

    Chromosome(f64, std::byte);

    void variateChromosome();

    std::byte getValue() const;

    f64 getDominance() const;

 private:
    std::byte m_value{};
    f64 m_dominance{};

    static std::default_random_engine engine;
};
