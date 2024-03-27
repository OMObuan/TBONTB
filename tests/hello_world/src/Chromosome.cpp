// Copyright [2024] <OMObuan>

#include <Chromosome.h>

#include <ctime>
#include <fstream>

std::ofstream &operator<<(std::ofstream &file, Chromosome const &chromosome) {
    file << chromosome.m_dominance << " "
         << static_cast<i32>(chromosome.m_value);
    return file;
}

std::ifstream &operator>>(std::ifstream &file, Chromosome &chromosome) {
    i16 inputByte{};
    file >> chromosome.m_dominance >> inputByte;
    chromosome.m_value = static_cast<std::byte>(inputByte);
    return file;
}

std::default_random_engine Chromosome::engine{
    static_cast<usize>(time(nullptr))};

Chromosome::Chromosome(f64 dominance, std::byte value)
    : m_value{value}, m_dominance{dominance} {}

Chromosome::Chromosome(bool random) {
    if (random) {
        variateChromosome();
    }
}

std::byte Chromosome::getValue() const { return m_value; }

f64 Chromosome::getDominance() const { return m_dominance; }

void Chromosome::variateChromosome() {
    static std::uniform_int_distribution<usize> randomValueGenerator{0, 1};
    static std::uniform_real_distribution<f64> randomDominanceGenerator{0., 1.};

    m_dominance = randomDominanceGenerator(Chromosome::engine);
    m_value = static_cast<std::byte>(randomValueGenerator(Chromosome::engine));
}
