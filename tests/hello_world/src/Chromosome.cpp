// Copyright [2024] <OMObuan>

#include <Chromosome.h>

#include <ctime>

std::default_random_engine Chromosome::engine{
    static_cast<usize>(time(nullptr))};

Chromosome::Chromosome(f64 dominance, std::byte value)
    : m_value{value}, m_dominance{dominance} {}

std::byte Chromosome::getValue() const { return m_value; }

f64 Chromosome::getDominance() const { return m_dominance; }

void Chromosome::variateChormosome() {
    static std::uniform_int_distribution<usize> randomValueGenerator{0, 1};
    static std::uniform_real_distribution<f64> randomDominanceGenerator{0., 1.};

    m_dominance = randomDominanceGenerator(Chromosome::engine);
    m_value = static_cast<std::byte>(randomValueGenerator(Chromosome::engine));
}
