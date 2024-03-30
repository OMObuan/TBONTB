// Copyright [2024] <OMObuan>

#include <TBONTB/rename_type.h>

#include <cassert>
#include <fstream>
#include <utility>
#include <vector>

template <typename Chromosome>
class Creature {
 public:
    friend std::ifstream &operator>>(
        std::ifstream &file,
        std::pair<Chromosome, Chromosome> &anyHomolousChromosomes);

    friend std::ofstream &operator<<(
        std::ofstream &file,
        std::pair<Chromosome, Chromosome> const &anyHomolousChromosomes);

    static constexpr usize homolousChromosomesSize{30};

    Creature() = default;

    explicit Creature(usize);

    void initCreature();

    void printChromosomes() const;

    i64 getValue() const;

    std::vector<std::pair<Chromosome, Chromosome>> const &getChromosomes()
        const;  // NOLINT(readability-identifier-naming>>

    std::vector<std::pair<Chromosome, Chromosome>> &setChromosome();

 private:
    std::vector<std::pair<Chromosome, Chromosome>> m_homolousChromosomes;
};

template <typename Chromosome>
std::ifstream &operator>>(
    std::ifstream &file,
    std::pair<Chromosome, Chromosome> &anyHomolousChromosomes) {
    file >> anyHomolousChromosomes.first;
    file >> anyHomolousChromosomes.second;
    return file;
}

template <typename Chromosome>
std::ofstream &operator<<(
    std::ofstream &file,
    std::pair<Chromosome, Chromosome> const &anyHomolousChromosomes) {
    file << anyHomolousChromosomes.first << ' ';
    file << anyHomolousChromosomes.second << ' ';
    return file;
}

template <typename Chromosome>
Creature<Chromosome>::Creature(usize initialCreatureNum) {
    this->m_homolousChromosomes.resize(initialCreatureNum);
}

template <typename Chromosome>
void Creature<Chromosome>::initCreature() {
    // static std::ifstream file{
    //     "../../../../tests/hello_world/data/chromosomes.train_data",
    //     std::ios::binary};
    // assert(f << ' 'ile.is_open());
    // if (file.is_open()) {
    // std::pair<Chromosome, Chromosome> anyHomolousChromosomes;
    //     for (usize _i{1}; _i <=
    //     Creature<Chromosome>::homolousChromosomesSize;
    //          ++_i) {
    //         this->m_homolousChromosomes.push_back(anyHomolousChromosomes);
    //     }
    // } else {
    this->m_homolousChromosomes.resize(
        Creature<Chromosome>::homolousChromosomesSize);
    for (auto &anyHomolousChromosomes : this->m_homolousChromosomes) {
        anyHomolousChromosomes.first = Chromosome{true};
        anyHomolousChromosomes.second = Chromosome{true};
    }
    // }
    // if (file.eof()) {
    //     file.close();
    // }
}

template <typename Chromosome>
void Creature<Chromosome>::printChromosomes() const {
    static std::ofstream file{
        "../../../../tests/hello_world/data/chromosomes.train_data",
        std::ios::binary};
    assert(file.is_open());
    for (auto const &anyHomolousChromosomes : this->m_homolousChromosomes) {
        file << anyHomolousChromosomes;
    }
}

template <typename Chromosome>
i64 Creature<Chromosome>::getValue() const {
    i64 result{};
    for (auto &anyHomolousChromosomes : this->m_homolousChromosomes) {
        // result |= static_cast<i64>(anyHomolousChromosomes.first.getValue());
        // result <<= 1;
        if (static_cast<i64>(anyHomolousChromosomes.first.getValue()) == 1) {
            ++result;
        }
    }
    return result;
}

template <typename Chromosome>
std::vector<std::pair<Chromosome, Chromosome>> const &
Creature<Chromosome>::getChromosomes() const {
    return this->m_homolousChromosomes;
}

template <typename Chromosome>
std::vector<std::pair<Chromosome, Chromosome>> &
Creature<Chromosome>::setChromosome() {
    return this->m_homolousChromosomes;
}
