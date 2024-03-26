// Copyright [2024] <OMObuan>

#include <rename_type.h>

template <typename Chromosome>
class Creature {
 public:
    Creature() = default;

 private:
    std::vector<std::pair<Chromosome, Chromosome>> m_homolousChromosomes;
};
