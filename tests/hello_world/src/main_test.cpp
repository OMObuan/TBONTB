// Copyright [2024] <OMObuan>

#include <BasicData.h>
#include <Chromosome.h>
#include <TBONTB/rename_type.h>
#include <gtest/gtest.h>

#include <Creature.hpp>
#include <TBONTB/GeneticAlgorithm.hpp>
#include <algorithm>
#include <bitset>
#include <filesystem>
#include <ranges>

TEST(HelloWorldUnitTest, HelloWorld) {
    std::filesystem::path current_path{
        std::filesystem::current_path(),
    };
    // std::cerr << current_path << std::endl;
    std::ofstream file{"../../../../tests/hello_world/data/result.train_data",
                       std::ios::binary};
    ASSERT_TRUE(file.is_open());
    GeneticAlgorithm<Creature, Chromosome> geneticAlgorithm{100};
    geneticAlgorithm.initCreatures(&Creature<Chromosome>::initCreature);
    usize loopCnt{100};
    for (usize _i{1}; _i <= loopCnt; ++_i) {
        geneticAlgorithm.birthNewCreatures<BasicData>(
            &Creature<Chromosome>::getChromosomes,
            &Creature<Chromosome>::setChromosome, &Chromosome::getDominance,
            &Chromosome::variateChromosome, &Creature<Chromosome>::getValue);
        geneticAlgorithm.eliminateCreatures<BasicData>(
            &Creature<Chromosome>::getValue);
        geneticAlgorithm.printCreatures(
            &Creature<Chromosome>::printChromosomes);
    }
    std::vector<usize> x{};
    std::ranges::sort(
        geneticAlgorithm.setPopulation(),
        [](Creature<Chromosome>& lft, Creature<Chromosome>& rgt) -> bool {
            return lft.getValue() > rgt.getValue();
        });
    ASSERT_GE(geneticAlgorithm.getPopulation()[0].getValue(), 28);
    // for (auto const& anyCreature : geneticAlgorithm.getPopulation()) {
    //     file << anyCreature.getValue() << '\n';
    // }
}
