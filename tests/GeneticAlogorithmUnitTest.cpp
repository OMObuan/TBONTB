// Copyright [2024] <OMObuan>

#include <gtest/gtest.h>

#include <GeneticAlogorithm.hpp>

namespace TEST_FOR_INIT_CREATURES {
template <typename Chromosome>
struct Creature {
 public:
    Chromosome val{-1};
    Creature() = default;
    void initCreatures() { this->val = Chromosome{0}; }
};

TEST(initCreaturesUnitTest, Initialization) {
    GeneticAlogorithm<Creature, i32> testObject(250);
    ASSERT_EQ(testObject.getCreatures().size(), 250);
    testObject.initCreatures(&Creature<i32>::initCreatures);
    for (auto &anyCreature : testObject.getCreatures()) {
        ASSERT_EQ(anyCreature.val, 0);
    }
}

}  // namespace TEST_FOR_INIT_CREATURES

namespace TEST_FOR_GETELIMATED_SIZE {

template <typename Chromosome>
struct Creature {};
struct ConstDate {
    static constexpr double eliminatedRate = 0.5;
};

TEST(GetElimatedSizeUnitTest, CalculationRight) {
    ASSERT_EQ(
        (GeneticAlogorithm<Creature, i32>::getEliminatedSize<ConstDate>(250)),
        (static_cast<usize>(125)));
    ASSERT_EQ(
        (GeneticAlogorithm<Creature, i32>::getEliminatedSize<ConstDate>(251)),
        (static_cast<usize>(126)));
}

}  // namespace TEST_FOR_GETELIMATED_SIZE

namespace TEST_FOR_ELIMINATE_CREATURE {

template <typename Chromosome>
struct Creature {
 private:
    static Chromosome nowVal;

 public:
    Chromosome val{};
    void initCreatures() { this->val = nowVal--; }
    Chromosome getValue() const { return val; }
};

template <typename Chromosome>
Chromosome Creature<Chromosome>::nowVal{100};

struct ConstDate {
    static constexpr double eliminatedRate = 0.75;
};

TEST(EliminateCreaturesUnitTest, EliminationRight) {
    GeneticAlogorithm<Creature, i32> testObject(100);
    testObject.initCreatures(&Creature<i32>::initCreatures);
    testObject.eliminateCreatures<ConstDate>(&Creature<i32>::getValue);
    ASSERT_EQ(testObject.getCreatures().size(), 25);
    // i32 nowVal = 100;
    for (auto &anyCreature : testObject.getCreatures()) {
        // std::cerr << anyCreature.getValue() << std::endl;
        ASSERT_GE(anyCreature.getValue(), 76);
        // --nowVal;
    }
}

}  // namespace TEST_FOR_ELIMINATE_CREATURE

namespace TEST_FOR_QUICK_NATRUE_NUMBER_POW {

template <typename Chromosome>
struct Creature {};

TEST(QuickNatureNumberPowUnitTest, CalculationRight) {}
}  // namespace TEST_FOR_QUICK_NATRUE_NUMBER_POW
