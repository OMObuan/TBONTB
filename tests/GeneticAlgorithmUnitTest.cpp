// Copyright [2024] <OMObuan>

#include <debug.h>
#include <gtest/gtest.h>
#include <math_calculate.h>

#include <GeneticAlgorithm.hpp>
#include <array>
#include <random>
#include <utility>

namespace TEST_FOR_INIT_CREATURES {
template <typename Chromosome>
struct Creature {
 public:
    Chromosome val{-1};
    Creature() = default;
    explicit Creature(usize);
    void initCreature() { this->val = Chromosome{0}; }
};

TEST(GeneticAlgorithmUnitTest, initCreaturesUnitTest) {
    GeneticAlgorithm<Creature, i32> testObject(250);
    ASSERT_EQ(testObject.getCreatures().size(), 250);
    testObject.initCreatures(&Creature<i32>::initCreature);
    for (auto &anyCreature : testObject.getCreatures()) {
        ASSERT_EQ(anyCreature.val, 0);
    }
}

}  // namespace TEST_FOR_INIT_CREATURES

namespace TEST_FOR_GETELIMATED_SIZE {

template <typename Chromosome>
struct Creature {
    Creature() = default;
    explicit Creature(usize) {}
};
struct BasicDate {
    static constexpr double eliminatedRate = 0.5;
};

TEST(GeneticAlgorithmUnitTest, GetElimatedSizeUnitTest) {
    ASSERT_EQ(
        (GeneticAlgorithm<Creature, i32>::getEliminatedSize<BasicDate>(250)),
        (static_cast<usize>(125)));
    ASSERT_EQ(
        (GeneticAlgorithm<Creature, i32>::getEliminatedSize<BasicDate>(251)),
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

struct BasicDate {
    static constexpr double eliminatedRate = 0.75;
};

TEST(GeneticAlgorithmUnitTest, EliminateCreaturesUnitTest) {
    GeneticAlgorithm<Creature, i32> testObject(100);
    testObject.initCreatures(&Creature<i32>::initCreatures);
    testObject.eliminateCreatures<BasicDate>(&Creature<i32>::getValue);
    ASSERT_EQ(testObject.getCreatures().size(), 25);
    // i32 nowVal = 100;
    for (auto &anyCreature : testObject.getCreatures()) {
        // std::cerr << anyCreature.getValue() << std::endl;
        ASSERT_GE(anyCreature.getValue(), 76);
        // --nowVal;
    }
}

}  // namespace TEST_FOR_ELIMINATE_CREATURE

namespace TEST_FOR_GET_VALUE_SUM {

template <typename Chromosome>
struct Creature {
 public:
    Chromosome getValue() const { return val; }

    void initCreature() { val = ++initialVal; }

    Chromosome val{};
    static Chromosome initialVal;
};

template <typename Chromosome>
Chromosome Creature<Chromosome>::initialVal{};

TEST(GeneticAlgorithmUnitTest, GetValueSumUnitTest) {
    GeneticAlgorithm<Creature, i32> testObject(100);
    testObject.initCreatures(&Creature<i32>::initCreature);
    ASSERT_EQ(testObject.getValueSum(&Creature<i32>::getValue), 5050);
}

}  // namespace TEST_FOR_GET_VALUE_SUM

namespace TEST_FOR_GENERATE_RAND_ARR {

template <typename Chromosome>
struct Creature {
 public:
    void initCreature() {
        if (initalID == 0) {
            val = Chromosome{100};
        } else {
            val = Chromosome{1};
        }
        ++initalID;
    }

    Chromosome getValue() const { return val; }

    Chromosome val{};
    static i32 initalID;
};

template <typename Chromosome>
i32 Creature<Chromosome>::initalID{};

struct BasicDate {
    static constexpr usize probabilityPrecision = 2;
};

TEST(GeneticAlgorithmUnitTest, GenerateRandArryUnitTest) {
    GeneticAlgorithm<Creature, i32> testObject(101);
    testObject.initCreatures(&Creature<i32>::initCreature);
    auto result{
        testObject.generateRandArray<BasicDate>(&Creature<i32>::getValue)};
    std::array<usize, 101> testArr{};
    for (auto anyElement : result) {
        ++testArr[anyElement];
    }
    // for (auto anyElement : testArr) {
    // std::cerr << anyElement << '\n';
    // }
    ASSERT_EQ(testArr[0], 50);
    for (usize i = 1; i < 101; ++i) {
        ASSERT_EQ(testArr[i], 1);
    }
}

}  // namespace TEST_FOR_GENERATE_RAND_ARR

namespace TEST_FOR_GENERATE_RAND_CREATURE {

template <typename Chromosome>
struct Creature {
 public:
    void initialValue() {
        if (initalID < 50) {
            val = 2;
        } else {
            val = 1;
        }
        ++initalID;
    }

    Chromosome getValue() const { return val; }

    Chromosome val{};
    static i32 initalID;
};

template <typename Chromosome>
i32 Creature<Chromosome>::initalID{};

struct BasicDate {
    static constexpr usize probabilityPrecision = 5;
};

TEST(GeneticAlgorithmUnitTest, GenerateRandCreatureUnitTest) {
    GeneticAlgorithm<Creature, i32> testObject(100);
    testObject.initCreatures(&Creature<i32>::initialValue);
    auto result{
        testObject.generateRandArray<BasicDate>(&Creature<i32>::getValue)};
    usize loopCnt{100000};
    usize twoCnt{};
    while (loopCnt--) {
        if (testObject.getRandCreature(result).val == 2) {
            ++twoCnt;
        }
    }
    ASSERT_LE(twoCnt, static_cast<usize>(
                          0.75 * MATH_CALCULATE::quickNatureNumberPow(
                                     10, BasicDate::probabilityPrecision)));
    ASSERT_GE(twoCnt, static_cast<usize>(
                          0.25 * MATH_CALCULATE::quickNatureNumberPow(
                                     10, BasicDate::probabilityPrecision)));
}

}  // namespace TEST_FOR_GENERATE_RAND_CREATURE

namespace TEST_FOR_GET_CHROMOSOME_FROM_BOOL {

template <typename>
struct Creature {
    explicit Creature(usize);
};

TEST(GeneticAlgorithmUnitTest, GetChromosomeFromBoolUnitTest) {
    ASSERT_EQ((GeneticAlgorithm<Creature, i32>::getChormosomeFromBool(
                  std::make_pair(0, 1), true)),
              0);
    ASSERT_EQ((GeneticAlgorithm<Creature, i32>::getChormosomeFromBool(
                  std::make_pair(0, 1), false)),
              1);
}

}  // namespace TEST_FOR_GET_CHROMOSOME_FROM_BOOL

namespace TEST_FOR_MIX_CHROMOSOME {

template <typename Chromosome>
struct Creature {
    explicit Creature(usize);
};

struct Chromosome {
 public:
    explicit Chromosome(i32 dominance, i32 val)
        : dominance{dominance}, val{val} {}

    bool operator==(const Chromosome &other) const {
        return (val == other.val && dominance == other.dominance);
    }

    i32 getDominance() const { return dominance; }

    i32 dominance{};
    i32 val{};
};

TEST(GeneticAlgorithmUnitTest, MixChromosomeUnitTest) {
    ASSERT_EQ(
        (GeneticAlgorithm<Creature, Chromosome>::mixChromosome(
            Chromosome{1, 0}, Chromosome{2, 1}, &Chromosome::getDominance)),
        std::make_pair(Chromosome{2, 1}, Chromosome{1, 0}));
}

}  // namespace TEST_FOR_MIX_CHROMOSOME

namespace TEST_FOR_VARIATE_CHROMOSOME {

template <typename Chromosome>
struct Creature {
    explicit Creature(usize);
};

struct Chromosome {
 public:
    void variateChromosome() {
        std::uniform_int_distribution<usize> randomValueGenerator{0, 1};
        std::uniform_real_distribution<f64> randomDominanceGenerator{0., 1.};
        value = randomValueGenerator(Chromosome::engine);
        dominance = randomDominanceGenerator(Chromosome::engine);
    }

    i32 value{};
    f64 dominance{};
    static std::default_random_engine engine;
};

std::default_random_engine Chromosome::engine{
    static_cast<usize>(time(nullptr))};

TEST(GeneticAlgorithmUnitTest, VariateChromosomeUnitTest) {
    usize loopCnt =
        static_cast<usize>(MATH_CALCULATE::quickNatureNumberPow(10, 5));
    usize zeroCnt{}, lessHalfCnt{};
    for (usize _i{1}; _i <= loopCnt; ++_i) {
        Chromosome anyChromosome{};
        GeneticAlgorithm<Creature, Chromosome>::variateChromosome(
            &anyChromosome, &Chromosome::variateChromosome);
        if (anyChromosome.value == 0) {
            ++zeroCnt;
        }
        if (MATH_CALCULATE::f64Less(anyChromosome.dominance, 0.5, 0.00001)) {
            ++lessHalfCnt;
        }
    }
    // _LOOK(std::cerr,
    //       typeid(MATH_CALCULATE::quickNatureNumberPow(10, 5)).name());
    // _LOOK(std::cerr,
    //       static_cast<usize>(MATH_CALCULATE::quickNatureNumberPow(10, 5)));
    // _LOOK(std::cerr, loopCnt);
    ASSERT_LE(zeroCnt, static_cast<usize>(0.6 * loopCnt));
    ASSERT_GE(zeroCnt, static_cast<usize>(0.4 * loopCnt));
    ASSERT_LE(lessHalfCnt, static_cast<usize>(0.6 * loopCnt));
    ASSERT_GE(lessHalfCnt, static_cast<usize>(0.4 * loopCnt));
}

}  // namespace TEST_FOR_VARIATE_CHROMOSOME

namespace TEST_FOR_MATE_CREATURE {

struct Chromosome {
 public:
    Chromosome() = default;

    Chromosome(f64 dominance, i32 val) : dominance{dominance}, value{val} {}

    f64 getDominance() const { return dominance; }

    i32 getValue() const { return value; }

    f64 dominance{};
    i32 value{};
};

template <typename Chromosome>
struct Creature {
 public:
    Creature() = default;

    explicit Creature(usize chromosomeSize)
        : homolousChromosomes{chromosomeSize} {}

    std::vector<std::pair<Chromosome, Chromosome>> const &getChromosomes() {
        return homolousChromosomes;
    }

    void setChromosome(usize id,
                       std::pair<Chromosome, Chromosome> anyChromosome) {
        homolousChromosomes[id] = anyChromosome;
    }

    void initCreature() { homolousChromosomes = {{{1., 1}, {0., 0}}}; }
    std::vector<std::pair<Chromosome, Chromosome>> homolousChromosomes{};
};

struct BasicDate {
    static constexpr f64 variationRate{};
};

TEST(GeneticAlgorithmUnitTest, MateCreatureUnitTest) {
    usize loop{static_cast<usize>(MATH_CALCULATE::quickNatureNumberPow(10, 5))};
    for (usize _i{1}; _i <= loop; ++_i) {
        GeneticAlgorithm<Creature, Chromosome>::mateCreature<BasicDate>(
            Creature<Chromosome>{})
    }
}

}  // namespace TEST_FOR_MATE_CREATURE
