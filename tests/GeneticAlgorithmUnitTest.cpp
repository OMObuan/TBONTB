// Copyright [2024] <OMObuan>

#include <debug.h>
#include <gtest/gtest.h>
#include <math_calculate.h>

#include <GeneticAlgorithm.hpp>
#include <array>
#include <cassert>
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
    ASSERT_EQ(testObject.getPopulation().size(), 250);
    testObject.initCreatures(&Creature<i32>::initCreature);
    for (auto &anyCreature : testObject.getPopulation()) {
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
    ASSERT_EQ(testObject.getPopulation().size(), 25);
    // i32 nowVal = 100;
    for (auto &anyCreature : testObject.getPopulation()) {
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

std::default_random_engine g_engine{static_cast<usize>(std::time(nullptr))};

TEST(GeneticAlgorithmUnitTest, GenerateRandArryUnitTest) {
    usize loopCnt{MATH_CALCULATE::quickNatureNumberPow(10, 5)}, zeroCnt{};
    GeneticAlgorithm<Creature, i32> testObject(101);
    testObject.initCreatures(&Creature<i32>::initCreature);
    auto result{testObject.generateRandArray(&Creature<i32>::getValue)};
    for (usize _i{1}; _i <= loopCnt; ++_i) {
        usize value{static_cast<usize>(result(g_engine))};
        if (value == 0) {
            ++zeroCnt;
        }
    }
    ASSERT_LE(zeroCnt, static_cast<usize>(loopCnt * 0.6));
    ASSERT_GE(zeroCnt, static_cast<usize>(loopCnt * 0.4));
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
    auto result{testObject.generateRandArray(&Creature<i32>::getValue)};
    usize loopCnt{100000};
    usize twoCnt{};
    while (loopCnt--) {
        if (testObject.getRandCreature(&result).val == 2) {
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

    void variateChromosome() {}

    f64 dominance{};
    i32 value{};
    static usize variateCnt;
};

usize Chromosome::variateCnt{};

template <typename Chromosome>
struct Creature {
 public:
    Creature() = default;

    explicit Creature(usize chromosomeSize)
        : homolousChromosomes{chromosomeSize} {}

    std::vector<std::pair<Chromosome, Chromosome>> const &getChromosomes()
        const {
        return homolousChromosomes;
    }

    std::vector<std::pair<Chromosome, Chromosome>> &setChromosome() {
        return homolousChromosomes;
    }

    void initCreatureForOnePair() {
        homolousChromosomes = {{{1., 1}, {0., 0}}};
    }

    void initCreatureForTowPair() {
        homolousChromosomes = {{{1., 1}, {0., 0}}, {{1., 1}, {0., 0}}};
    }

    std::vector<std::pair<Chromosome, Chromosome>> homolousChromosomes{};
};

struct BasicDate {
    static constexpr f64 variationRate{};
};

TEST(GeneticAlgorithmUnitTest, MateCreatureUnitTest) {
    usize loop{static_cast<usize>(MATH_CALCULATE::quickNatureNumberPow(10, 5))},
        dominanceCnt{}, towDominanceCnt{}, formerOneDominanceCnt{},
        latterOneDominanceCnt{};
    // usize loop{1};
    Creature<Chromosome> hybrid{};
    hybrid.initCreatureForOnePair();
    for (usize _i{1}; _i <= loop; ++_i) {
        auto result =
            GeneticAlgorithm<Creature, Chromosome>::mateCreature<BasicDate>(
                hybrid, hybrid, &Creature<Chromosome>::getChromosomes,
                &Creature<Chromosome>::setChromosome, &Chromosome::getDominance,
                &Chromosome::variateChromosome);
        if (result.getChromosomes()[0].first.getValue() == 1) {
            ++dominanceCnt;
        }
    }
    ASSERT_LE(dominanceCnt, static_cast<usize>(loop * 0.8));
    ASSERT_GE(dominanceCnt, static_cast<usize>(loop * 0.7));
    hybrid.initCreatureForTowPair();
    for (usize _i{1}; _i <= loop; ++_i) {
        auto result =
            GeneticAlgorithm<Creature, Chromosome>::mateCreature<BasicDate>(
                hybrid, hybrid, &Creature<Chromosome>::getChromosomes,
                &Creature<Chromosome>::setChromosome, &Chromosome::getDominance,
                &Chromosome::variateChromosome);
        if (result.getChromosomes()[0].first.getValue() == 1 &&
            result.getChromosomes()[1].first.getValue() == 1) {
            ++towDominanceCnt;
        } else if (result.getChromosomes()[0].first.getValue() == 1 &&
                   result.getChromosomes()[1].first.getValue() == 0) {
            ++formerOneDominanceCnt;
        } else if (result.getChromosomes()[0].first.getValue() == 0 &&
                   result.getChromosomes()[1].first.getValue() == 1) {
            ++latterOneDominanceCnt;
        }
    }
    ASSERT_LE(towDominanceCnt, static_cast<usize>(loop * 0.6));
    ASSERT_GE(towDominanceCnt, static_cast<usize>(loop * 0.5));
    ASSERT_LE(formerOneDominanceCnt, static_cast<usize>(loop * 0.2));
    ASSERT_GE(formerOneDominanceCnt, static_cast<usize>(loop * 0.1));
    ASSERT_LE(latterOneDominanceCnt, static_cast<usize>(loop * 0.2));
    ASSERT_GE(latterOneDominanceCnt, static_cast<usize>(loop * 0.1));
}

}  // namespace TEST_FOR_MATE_CREATURE

namespace TEST_FOR_BIRTH_NEW_CREATURES {

struct Chromosome {
 public:
    Chromosome() = default;

    Chromosome(f64 dominance, i32 val) : dominance{dominance}, value{val} {}

    f64 getDominance() const { return dominance; }

    i32 getValue() const { return value; }

    void variateChromosome() {
        dominance = 0.;
        value = 0;
        // if (MATH_CALCULATE::f64Equal(dominance, 1., 0.0001)) {
        //     dominance = 0.;
        //     ++variateCnt;
        // }
    }

    f64 dominance{};
    i32 value{};
    static usize variateCnt;
};

usize Chromosome::variateCnt{};

template <typename Chromosome>
struct Creature {
 public:
    Creature() = default;

    explicit Creature(usize chromosomeSize)
        : homolousChromosomes{chromosomeSize} {}

    Creature(Creature<Chromosome> const &creature) {
        homolousChromosomes = creature.getChromosomes();
    }

    Creature(Creature<Chromosome> &&creature) {
        homolousChromosomes = std::move(creature.getChromosomes());
    }

    std::vector<std::pair<Chromosome, Chromosome>> const &getChromosomes()
        const {
        return homolousChromosomes;
    }

    std::vector<std::pair<Chromosome, Chromosome>> &setChromosome() {
        return homolousChromosomes;
    }

    void initCreature() { homolousChromosomes = {{{1., 1}, {0., 0}}}; }

    i32 getValue() const { return 1; }

    std::vector<std::pair<Chromosome, Chromosome>> homolousChromosomes{};
};

struct BasicDate {
    static constexpr f64 variationRate{0.5};
    static constexpr f64 increasedCreatureSize{1.};
};

TEST(GeneticAlgorithmUnitTest, BirthNewCreauresUnitTest) {
    GeneticAlgorithm<Creature, Chromosome> testObject{
        MATH_CALCULATE::quickNatureNumberPow(10, 5)};
    testObject.initCreatures(&Creature<Chromosome>::initCreature);
    testObject.birthNewCreatures<BasicDate>(
        &Creature<Chromosome>::getChromosomes,
        &Creature<Chromosome>::setChromosome, &Chromosome::getDominance,
        &Chromosome::variateChromosome, &Creature<Chromosome>::getValue);
    usize populationSize{testObject.getPopulation().size()}, oneCnt{};
    // _LOOK(std::cerr, populationSize);
    // _LOOK(std::cerr, Chromosome::variateCnt);
    // usize AA{};
    for (usize _i{MATH_CALCULATE::quickNatureNumberPow(10, 5)};
         _i < populationSize; ++_i) {
        ASSERT_TRUE(MATH_CALCULATE::f64Equal(
            testObject.getPopulation()[_i].getChromosomes()[0].first.getValue(),
            testObject.getPopulation()[_i]
                .getChromosomes()[0]
                .first.getDominance(),
            0.00001));
        if (testObject.getPopulation()[_i]
                .getChromosomes()[0]
                .first.getValue() == 1) {
            ++oneCnt;
        }
    }
    // _LOOK(std::cerr, AA);
    ASSERT_LE(oneCnt, static_cast<usize>(8. / 16. / 2. * populationSize));
    ASSERT_GE(oneCnt, static_cast<usize>(6. / 16. / 2. * populationSize));
}

}  // namespace TEST_FOR_BIRTH_NEW_CREATURES
