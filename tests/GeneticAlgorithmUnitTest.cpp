// Copyright [2024] <OMObuan>

#include <TBONTB/debug.h>
#include <TBONTB/math_calculate.h>
#include <gtest/gtest.h>

#include <TBONTB/GeneticAlgorithm.hpp>
#include <random>
#include <utility>

namespace TEST_FOR_INIT_CREATURES {
template <typename Chromosome>
struct Creature {
 public:
    Creature() = default;

    explicit Creature(usize);

    void initCreature();

    Chromosome getValue() const;

 private:
    Chromosome m_val{-1};
};

template <typename Chromosome>
void Creature<Chromosome>::initCreature() {
    this->m_val = Chromosome{0};
}

template <typename Chromosome>
Chromosome Creature<Chromosome>::getValue() const {
    return m_val;
}

TEST(GeneticAlgorithmUnitTest, initCreaturesUnitTest) {
    GeneticAlgorithm<Creature, i32> testObject(250);
    ASSERT_EQ(testObject.getPopulation().size(), 250);
    testObject.initCreatures(&Creature<i32>::initCreature);
    for (auto &anyCreature : testObject.getPopulation()) {
        ASSERT_EQ(anyCreature.getValue(), 0);
    }
}

}  // namespace TEST_FOR_INIT_CREATURES

namespace TEST_FOR_GETELIMATED_SIZE {

template <typename Chromosome>
struct Creature {
    Creature() = default;

    explicit Creature(usize) {}
};

struct BasicData {
    static constexpr double eliminatedRate{0.5};
};

TEST(GeneticAlgorithmUnitTest, GetElimatedSizeUnitTest) {
    ASSERT_EQ(
        (GeneticAlgorithm<Creature, i32>::getEliminatedSize<BasicData>(250)),
        (static_cast<usize>(125)));
    ASSERT_EQ(
        (GeneticAlgorithm<Creature, i32>::getEliminatedSize<BasicData>(251)),
        (static_cast<usize>(126)));
}

}  // namespace TEST_FOR_GETELIMATED_SIZE

namespace TEST_FOR_ELIMINATE_CREATURE {

template <typename Chromosome>
struct Creature {
 public:
    Creature() = default;

    explicit Creature(usize);

    void initCreatures();

    Chromosome getValue() const;

 private:
    static Chromosome m_nowVal;

    Chromosome m_val{};
};

template <typename Chromosome>
void Creature<Chromosome>::initCreatures() {
    this->m_val = m_nowVal--;
}

template <typename Chromosome>
Chromosome Creature<Chromosome>::getValue() const {
    return m_val;
}

template <typename Chromosome>
Chromosome Creature<Chromosome>::m_nowVal{100};

struct BasicData {
    static constexpr double eliminatedRate = 0.75;
};

TEST(GeneticAlgorithmUnitTest, EliminateCreaturesUnitTest) {
    GeneticAlgorithm<Creature, i32> testObject(100);
    testObject.initCreatures(&Creature<i32>::initCreatures);
    testObject.eliminateCreatures<BasicData>(&Creature<i32>::getValue);
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
    Creature() = default;

    explicit Creature(usize);

    Chromosome getValue() const;

    void initCreature();

 private:
    Chromosome m_val{};

    static Chromosome m_initialVal;
};

template <typename Chromosome>
Chromosome Creature<Chromosome>::m_initialVal{};

template <typename Chromosome>
Chromosome Creature<Chromosome>::getValue() const {
    return m_val;
}

template <typename Chromosome>
void Creature<Chromosome>::initCreature() {
    m_val = ++m_initialVal;
}

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
    Creature() = default;

    explicit Creature(usize) {}

    void initCreature();

    Chromosome getValue() const;

 private:
    Chromosome m_val{};

    static i32 m_initalID;
};

template <typename Chromosome>
i32 Creature<Chromosome>::m_initalID{};

template <typename Chromosome>
void Creature<Chromosome>::initCreature() {
    if (Creature<Chromosome>::m_initalID == 0) {
        this->m_val = Chromosome{100};
    } else {
        this->m_val = Chromosome{1};
    }
    ++Creature<Chromosome>::m_initalID;
}

template <typename Chromosome>
Chromosome Creature<Chromosome>::getValue() const {
    return this->m_val;
}

struct BasicData {
    static constexpr usize probabilityPrecision{2};
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
    Creature() = default;

    explicit Creature(usize);

    void initialValue();

    Chromosome getValue() const;

 private:
    Chromosome m_val{};

    static i32 m_initalID;
};

template <typename Chromosome>
i32 Creature<Chromosome>::m_initalID{};

template <typename Chromosome>
void Creature<Chromosome>::initialValue() {
    if (m_initalID < 50) {
        this->m_val = 2;
    } else {
        this->m_val = 1;
    }
    ++Creature<Chromosome>::m_initalID;
}

template <typename Chromosome>
Chromosome Creature<Chromosome>::getValue() const {
    return this->m_val;
}

struct BasicData {
    static constexpr usize probabilityPrecision = 5;
};

TEST(GeneticAlgorithmUnitTest, GenerateRandCreatureUnitTest) {
    GeneticAlgorithm<Creature, i32> testObject(100);
    testObject.initCreatures(&Creature<i32>::initialValue);
    auto result{testObject.generateRandArray(&Creature<i32>::getValue)};
    usize loopCnt{100000};
    usize twoCnt{};
    while (loopCnt--) {
        if (testObject.getRandCreature(&result).getValue() == 2) {
            ++twoCnt;
        }
    }
    ASSERT_LE(twoCnt, static_cast<usize>(
                          0.75 * MATH_CALCULATE::quickNatureNumberPow(
                                     10, BasicData::probabilityPrecision)));
    ASSERT_GE(twoCnt, static_cast<usize>(
                          0.25 * MATH_CALCULATE::quickNatureNumberPow(
                                     10, BasicData::probabilityPrecision)));
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
    explicit Chromosome(i32, i32);

    bool operator==(const Chromosome &) const;

    i32 getDominance() const;

 private:
    i32 m_dominance{};

    i32 m_val{};
};

Chromosome::Chromosome(i32 dominance, i32 val)
    : m_dominance{dominance}, m_val{val} {}

bool Chromosome::operator==(const Chromosome &other) const {
    return (this->m_val == other.m_val &&
            this->m_dominance == other.m_dominance);
}

i32 Chromosome::getDominance() const { return this->m_dominance; }

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
    void variateChromosome();

    i32 getValue() const;

    f64 getDominance() const;

 private:
    i32 m_value{};

    f64 m_dominance{};

    static std::default_random_engine m_engine;
};

std::default_random_engine Chromosome::m_engine{
    static_cast<usize>(time(nullptr))};

void Chromosome::variateChromosome() {
    std::uniform_int_distribution<usize> randomValueGenerator{0, 1};
    std::uniform_real_distribution<f64> randomDominanceGenerator{0., 1.};
    this->m_value = randomValueGenerator(Chromosome::m_engine);
    this->m_dominance = randomDominanceGenerator(Chromosome::m_engine);
}

i32 Chromosome::getValue() const { return this->m_value; }

f64 Chromosome::getDominance() const { return this->m_dominance; }

TEST(GeneticAlgorithmUnitTest, VariateChromosomeUnitTest) {
    usize loopCnt =
        static_cast<usize>(MATH_CALCULATE::quickNatureNumberPow(10, 5));
    usize zeroCnt{}, lessHalfCnt{};
    for (usize _i{1}; _i <= loopCnt; ++_i) {
        Chromosome anyChromosome{};
        GeneticAlgorithm<Creature, Chromosome>::variateChromosome(
            &anyChromosome, &Chromosome::variateChromosome);
        if (anyChromosome.getValue() == 0) {
            ++zeroCnt;
        }
        if (MATH_CALCULATE::f64Less(anyChromosome.getDominance(), 0.5,
                                    0.00001)) {
            ++lessHalfCnt;
        }
    }

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

    Chromosome(f64, i32);

    f64 getDominance() const;

    i32 getValue() const;

    void variateChromosome();

 private:
    f64 m_dominance{};

    i32 m_value{};

    static usize m_variateCnt;
};

usize Chromosome::m_variateCnt{};

Chromosome::Chromosome(f64 dominance, i32 val)
    : m_dominance{dominance}, m_value{val} {}

f64 Chromosome::getDominance() const { return this->m_dominance; }

i32 Chromosome::getValue() const { return this->m_value; }

void Chromosome::variateChromosome() {}

template <typename Chromosome>
struct Creature {
 public:
    Creature() = default;

    explicit Creature(usize);

    std::vector<std::pair<Chromosome, Chromosome>> const &getChromosomes()
        const;

    std::vector<std::pair<Chromosome, Chromosome>> &setChromosome();

    void initCreatureForOnePair();

    void initCreatureForTowPair();

 private:
    std::vector<std::pair<Chromosome, Chromosome>> m_homolousChromosomes{};
};

template <typename Chromosome>
Creature<Chromosome>::Creature(usize chromosomeSize)
    : m_homolousChromosomes{chromosomeSize} {}

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

template <typename Chromosome>
void Creature<Chromosome>::initCreatureForOnePair() {
    this->m_homolousChromosomes = {{{1., 1}, {0., 0}}};
}

template <typename Chromosome>
void Creature<Chromosome>::initCreatureForTowPair() {
    this->m_homolousChromosomes = {{{1., 1}, {0., 0}}, {{1., 1}, {0., 0}}};
}

struct BasicData {
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
            GeneticAlgorithm<Creature, Chromosome>::mateCreature<BasicData>(
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
            GeneticAlgorithm<Creature, Chromosome>::mateCreature<BasicData>(
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

    Chromosome(f64, i32);

    f64 getDominance() const;

    i32 getValue() const;

    void variateChromosome();

 private:
    f64 m_dominance{};

    i32 m_value{};

    static usize m_variateCnt;
};

usize Chromosome::m_variateCnt{};

Chromosome::Chromosome(f64 dominance, i32 val)
    : m_dominance{dominance}, m_value{val} {}

f64 Chromosome::getDominance() const { return this->m_dominance; }

i32 Chromosome::getValue() const { return this->m_value; }

void Chromosome::variateChromosome() {
    this->m_dominance = 0.;
    this->m_value = 0;
}

template <typename Chromosome>
struct Creature {
 public:
    Creature() = default;

    explicit Creature(usize);

    Creature(Creature<Chromosome> const &);

    Creature(Creature<Chromosome> &&);

    std::vector<std::pair<Chromosome, Chromosome>> const &getChromosomes()
        const;

    std::vector<std::pair<Chromosome, Chromosome>> &setChromosome();

    void initCreature();

    i32 getValue() const;

 private:
    std::vector<std::pair<Chromosome, Chromosome>> m_homolousChromosomes{};
};

template <typename Chromosome>
Creature<Chromosome>::Creature(usize chromosomeSize)
    : m_homolousChromosomes{chromosomeSize} {}

template <typename Chromosome>
Creature<Chromosome>::Creature(Creature<Chromosome> const &creature) {
    this->m_homolousChromosomes = creature.getChromosomes();
}

template <typename Chromosome>
Creature<Chromosome>::Creature(Creature<Chromosome> &&creature) {
    this->m_homolousChromosomes = std::move(creature.getChromosomes());
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

template <typename Chromosome>
void Creature<Chromosome>::initCreature() {
    this->m_homolousChromosomes = {{{1., 1}, {0., 0}}};
}

template <typename Chromosome>
i32 Creature<Chromosome>::getValue() const {
    return 1;
}

struct BasicData {
    static constexpr f64 variationRate{0.5};

    static constexpr f64 increasedCreatureRate{1.};
};

TEST(GeneticAlgorithmUnitTest, BirthNewCreauresUnitTest) {
    GeneticAlgorithm<Creature, Chromosome> testObject{
        MATH_CALCULATE::quickNatureNumberPow(10, 5)};
    testObject.initCreatures(&Creature<Chromosome>::initCreature);
    testObject.birthNewCreatures<BasicData>(
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
