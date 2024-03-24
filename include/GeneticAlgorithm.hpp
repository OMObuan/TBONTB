// Copyright [2024] <OMObuan>

#pragma once
#include <gtest/gtest.h>
#include <math_calculate.h>
#include <rename_type.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
class GeneticAlgorithm {
 public:
    GeneticAlgorithm() = default;

    explicit GeneticAlgorithm(usize initialCreatureNum);

    template <typename _FuncPointer>
        requires requires(Creature<Chromosome>* anyCreature,
                          _FuncPointer initChromosome) {
            (anyCreature->*initChromosome)();
        }
    void initCreatures(_FuncPointer);

    template <typename _BasicDate>
        requires requires() { _BasicDate::eliminatedRate; }
    static constexpr usize getEliminatedSize(usize);

    template <typename _BasicDate, typename _FuncPointer>
        requires requires(Creature<Chromosome> const& lftOrRgt,
                          _FuncPointer getValueFunction) {
            _BasicDate::eliminatedRate;
            (lftOrRgt.*getValueFunction)() < (lftOrRgt.*getValueFunction)();
        }
    void eliminateCreatures(_FuncPointer);

    template <typename _FuncPointer>
        requires requires(_FuncPointer getValueFunction,
                          Creature<Chromosome>* anyCreature) {
            (anyCreature->*getValueFunction)() +
                (anyCreature->*getValueFunction)();
        }
    auto getValueSum(_FuncPointer) const;

    template <typename _BasicDate, typename _FuncPointer>
        requires requires(_FuncPointer getValueFunction,
                          const Creature<Chromosome>& anyCreature) {
            static_cast<f64>((anyCreature.*getValueFunction)() /
                             (anyCreature.*getValueFunction)());
            _BasicDate::probabilityPrecision;
        }
    std::vector<usize> generateRandArray(_FuncPointer) const;

    Creature<Chromosome> getRandCreature(std::vector<usize> const&) const;

    static Chromosome getChormosomeFromBool(
        std::pair<Chromosome, Chromosome> const&, bool);

    template <typename _FuncPointer>
        requires requires(Chromosome const& anySide,
                          _FuncPointer getDominance) {
            (anySide.*getDominance)();
        }
    static std::pair<Chromosome, Chromosome> mixChromosome(Chromosome const&,
                                                           Chromosome const&,
                                                           _FuncPointer);

    template <typename _FuncPointer>
        requires requires(Chromosome* variatedChromosome,
                          _FuncPointer variateFunction) {
            (variatedChromosome->*variateFunction)();
        }
    static void variateChromosome(Chromosome*, _FuncPointer);

    template <typename _BasicDate, typename _FuncPointer,
              typename _FuncPointer2, typename _FuncPointer3,
              typename _FuncPointer4>
        requires(requires(Creature<Chromosome> const& c_anyCreature,
                          Creature<Chromosome>* anyCreature,
                          _FuncPointer getChormosome,
                          _FuncPointer2 setChromosome,
                          Chromosome const& c_anySide, Chromosome* anySide,
                          _FuncPointer3 getDominance,
                          _FuncPointer4 variateFunction) {
            (anyCreature->*getChormosome)();
            (anyCreature->*setChromosome)(0,
                                          std::make_pair(c_anySide, c_anySide));
            (c_anySide.*getDominance)();
            (anySide->*variateFunction)();
        } && std::is_same_v<decltype((std::declval<Creature<Chromosome>>().*
                                      std::declval<_FuncPointer>)()),
                            std::vector<std::pair<Chromosome, Chromosome>>>)
    Creature<Chromosome> mateCreature(const Creature<Chromosome>& father,
                                      const Creature<Chromosome>& mother,
                                      _FuncPointer getChormosome,
                                      _FuncPointer2 setChromosome,
                                      _FuncPointer3 getDominance,
                                      _FuncPointer4 variateFunction);

    template <typename _BasicDate, typename _FuncPointer>
        requires requires(_FuncPointer getValueFunction,
                          Creature<Chromosome> const& anyCreature) {
            _BasicDate::probabilityPrecision;
            static_cast<f64>((anyCreature.*getValueFunction)() /
                             (anyCreature.*getValueFunction)());
            anyCreature.*fuckCreature(anyCreature);
        }
    void birthNewCreature(_FuncPointer);

    std::vector<Creature<Chromosome>> const& getCreatures() const;

 private:
    std::vector<Creature<Chromosome>> m_creatures{};
    static std::default_random_engine m_randomEngine;
};

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
std::default_random_engine
    GeneticAlgorithm<Creature, Chromosome>::m_randomEngine{
        static_cast<usize>(std::time(nullptr))};

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
GeneticAlgorithm<Creature, Chromosome>::GeneticAlgorithm(
    usize initialCreatureNum)
    : m_creatures{initialCreatureNum} {}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _FuncPointer>
    requires requires(Creature<Chromosome>* anyCreature,
                      _FuncPointer initChromosome) {
        (anyCreature->*initChromosome)();
    }
void GeneticAlgorithm<Creature, Chromosome>::initCreatures(
    _FuncPointer initChromosome) {
    for (auto& anyCreature : this->m_creatures) {
        (anyCreature.*initChromosome)();
    }
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicDate>
    requires requires() { _BasicDate::eliminatedRate; }
constexpr usize GeneticAlgorithm<Creature, Chromosome>::getEliminatedSize(
    usize size) {
    return static_cast<usize>(
        std::ceil(size * (1. - _BasicDate::eliminatedRate)));
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicDate, typename _FuncPointer>
    requires requires(Creature<Chromosome> const& lftOrRgt,
                      _FuncPointer getValueFunction) {
        _BasicDate::eliminatedRate;
        (lftOrRgt.*getValueFunction)() < (lftOrRgt.*getValueFunction)();
    }
void GeneticAlgorithm<Creature, Chromosome>::eliminateCreatures(
    _FuncPointer getValueFunction) {
    if (MATH_CALCULATE::f64Equal(_BasicDate::eliminatedRate, 1., 1e-15)) {
        this->m_creatures.clear();
        return void();
    }
    usize dividePos =
        this->getEliminatedSize<_BasicDate>(this->m_creatures.size());
    std::ranges::nth_element(
        this->m_creatures, this->m_creatures.begin() + dividePos,
        [getValueFunction](Creature<Chromosome>& lft,
                           Creature<Chromosome>& rgt) -> bool {
            return (rgt.*getValueFunction)() < (lft.*getValueFunction)();
        });
    // std::iter_swap(this->m_creatures.begin(),
    //                this->m_creatures.begin() + dividePos );
    this->m_creatures.erase(this->m_creatures.begin() + dividePos,
                            this->m_creatures.end());
    std::reverse(this->m_creatures.begin(), this->m_creatures.end());
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _FuncPointer>
    requires requires(_FuncPointer getValueFunction,
                      Creature<Chromosome>* anyCreature) {
        (anyCreature->*getValueFunction)() + (anyCreature->*getValueFunction)();
    }
auto GeneticAlgorithm<Creature, Chromosome>::getValueSum(
    _FuncPointer getValueFunction) const {
    decltype((std::declval<Creature<Chromosome>>().*getValueFunction)()) sum{};
    for (auto& anyCreature : this->m_creatures) {
        sum = sum + (anyCreature.*getValueFunction)();
    }
    return sum;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicDate, typename _FuncPointer>
    requires requires(_FuncPointer getValueFunction,
                      const Creature<Chromosome>& anyCreature) {
        static_cast<f64>((anyCreature.*getValueFunction)() /
                         (anyCreature.*getValueFunction)());
        _BasicDate::probabilityPrecision;
    }
std::vector<usize> GeneticAlgorithm<Creature, Chromosome>::generateRandArray(
    _FuncPointer getValueFunction) const {
    std::vector<usize> result{};
    usize randSize{MATH_CALCULATE::quickNatureNumberPow(
        10, _BasicDate::probabilityPrecision)};
    auto valueSum{getValueSum(getValueFunction)};
    if (valueSum == 0) {
        return result;
    }
    usize nowID{};
    for (auto& anyCreature : this->m_creatures) {
        f64 bit{static_cast<f64>((anyCreature.*getValueFunction)()) /
                static_cast<f64>(valueSum)};
        // std::cerr << bit << '\n';
        usize creatrueSize{static_cast<usize>(std::ceil(bit * randSize))};
        // std::cerr << creatrueSize << '\n';
        for (usize _i{1}; _i <= creatrueSize; ++_i) {
            result.push_back(nowID);
        }
        ++nowID;
    }
    return result;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
Creature<Chromosome> GeneticAlgorithm<Creature, Chromosome>::getRandCreature(
    std::vector<usize> const& randPos) const {
    std::uniform_int_distribution<usize> randGenerator(0, randPos.size() - 1);
    return this->m_creatures[randPos[randGenerator(this->m_randomEngine)]];
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
Chromosome GeneticAlgorithm<Creature, Chromosome>::getChormosomeFromBool(
    std::pair<Chromosome, Chromosome> const& homologousChromosome,
    bool gamete) {
    return gamete ? homologousChromosome.first : homologousChromosome.second;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _FuncPointer>
    requires requires(Chromosome const& anySide, _FuncPointer getDominance) {
        (anySide.*getDominance)();
    }
std::pair<Chromosome, Chromosome>
GeneticAlgorithm<Creature, Chromosome>::mixChromosome(
    Chromosome const& fatherSide, Chromosome const& motherSide,
    _FuncPointer getDominance) {
    if ((fatherSide.*getDominance)() >= (motherSide.*getDominance)()) {
        return std::make_pair(fatherSide, motherSide);
    } else {
        return std::make_pair(motherSide, fatherSide);
    }
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _FuncPointer>
    requires requires(Chromosome* variatedChromosome,
                      _FuncPointer variateFunction) {
        (variatedChromosome->*variateFunction)();
    }
void GeneticAlgorithm<Creature, Chromosome>::variateChromosome(
    Chromosome* variatedChromosome, _FuncPointer variateChromosome) {
    (variatedChromosome->*variateChromosome)();
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicDate, typename _FuncPointer, typename _FuncPointer2,
          typename _FuncPointer3, typename _FuncPointer4>
    requires(requires(Creature<Chromosome> const& c_anyCreature,
                      Creature<Chromosome>* anyCreature,
                      _FuncPointer getChormosome, _FuncPointer2 setChromosome,
                      Chromosome const& c_anySide, Chromosome* anySide,
                      _FuncPointer3 getDominance,
                      _FuncPointer4 variateFunction) {
        (anyCreature->*getChormosome)();
        (anyCreature->*setChromosome)(0, std::make_pair(c_anySide, c_anySide));
        (c_anySide.*getDominance)();
        (anySide->*variateFunction)();
    } && std::is_same_v<decltype((std::declval<Creature<Chromosome>>().*
                                  std::declval<_FuncPointer>)()),
                        std::vector<std::pair<Chromosome, Chromosome>>>)
Creature<Chromosome> GeneticAlgorithm<Creature, Chromosome>::mateCreature(
    const Creature<Chromosome>& father, const Creature<Chromosome>& mother,
    _FuncPointer getChormosome, _FuncPointer2 setChromosome,
    _FuncPointer3 getDominance, _FuncPointer4 variateFunction) {
    std::uniform_int_distribution<usize> randChromosomeGenerator{0, 1};
    std::uniform_real_distribution<f64> randVariationGenerator{0, 1};
    usize chromosomeSize = (father.*getChormosome()).size();
    Creature<Chromosome> result{chromosomeSize};
    for (usize _i{}; _i < chromosomeSize; ++_i) {
        bool fatherChromosome{randChromosomeGenerator(m_randomEngine)},
            motherChromosome{randChromosomeGenerator(m_randomEngine)};
        auto fatherSideChromosome{getChormosomeFromBool(
            father.getChormosome()[_i], fatherChromosome)},
            motherSideChromosome{getChormosomeFromBool(
                mother.getChormosome()[_i], motherChromosome)};
        if (randVariationGenerator(m_randomEngine) <
            _BasicDate::variationRate) {
            variateChromosome(&fatherSideChromosome, variateFunction);
        }
        if (randVariationGenerator(m_randomEngine) <
            _BasicDate::variationRate) {
            variateChromosome(&motherSideChromosome, variateFunction);
        }
        (result.*setChromosome)(
            _i, mixChromosome(fatherSideChromosome, motherSideChromosome,
                              getDominance));
    }
    return result;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicDate, typename _FuncPointer>
    requires requires(_FuncPointer getValueFunction,
                      Creature<Chromosome> const& anyCreature) {
        _BasicDate::probabilityPrecision;
        static_cast<f64>((anyCreature.*getValueFunction)() /
                         (anyCreature.*getValueFunction)());
        anyCreature.*fuckCreature(anyCreature);
    }
void GeneticAlgorithm<Creature, Chromosome>::birthNewCreature(
    _FuncPointer getValueFunction) {
    std::vector<usize> randPos{
        generateRandArray(_BasicDate::probabilityPrecision, getValueFunction)};
    usize addCreatureSize =
        std::ceil(this->m_creatures.size() * _BasicDate::increasedCreatureSize);
    for (usize _i = 1; _i <= addCreatureSize; ++_i) {
        Creature<Chromosome> father{getRandCreature(randPos)},
            mother{getRandCreature(randPos)},
            son{GeneticAlgorithm::mateCreature(father, mother,
                                               getValueFunction)};
        m_creatures.push_back(son);
    }
}
template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
std::vector<Creature<Chromosome>> const&
GeneticAlgorithm<Creature, Chromosome>::getCreatures() const {
    return m_creatures;
}
