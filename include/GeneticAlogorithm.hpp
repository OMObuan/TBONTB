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
class GeneticAlogorithm {
 public:
    GeneticAlogorithm() = default;

    explicit GeneticAlogorithm(usize initialCreatureNum);

    template <typename _FuncPointer>
        requires requires(Creature<Chromosome>* anyCreature,
                          _FuncPointer initChromosome) {
            (anyCreature->*initChromosome)();
        }
    void initCreatures(_FuncPointer);

    template <typename _ConstDate>
        requires requires() { _ConstDate::eliminatedRate; }
    static constexpr usize getEliminatedSize(usize);

    template <typename _ConstDate, typename _FuncPointer>
        requires requires(Creature<Chromosome> const& lftOrRgt,
                          _FuncPointer getValueFunction) {
            _ConstDate::eliminatedRate;
            (lftOrRgt.*getValueFunction)() < (lftOrRgt.*getValueFunction)();
        }
    void eliminateCreatures(_FuncPointer);

    template <typename _FuncPointer>
        requires requires(_FuncPointer const& getValueFunction,
                          Creature<Chromosome> const& anyCreature) {
            anyCreature.*getValueFunction() + anyCreature.getValueFunction();
        }
    auto getValueSum(_FuncPointer) const;

    template <typename _ConstDate, typename _FuncPointer>
        requires requires(const _ConstDate& constValues,
                          const _FuncPointer& getValueFunction,
                          const Creature<Chromosome>& anyCreature) {
            static_cast<double>(anyCreature.*getValueFunction() /
                                anyCreature.getValueFunction());
            constValues.probabilityPrecision;
        }
    std::vector<usize> generateRandArr(_ConstDate const&, _FuncPointer) const;

    Creature<Chromosome> getRandCreature(std::vector<usize> const&) const;

    static Chromosome getChormosomeFromBool(
        std::pair<Chromosome, Chromosome> const&, bool);

    template <typename _FuncPointer>
        requires requires(Chromosome const& anySide,
                          _FuncPointer getDominance) {
            anySide.*getDominance();
        }
    static std::pair<Chromosome, Chromosome> mixChromosome(Chromosome const&,
                                                           Chromosome const&,
                                                           _FuncPointer);

    template <typename _ConstDate, typename _FuncPointer>
    static void variateChromosome(
        _ConstDate const& constValues,
        std::pair<Chromosome, Chromosome> ChromosomePair,
        _FuncPointer variateFunction);

    template <typename _ConstDate, typename _FuncPointer,
              typename _FuncPointer2>
        requires(requires(Creature<Chromosome> anyCreature,
                          _FuncPointer getChormosome, Chromosome const& anySide,
                          _FuncPointer2 getDominance) {
            anyCreature.*getChormosome();
            anySide.*getDominance();
        } && std::is_same_v<decltype(std::declval<Creature<Chromosome>>().*
                                     std::declval<_FuncPointer>()),
                            std::vector<std::pair<Chromosome, Chromosome>>>)
    static Creature<Chromosome> fuckCreatrue(_ConstDate const&,
                                             const Creature<Chromosome>&,
                                             const Creature<Chromosome>&,
                                             _FuncPointer, _FuncPointer2);

    template <typename _ConstDate, typename _FuncPointer,
              typename _FuncPointer2>
        requires requires(_ConstDate const& constValues,
                          _FuncPointer getValueFunction,
                          Creature<Chromosome> const& anyCreature,
                          _FuncPointer2 fuckCreatrue) {
            constValues.probabilityPrecision;
            static_cast<double>(anyCreature.*getValueFunction() /
                                anyCreature.getValueFunction());
            anyCreature.*fuckCreatrue(anyCreature);
        }
    void birthNewCreature(_ConstDate const&, _FuncPointer, _FuncPointer2);

    std::vector<Creature<Chromosome>> getCreatures() const {
        return m_creatures;
    }

 private:
    std::vector<Creature<Chromosome>> m_creatures;
    static std::default_random_engine m_randomEngine;
};

template <template <typename Chromosome> class Creature, typename Chromosome>
GeneticAlogorithm<Creature, Chromosome>::GeneticAlogorithm(
    usize initialCreatureNum)
    : m_creatures{initialCreatureNum} {}

template <template <typename Chromosome> class Creature, typename Chromosome>
std::default_random_engine
    GeneticAlogorithm<Creature, Chromosome>::m_randomEngine{};

template <template <typename Chromosome> class Creature, typename Chromosome>
template <typename _FuncPointer>
    requires requires(Creature<Chromosome>* anyCreature,
                      _FuncPointer initChromosome) {
        (anyCreature->*initChromosome)();
    }
void GeneticAlogorithm<Creature, Chromosome>::initCreatures(
    _FuncPointer initChromosome) {
    for (auto& anyCreature : this->m_creatures) {
        (anyCreature.*initChromosome)();
    }
}

template <template <typename Chromosome> class Creature, typename Chromosome>
template <typename _ConstDate>
    requires requires() { _ConstDate::eliminatedRate; }
constexpr usize GeneticAlogorithm<Creature, Chromosome>::getEliminatedSize(
    usize size) {
    return static_cast<usize>(
        std::ceil(size * (1. - _ConstDate::eliminatedRate)));
}

template <template <typename Chromosome> class Creature, typename Chromosome>
template <typename _ConstDate, typename _FuncPointer>
    requires requires(Creature<Chromosome> const& lftOrRgt,
                      _FuncPointer getValueFunction) {
        _ConstDate::eliminatedRate;
        (lftOrRgt.*getValueFunction)() < (lftOrRgt.*getValueFunction)();
    }
void GeneticAlogorithm<Creature, Chromosome>::eliminateCreatures(
    _FuncPointer getValueFunction) {
    if (MATH_CALCULATE::doubleEqual(_ConstDate::eliminatedRate, 1., 1e-15)) {
        this->m_creatures.clear();
        return void();
    }
    usize dividePos =
        this->getEliminatedSize<_ConstDate>(this->m_creatures.size());
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
template <typename _FuncPointer>
    requires requires(_FuncPointer const& getValueFunction,
                      Creature<Chromosome> const& anyCreature) {
        anyCreature.*getValueFunction() + anyCreature.getValueFunction();
    }
auto GeneticAlogorithm<Creature, Chromosome>::getValueSum(
    _FuncPointer getValueFunction) const {
    decltype(std::declval<Creature<Chromosome>>().*getValueFunction()) sum{};
    for (auto& anyCreature : this->m_creatures) {
        sum = sum + anyCreature.*getValueFunction();
    }
    return sum;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
template <typename _ConstDate, typename _FuncPointer>
    requires requires(const _ConstDate& constValues,
                      const _FuncPointer& getValueFunction,
                      const Creature<Chromosome>& anyCreature) {
        static_cast<double>(anyCreature.*getValueFunction() /
                            anyCreature.getValueFunction());
        constValues.probabilityPrecision;
    }
std::vector<usize> GeneticAlogorithm<Creature, Chromosome>::generateRandArr(
    _ConstDate const& constValues, _FuncPointer getValueFunction) const {
    std::vector<usize> result{};
    usize randSize{MATH_CALCULATE::quickNatureNumberPow(
        10, constValues.probabilityPrecision)};
    auto valueSum{getValueSum(getValueFunction)};
    usize nowID{};
    for (auto& anyCreature : this->m_creatures) {
        double bit{static_cast<double>(anyCreature.*getValueFunction()) /
                   static_cast<double>(valueSum)};
        usize creatrueSize{static_cast<usize>(std::ceil(bit * randSize))};
        for (usize _i{1}; _i <= creatrueSize; ++_i) {
            result.push_back(nowID);
        }
        ++nowID;
    }
    return result;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
Creature<Chromosome> GeneticAlogorithm<Creature, Chromosome>::getRandCreature(
    std::vector<usize> const& randPos) const {
    std::uniform_int_distribution<usize> randGenerator(0, randPos.size() - 1);
    return randPos[randGenerator(this->m_randomEngine)];
}

template <template <typename Chromosome> class Creature, typename Chromosome>
Chromosome GeneticAlogorithm<Creature, Chromosome>::getChormosomeFromBool(
    std::pair<Chromosome, Chromosome> const& homologousChromosome,
    bool gamete) {
    return gamete ? homologousChromosome.first : homologousChromosome.second;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
template <typename _FuncPointer>
    requires requires(Chromosome const& anySide, _FuncPointer getDominance) {
        anySide.*getDominance();
    }
std::pair<Chromosome, Chromosome>
GeneticAlogorithm<Creature, Chromosome>::mixChromosome(
    Chromosome const& fatherSide, Chromosome const& motherSide,
    _FuncPointer getDominance) {
    if (fatherSide.*getDominance() >= motherSide.*getDominance()) {
        return std::make_pair(fatherSide, motherSide);
    } else {
        return std::make_pair(motherSide, fatherSide);
    }
}

template <template <typename Chromosome> class Creature, typename Chromosome>
template <typename _ConstDate, typename _FuncPointer>
void GeneticAlogorithm<Creature, Chromosome>::variateChromosome(
    _ConstDate const& constValues,
    std::pair<Chromosome, Chromosome> ChromosomePair,
    _FuncPointer variateFunction) {
    ChromosomePair.first.*variateFunction();
    ChromosomePair.second.*variateFunction();
}

template <template <typename Chromosome> class Creature, typename Chromosome>
template <typename _ConstDate, typename _FuncPointer, typename _FuncPointer2>
    requires(requires(Creature<Chromosome> anyCreature,
                      _FuncPointer getChormosome, Chromosome const& anySide,
                      _FuncPointer2 getDominance) {
        anyCreature.*getChormosome();
        anySide.*getDominance();
    } && std::is_same_v<decltype(std::declval<Creature<Chromosome>>().*
                                 std::declval<_FuncPointer>()),
                        std::vector<std::pair<Chromosome, Chromosome>>>)
Creature<Chromosome> GeneticAlogorithm<Creature, Chromosome>::fuckCreatrue(
    _ConstDate const& constValues, const Creature<Chromosome>& father,
    const Creature<Chromosome>& mother, _FuncPointer getChormosome,
    _FuncPointer2 getDominance) {
    std::uniform_int_distribution<usize> randChromosomeGenerator{0, 1};
    usize chromosomeSize = father.*getChormosome().size();
    for (usize _i{}; _i < chromosomeSize; ++_i) {
        bool fatherChromosome{randChromosomeGenerator(m_randomEngine)},
            motherChromosome{randChromosomeGenerator(m_randomEngine)};
        auto fatherSideChromosome{
            getChormosomeFromBool(father, fatherChromosome)},
            motherSideChromosome{
                getChormosomeFromBool(mother, motherChromosome)};
        variateChromosome(constValues, fatherSideChromosome,
                          motherSideChromosome);
        mixChromosome(fatherSideChromosome, motherSideChromosome, getDominance);
    }
}

template <template <typename Chromosome> class Creature, typename Chromosome>
template <typename _ConstDate, typename _FuncPointer, typename _FuncPointer2>
    requires requires(_ConstDate const& constValues,
                      _FuncPointer getValueFunction,
                      Creature<Chromosome> const& anyCreature,
                      _FuncPointer2 fuckCreatrue) {
        constValues.probabilityPrecision;
        static_cast<double>(anyCreature.*getValueFunction() /
                            anyCreature.getValueFunction());
        anyCreature.*fuckCreatrue(anyCreature);
    }
void GeneticAlogorithm<Creature, Chromosome>::birthNewCreature(
    _ConstDate const& constValues, _FuncPointer getValueFunction,
    _FuncPointer2 fuckCreatrue) {
    std::vector<usize> randPos{
        generateRandArr(constValues.probabilityPrecision, getValueFunction)};
    usize addCreatureSize =
        std::ceil(this->m_creatures.size() * constValues.increasedCreatureSize);
    for (usize _i = 1; _i <= addCreatureSize; ++_i) {
        Creature<Chromosome> father{getRandCreature(randPos)},
            mother{getRandCreature(randPos)},
            son = father.*fuckCreatrue(mother);
        m_creatures.push_back(son);
    }
}
