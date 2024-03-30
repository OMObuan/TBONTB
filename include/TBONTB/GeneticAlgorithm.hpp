// Copyright [2024] <OMObuan>

#pragma once
#include <TBONTB/debug.h>
#include <TBONTB/math_calculate.h>
#include <TBONTB/rename_type.h>

#include <algorithm>
#include <cmath>
#include <ctime>
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

    template <typename _FuncPointer>
        requires requires(Creature<Chromosome> const& anyCreature,
                          _FuncPointer printChromosomes) {
            (anyCreature.*printChromosomes)();
        }
    void printCreatures(_FuncPointer) const;

    template <typename _BasicData>
        requires requires() { _BasicData::eliminatedRate; }
    static constexpr usize getEliminatedSize(usize);

    template <typename _BasicData, typename _FuncPointer>
        requires requires(Creature<Chromosome> const& lftOrRgt,
                          _FuncPointer getValueFunction) {
            _BasicData::eliminatedRate;
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

    template <typename _FuncPointer>
        requires requires(_FuncPointer getValueFunction,
                          Creature<Chromosome> const& anyCreature) {
            static_cast<f64>((anyCreature.*getValueFunction)() /
                             (anyCreature.*getValueFunction)());
        }
    std::discrete_distribution<> generateRandArray(_FuncPointer) const;

    Creature<Chromosome> getRandCreature(std::discrete_distribution<>*);

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

    template <typename _BasicData, typename _FuncPointer,
              typename _FuncPointer2, typename _FuncPointer3,
              typename _FuncPointer4>
        requires(
            requires(Creature<Chromosome> const& c_anyCreature,
                     Creature<Chromosome>* anyCreature,
                     _FuncPointer getChormosome, _FuncPointer2 setChromosome,
                     Chromosome const& c_anySide, Chromosome* anySide,
                     _FuncPointer3 getDominance,
                     _FuncPointer4 variateFunction) {
                (anyCreature->*getChormosome)();
                (c_anySide.*getDominance)();
                (anySide->*variateFunction)();
            } &&
            std::is_same_v<std::decay_t<decltype(((Creature<Chromosome>{}).*
                                                  (_FuncPointer{}))())>,
                           std::vector<std::pair<Chromosome, Chromosome>>> &&
            std::is_same_v<
                std::decay_t<decltype(((std::decay_t<Creature<Chromosome>>()).*
                                       (_FuncPointer2{}))())>,
                std::vector<std::pair<Chromosome, Chromosome>>>)
    static Creature<Chromosome> mateCreature(const Creature<Chromosome>&,
                                             const Creature<Chromosome>&,
                                             _FuncPointer, _FuncPointer2,
                                             _FuncPointer3, _FuncPointer4);

    template <typename _BasicData, typename _FuncPointer,
              typename _FuncPointer2, typename _FuncPointer3,
              typename _FuncPointer4, typename _FuncPointer5>
        requires(
            requires(Creature<Chromosome> const& c_anyCreature,
                     Creature<Chromosome>* anyCreature,
                     _FuncPointer getChormosome, _FuncPointer2 setChromosome,
                     Chromosome const& c_anySide, Chromosome* anySide,
                     _FuncPointer3 getDominance,
                     _FuncPointer4 variateFunction) {
                (anyCreature->*getChormosome)();
                (c_anySide.*getDominance)();
                (anySide->*variateFunction)();
                _BasicData::increasedCreatureRate;
            } &&
            std::is_same_v<std::decay_t<decltype(((Creature<Chromosome>{}).*
                                                  (_FuncPointer{}))())>,
                           std::vector<std::pair<Chromosome, Chromosome>>> &&
            std::is_same_v<
                std::decay_t<decltype(((std::decay_t<Creature<Chromosome>>()).*
                                       (_FuncPointer2{}))())>,
                std::vector<std::pair<Chromosome, Chromosome>>> &&
            requires(_FuncPointer5 getValueFunction,
                     const Creature<Chromosome>& anyCreature) {
                static_cast<f64>((anyCreature.*getValueFunction)() /
                                 (anyCreature.*getValueFunction)());
            })
    void birthNewCreatures(_FuncPointer, _FuncPointer2, _FuncPointer3,
                           _FuncPointer4, _FuncPointer5);

    std::vector<Creature<Chromosome>> const& getPopulation() const;

    std::vector<Creature<Chromosome>>& setPopulation();

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
template <typename _FuncPointer>
    requires requires(Creature<Chromosome> const& anyCreature,
                      _FuncPointer printChromosomes) {
        (anyCreature.*printChromosomes)();
    }
void GeneticAlgorithm<Creature, Chromosome>::printCreatures(
    _FuncPointer printChromosome) const {
    for (auto& anyCreature : this->m_creatures) {
        (anyCreature.*printChromosome)();
    }
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicData>
    requires requires() { _BasicData::eliminatedRate; }
constexpr usize GeneticAlgorithm<Creature, Chromosome>::getEliminatedSize(
    usize size) {
    return static_cast<usize>(
        std::ceil(size * (1. - _BasicData::eliminatedRate)));
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicData, typename _FuncPointer>
    requires requires(Creature<Chromosome> const& lftOrRgt,
                      _FuncPointer getValueFunction) {
        _BasicData::eliminatedRate;
        (lftOrRgt.*getValueFunction)() < (lftOrRgt.*getValueFunction)();
    }
void GeneticAlgorithm<Creature, Chromosome>::eliminateCreatures(
    _FuncPointer getValueFunction) {
    if (MATH_CALCULATE::f64Equal(_BasicData::eliminatedRate, 1., 1e-15)) {
        this->m_creatures.clear();
        return void();
    }
    usize dividePos =
        this->getEliminatedSize<_BasicData>(this->m_creatures.size());
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
template <typename _FuncPointer>
    requires requires(_FuncPointer getValueFunction,
                      Creature<Chromosome> const& anyCreature) {
        static_cast<f64>((anyCreature.*getValueFunction)() /
                         (anyCreature.*getValueFunction)());
    }
std::discrete_distribution<>
GeneticAlgorithm<Creature, Chromosome>::generateRandArray(
    _FuncPointer getValueFunction) const {
    std::vector<f64> everyCreautreValue;
    for (auto& anyCreature : this->m_creatures) {
        everyCreautreValue.push_back((anyCreature.*getValueFunction)());
    }
    return std::discrete_distribution<>(everyCreautreValue.begin(),
                                        everyCreautreValue.end());
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
Creature<Chromosome> GeneticAlgorithm<Creature, Chromosome>::getRandCreature(
    std::discrete_distribution<>* randPos) {
    return this->m_creatures[static_cast<usize>(
        (*randPos)(GeneticAlgorithm::m_randomEngine))];
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
    Chromosome* variatedChromosome, _FuncPointer variateOperate) {
    (variatedChromosome->*variateOperate)();
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicData, typename _FuncPointer, typename _FuncPointer2,
          typename _FuncPointer3, typename _FuncPointer4>
    requires(requires(Creature<Chromosome> const& c_anyCreature,
                      Creature<Chromosome>* anyCreature,
                      _FuncPointer getChormosome, _FuncPointer2 setChromosome,
                      Chromosome const& c_anySide, Chromosome* anySide,
                      _FuncPointer3 getDominance,
                      _FuncPointer4 variateFunction) {
        (anyCreature->*getChormosome)();
        (c_anySide.*getDominance)();
        (anySide->*variateFunction)();
    } &&
             std::is_same_v<std::decay_t<decltype(((Creature<Chromosome>{}).*
                                                   (_FuncPointer{}))())>,
                            std::vector<std::pair<Chromosome, Chromosome>>> &&
             std::is_same_v<
                 std::decay_t<decltype(((std::decay_t<Creature<Chromosome>>()).*
                                        (_FuncPointer2{}))())>,
                 std::vector<std::pair<Chromosome, Chromosome>>>)
Creature<Chromosome> GeneticAlgorithm<Creature, Chromosome>::mateCreature(
    Creature<Chromosome> const& father, Creature<Chromosome> const& mother,
    _FuncPointer getChormosomes, _FuncPointer2 setChromosomes,
    _FuncPointer3 getDominance, _FuncPointer4 variateOperate) {
    std::uniform_int_distribution<usize> randChromosomeGenerator{0, 1};
    std::uniform_real_distribution<f64> randVariationGenerator{0, 1};
    usize chromosomeSize = (father.*getChormosomes)().size();
    // _LOOK(std::cerr, chromosomeSize);
    Creature<Chromosome> result{chromosomeSize};
    for (usize _i{}; _i < chromosomeSize; ++_i) {
        bool fatherChromosome{randChromosomeGenerator(m_randomEngine)},
            motherChromosome{randChromosomeGenerator(m_randomEngine)};
        auto fatherSideChromosome{getChormosomeFromBool(
            (father.*getChormosomes)()[_i], fatherChromosome)},
            motherSideChromosome{getChormosomeFromBool(
                (mother.*getChormosomes)()[_i], motherChromosome)};
        if (randVariationGenerator(m_randomEngine) <
            _BasicData::variationRate) {
            variateChromosome(&fatherSideChromosome, variateOperate);
        }
        if (randVariationGenerator(m_randomEngine) <
            _BasicData::variationRate) {
            variateChromosome(&motherSideChromosome, variateOperate);
        }
        (result.*setChromosomes)()[_i] = mixChromosome(
            fatherSideChromosome, motherSideChromosome, getDominance);
    }
    // _LOOK(std::cerr, (result.*getChormosomes)().size());
    return result;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
template <typename _BasicData, typename _FuncPointer, typename _FuncPointer2,
          typename _FuncPointer3, typename _FuncPointer4,
          typename _FuncPointer5>
    requires(
        requires(Creature<Chromosome> const& c_anyCreature,
                 Creature<Chromosome>* anyCreature, _FuncPointer getChormosome,
                 _FuncPointer2 setChromosome, Chromosome const& c_anySide,
                 Chromosome* anySide, _FuncPointer3 getDominance,
                 _FuncPointer4 variateFunction) {
            (anyCreature->*getChormosome)();
            (c_anySide.*getDominance)();
            (anySide->*variateFunction)();
            _BasicData::increasedCreatureRate;
        } &&
        std::is_same_v<std::decay_t<decltype(((Creature<Chromosome>{}).*
                                              (_FuncPointer{}))())>,
                       std::vector<std::pair<Chromosome, Chromosome>>> &&
        std::is_same_v<
            std::decay_t<decltype(((std::decay_t<Creature<Chromosome>>()).*
                                   (_FuncPointer2{}))())>,
            std::vector<std::pair<Chromosome, Chromosome>>> &&
        requires(_FuncPointer5 getValueFunction,
                 const Creature<Chromosome>& anyCreature) {
            static_cast<f64>((anyCreature.*getValueFunction)() /
                             (anyCreature.*getValueFunction)());
        })
void GeneticAlgorithm<Creature, Chromosome>::birthNewCreatures(
    _FuncPointer getChormosomes, _FuncPointer2 setChromosomes,
    _FuncPointer3 getDominance, _FuncPointer4 variateOperate,
    _FuncPointer5 getValueFunction) {
    // generateRandArray(getValueFunction);
    auto randPosGenerator{generateRandArray(getValueFunction)};
    usize addCreatureSize =
        std::ceil(this->m_creatures.size() * _BasicData::increasedCreatureRate);
    // _LOOK(std::cerr, addCreatureSize);
    for (usize _i = 1; _i <= addCreatureSize; ++_i) {
        // getRandCreature(randPosGenerator);
        Creature<Chromosome> father{getRandCreature(&randPosGenerator)},
            mother{getRandCreature(&randPosGenerator)},
            son{GeneticAlgorithm::mateCreature<_BasicData>(
                father, mother, getChormosomes, setChromosomes, getDominance,
                variateOperate)};
        this->m_creatures.push_back(son);
        // _LOOK(std::cerr, _i);
        // _LOOK(std::cerr, (_i <= addCreatureSize));
    }
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
std::vector<Creature<Chromosome>> const&
GeneticAlgorithm<Creature, Chromosome>::getPopulation() const {
    return m_creatures;
}

template <template <typename Chromosome> class Creature, typename Chromosome>
    requires requires() { Creature<Chromosome>{usize{}}; }
std::vector<Creature<Chromosome>>&
GeneticAlgorithm<Creature, Chromosome>::setPopulation() {
    return this->m_creatures;
}
