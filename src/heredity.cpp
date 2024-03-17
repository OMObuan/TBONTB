// Copyright [2024] <OMObuan>
#include <algorithm>
#include <cmath>
#include <rename_type.hpp>
#include <type_traits>
#include <utility>
#include <vector>

template <template <typename Chromosome> class Creature, typename Chromosome>
class GeneticAlogorithm {
    template <typename _T>
        requires requires(_T initChromosome, Creature<Chromosome> anyCreature) {
            anyCreature.*initChromosome();
        }
    void initCreatures(_T initChromosome) {
        for (auto& anyCreature : this->m_creatures) {
            anyCreature.*initChromosome();
        }
    }

    static usize getSize(usize size, f64 elimateRate) {
        return static_cast<usize>(std::ceil(size * elimateRate));
    }

    template <typename _T>
        requires requires(Creature<Chromosome> lftOrRgt, _T getValueFunction) {
            lftOrRgt.*getValueFunction() < lftOrRgt.*getValueFunction();
        }
    void eliminateCreatures(f64 eliminateRate, _T getValueFunction) {
        usize dividePos = getSize(this->m_creatures.size(), eliminateRate) - 1;
        std::ranges::nth_element(
            this->m_creatures, dividePos,
            [getValueFunction](Creature<Chromosome>& lft,
                               Creature<Chromosome>& rgt) -> bool {
                return lft.*getValueFunction() < rgt.*getValueFunction();
            });
        std::erase(this->m_creatures.begin() + dividePos,
                   this->m_creatures.end());
    }

    template <typename _T1, typename _T2,
              typename _Result = decltype(true ? std::declval<_T1>()
                                               : std::declval<_T2>())>
        requires((std::is_same_v<_T2, i32> || std::is_same_v<_T2, usize>) &&
                 requires(_Result result) { result* result; })
    auto quickNatureNumberPow(const _T1& baseNumber, const _T2& indexNumber) {
        _Result result{1}, nowBaseNumber{baseNumber};
        while (indexNumber) {
            if (indexNumber % 2 == 0) {
                nowBaseNumber = nowBaseNumber * nowBaseNumber;
                indexNumber /= 2;
            } else {
                result *= nowBaseNumber;
                --indexNumber;
            }
        }
        return result;
    }

    void birthNewCreature(usize probabilityPrecision) {
        std::vector<usize> randPos{};
        usize randSize = quickNatureNumberPow(10, probabilityPrecision);
        for (auto& anyCreature : this->m_creatures) {
        }
    }

 private:
    std::vector<Creature<Chromosome>> m_creatures;
};
