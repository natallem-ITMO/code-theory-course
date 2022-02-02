//
// Created by natallem on 10/10/21.
//

#ifndef POLARCODING_POLARDECODER_HPP
#define POLARCODING_POLARDECODER_HPP

#include <vector>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <cstdint>

namespace util{
    template<typename T>
    static T getMedian(std::vector<T> &A) {
        std::sort(A.begin(), A.end());
        int median_index = (A.size() - 1) / 2;
        T ret = A[median_index];
        A.clear();
        return ret;
    }

    template<typename T>
    static int partition(std::vector<T> &A, int l, int r, T x) {
        int i;
        for (i = l; i < r; i++)
            if (A[i] == x)
                break;
        std::swap(A[i], A[r]);

        i = l;
        for (int j = l; j <= r - 1; j++) {
            if (A[j] <= x) {
                std::swap(A[i], A[j]);
                i++;
            }
        }
        std::swap(A[i], A[r]);
        return i;
    }

    template<typename T>
    static T select(std::vector<T> &A, int beg, int end, int i_th_smallest) {
        assert(end > beg);
        std::vector<T> medians_5_groups, current_group;
        for (int index = beg; index < end; ++index) {
            current_group.push_back(A[index]);
            if (current_group.size() == 5) {
                medians_5_groups.push_back(getMedian(current_group));
            }
        }
        if (!current_group.empty()) {
            medians_5_groups.push_back(getMedian(current_group));
        }
        T x = (medians_5_groups.size() == 1)
              ? medians_5_groups[0]
              : select(medians_5_groups, 0, medians_5_groups.size(), medians_5_groups.size() / 2);
        auto x_index = partition(A, beg, end - 1, x);
        if (x_index - beg == i_th_smallest - 1) {
            return x;
        } else if (x_index - beg > i_th_smallest - 1) {
            return select(A, beg, x_index, i_th_smallest);
        }
        return select(A, x_index + 1, end, i_th_smallest - x_index + beg - 1);
    }
}

class PolarDecoder {
public:

    std::vector<uint8_t> SCL_decoder(int L_input, int m_input, std::vector<double> W_y_1, std::vector<int> &frozen) {
        L = L_input;
        m = m_input;
        frozen_bits_num = frozen;
        int n = (1 << m);
        initializeDataStructures();
        int l = assignInitialPath();
        auto &P_0 = getArrayPointer_P(0, l);
        for (int b = 0; b < n; ++b) {
            P_0[2 * b] = W_y_1[b] == 1 ? 0 : 1;
            P_0[2 * b + 1] = 1 - P_0[2 * b];
        }
        for (int phi = 0; phi < n; ++phi) {
            recursivelyCalcP(m, phi);
            if (isFrozen(phi)) {
                continuePaths_FrozenBit(phi);
            } else {
                continuePaths_UnfrozenBit(phi);
            }
            if (phi % 2 == 1) {
                recursivelyUpdateC(m, phi);
            }
        }
        l = findMostProbablePath();
        std::vector<uint8_t> &c_0 = arrayPointer_I[l];
        std::vector<uint8_t> decoded_info_bits((1 << m) - frozen_bits_num.size());
        for (int beta = 0; beta < decoded_info_bits.size(); ++beta)
            decoded_info_bits[beta] = c_0[beta];
        return decoded_info_bits;
    }

private:

    void initializeDataStructures() {
        inactivePathIndices.clear();
        inactivePathIndices.resize(L);
        std::iota(inactivePathIndices.begin(), inactivePathIndices.end(), 0);
        arrayPointer_P.resize(m + 1, std::vector<std::vector<double>>(L));
        arrayPointer_C.resize(m + 1, std::vector<std::vector<uint8_t>>(L));
        arrayPointer_I.resize(L, std::vector<uint8_t>(1 << m, 0));
        arrayReferenceCount.resize(m + 1, std::vector<int>(L, 0));
        activePath.resize(L, false);
        inactiveArrayIndices.clear();
        inactiveArrayIndices.resize(m + 1);
        pathIndexToArrayIndex.resize(m + 1, std::vector<int>(L));
        for (int l = 0; l <= m; ++l) {
            for (int s = 0; s < L; s++) {
                arrayPointer_P[l][s] = std::vector<double>(1 << (m - l + 1));
                arrayPointer_C[l][s] = std::vector<uint8_t>(1 << (m - l + 1));
                inactiveArrayIndices[l].push_back(s);
            }
        }
    }

    int assignInitialPath() {
        int l = inactivePathIndices.back();
        inactivePathIndices.pop_back();
        activePath[l] = true;
        for (int la = 0; la <= m; ++la) {
            int s = inactiveArrayIndices[la].back();
            inactiveArrayIndices[la].pop_back();
            pathIndexToArrayIndex[la][l] = s;
            arrayReferenceCount[la][s] = 1;
        }
        return l;
    }

    int clonePath(int l) {
        int l_ = inactivePathIndices.back();
        inactivePathIndices.pop_back();
        activePath[l_] = true;
        for (int la = 0; la <= m; ++la) {
            int s = pathIndexToArrayIndex[la][l];
            pathIndexToArrayIndex[la][l_] = s;
            ++arrayReferenceCount[la][s];
        }
        return l_;
    }

    void killPath(int l) {
        activePath[l] = false;
        inactivePathIndices.push_back(l);
        for (int la = 0; la <= m; ++la) {
            int s = pathIndexToArrayIndex[la][l];
            --arrayReferenceCount[la][s];
            if (arrayReferenceCount[la][s] == 0) {
                inactiveArrayIndices[la].push_back(s);
            }
        }
    }

    std::vector<double> &getArrayPointer_P(int la, int l) {
        int s = pathIndexToArrayIndex[la][l];
        int s_ = s;
        if (arrayReferenceCount[la][s] != 1) {
            s_ = inactiveArrayIndices[la].back();
            inactiveArrayIndices[la].pop_back();
            arrayPointer_C[la][s_] = arrayPointer_C[la][s];
            arrayPointer_P[la][s_] = arrayPointer_P[la][s];
            --arrayReferenceCount[la][s];
            arrayReferenceCount[la][s_] = 1;
            pathIndexToArrayIndex[la][l] = s_;
        }
        return arrayPointer_P[la][s_];
    }

    std::vector<uint8_t> &getArrayPointer_C(int la, int l) {
//        std::cout << "la=" << la << " l=" << l << "\n";
        int s = pathIndexToArrayIndex[la][l];
        int s_ = s;
        if (la == 3 && l == 4 && arrayReferenceCount[0][4] == 4) {
            int x = 10;
        }
        if (arrayReferenceCount[la][s] != 1) {
            assert(!inactiveArrayIndices[la].empty());
            s_ = inactiveArrayIndices[la].back();
            inactiveArrayIndices[la].pop_back();
            arrayPointer_C[la][s_] = arrayPointer_C[la][s];
            arrayPointer_P[la][s_] = arrayPointer_P[la][s];
            --arrayReferenceCount[la][s];
            arrayReferenceCount[la][s_] = 1;
            pathIndexToArrayIndex[la][l] = s_;
        }
        return arrayPointer_C[la][s_];
    }

    void recursivelyCalcP(int la, int phi) {
        if (la == 0) {
            return;
        }
        int psi = phi >> 1;
        if (phi % 2 == 0) {
            recursivelyCalcP(la - 1, psi);
        }
        double a = 0;
        for (int l = 0; l < L; ++l) {
            if (pathIndexInactive(l)) {
                continue;
            }
            std::vector<double> &P_la = getArrayPointer_P(la, l);
            std::vector<double> &P_la_minus_1 = getArrayPointer_P(la - 1, l);
            std::vector<uint8_t> &C_la = getArrayPointer_C(la, l);
            for (int b = 0; b < (1 << (m - la)); ++b) {
                if (phi % 2 == 0) {
                    P_la[2 * b] = 0.5 * (P_la_minus_1[2 * (2 * b)] * P_la_minus_1[2 * (2 * b + 1)]
                                          + P_la_minus_1[2 * (2 * b) + 1] * P_la_minus_1[2 * (2 * b + 1) + 1]);
                    P_la[2 * b + 1] = 0.5 * (P_la_minus_1[2 * (2 * b) + 1] * P_la_minus_1[2 * (2 * b + 1)]
                                              + P_la_minus_1[2 * (2 * b)] * P_la_minus_1[2 * (2 * b + 1) + 1]);
                } else {
                    uint8_t u = C_la[2 * b];
                    P_la[2 * b] = 0.5 * P_la_minus_1[2 * (2 * b) + (u % 2)] * P_la_minus_1[2 * (2 * b + 1)];
                    P_la[2 * b + 1] =
                            0.5 * P_la_minus_1[2 * (2 * b) + ((u + 1) % 2)] * P_la_minus_1[2 * (2 * b + 1) + 1];
                }
                a = std::max(a, std::max(P_la[2 * b], P_la[2 * b + 1]));
            }
        }
        for (int l = 0; l < L; ++l) {
            if (pathIndexInactive(l)) {
                continue;
            }
            auto &P_la = getArrayPointer_P(la, l);
            for (int b = 0; b < (1 << (m - la)); ++b) {
                P_la[2 * b] /= a;
                P_la[2 * b + 1] /= a;
            }
        }
    }

    bool pathIndexInactive(int l) {
        return !activePath[l];
    }

    void recursivelyUpdateC(int la, int phi) {
        assert(phi % 2 == 1);
        int psi = phi >> 1;
        for (int l = 0; l < L; ++l) {
            if (pathIndexInactive(l)) {
                continue;
            }
            auto &c_lambda = getArrayPointer_C(la, l);
            auto &c_lambda_1 = getArrayPointer_C(la - 1, l);

            for (int beta = 0; beta < (1 << (m - la)); ++beta) {
                c_lambda_1[2 * (2 * beta) + (psi % 2)] = (uint8_t) ((c_lambda[2 * beta] + c_lambda[2 * beta + 1]) % 2);
                c_lambda_1[2 * (2 * beta + 1) + (psi % 2)] = c_lambda[2 * beta + 1];
            }
        }
        if (psi % 2 == 1) {
            recursivelyUpdateC(la - 1, psi);
        }
    }

    void continuePaths_FrozenBit(int phi) {
        for (int l = 0; l < L; ++l) {
            if (pathIndexInactive(l)) {
                continue;
            }
            auto &C_m = getArrayPointer_C(m, l);
            C_m[phi % 2] = 0;
        }
    }

    void continuePaths_UnfrozenBit(int phi) {
        std::vector<std::vector<double>> probForks(L, std::vector<double>(2));
        int i = 0;
        for (int l = 0; l < L; ++l) {
            if (pathIndexInactive(l)) {
                probForks[l][0] = -1;
                probForks[l][1] = -1;
            } else {
                auto &P_m = getArrayPointer_P(m, l);
                probForks[l][0] = P_m[0];
                probForks[l][1] = P_m[1];
                ++i;
            }
        }
        int ro = std::min(2 * i, L);
        std::vector<std::vector<bool>> contForks(L, std::vector<bool>(2, false));
        populate_contForks(contForks, probForks, ro);
        for (int l = 0; l < L; ++l) {
            if (pathIndexInactive(l)) {
                continue;
            }
            if (!contForks[l][0] && !contForks[l][1]) {
                killPath(l);
            }
        }
        for (int l = 0; l < L; ++l) {
            if (!contForks[l][0] && !contForks[l][1]) {
                continue;
            }
            if (contForks[l][0] && contForks[l][1]) {
                getArrayPointer_C(m, l)[phi % 2] = 0;
                int l_ = clonePath(l);
                getArrayPointer_C(m, l_)[phi % 2] = 1;
                arrayPointer_I[l_] = arrayPointer_I[l];
                arrayPointer_I[l][phi] = 0;
                arrayPointer_I[l_][phi] = 1;
            } else {
                if (contForks[l][0]) {
                    getArrayPointer_C(m, l)[phi % 2] = 0;
                    arrayPointer_I[l][phi] = 0;
                } else {
                    getArrayPointer_C(m, l)[phi % 2] = 1;
                    arrayPointer_I[l][phi] = 1;
                }
            }
        }

    }

    static void
    populate_contForks(std::vector<std::vector<bool>> &contForks, std::vector<std::vector<double>> &probForks, int r) {
        std::vector<double> allProbForks;
        for (auto &v: probForks) {
            for (auto d: v) {
                allProbForks.push_back(d);
            }
        }
        int n = allProbForks.size();
        int index = n - r + 1;
        double first_r_value = util::select(allProbForks, 0, allProbForks.size(), index);
        int contForksNotZero = 0;
        for (int i = 0; i < contForks.size(); ++i) {
            for (int j = 0; j < contForks[i].size(); ++j) {
                contForks[i][j] = probForks[i][j] > first_r_value;
                contForksNotZero += contForks[i][j];
                assert(contForksNotZero < r);
            }
        }
        for (int i = 0; i < contForks.size(); ++i) {
            for (int j = 0; j < contForks[i].size(); ++j) {
                if (contForksNotZero >= r) {
                    return;
                }
                if (probForks[i][j] == first_r_value) {
                    contForks[i][j] = true;
                    ++contForksNotZero;
                }
            }
        }
    }

    int findMostProbablePath() {
        int l_ = 0;
        double p_ = 0;
        for (int l = 0; l < L; ++l) {
            if (pathIndexInactive(l)) {
                continue;
            }
            auto &C_m = getArrayPointer_C(m, l);
            auto &P_m = getArrayPointer_P(m, l);
            if (p_ < P_m[C_m[1]]) {
                l_ = l;
                p_ = P_m[C_m[1]];
            }
        }
        return l_;
    }

    bool isFrozen(int phi) {
        return std::find(frozen_bits_num.begin(), frozen_bits_num.end(), phi) != frozen_bits_num.end();
    }

    int L;
    int m;
    std::vector<int> inactivePathIndices;
    std::vector<std::vector<std::vector<double>>> arrayPointer_P;
    std::vector<std::vector<std::vector<uint8_t>>> arrayPointer_C;
    std::vector<std::vector<uint8_t>> arrayPointer_I;
    std::vector<std::vector<int>> arrayReferenceCount;
    std::vector<bool> activePath;
    std::vector<std::vector<int>> inactiveArrayIndices;
    std::vector<std::vector<int>> pathIndexToArrayIndex;
    std::vector<int> frozen_bits_num;
};


#endif //POLARCODING_POLARDECODER_HPP
