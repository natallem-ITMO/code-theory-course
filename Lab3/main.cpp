#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <algorithm>
#include <random>
#include <set>
#include <iomanip>

using namespace std;

constexpr bool verbose_debug = false;

int max_degree_eq_zero; //  == n

int add_alpha_degrees(int a, int b) {
    if (a == max_degree_eq_zero || b == max_degree_eq_zero) {
        return max_degree_eq_zero;
    }
    return (a + b) % max_degree_eq_zero;
}

int multiply_alpha_degrees(int a, int b) {
    if (a == max_degree_eq_zero || b == max_degree_eq_zero) {
        return max_degree_eq_zero;
    }
    return (a * b) % max_degree_eq_zero;
}

namespace show {

    void show_polinom(int p) {
        int t = 1;
        vector<int> res;
        while (t <= p) {
            res.push_back((p & t) != 0);
            t <<= 1;
        }
        bool first = true;
        for (int i = res.size() - 1; i >= 0; --i) {
            if (res[i]) {
                if (!first) {
                    cout << " + ";
                } else {
                    first = false;
                }
                cout << "x^" << i;
            }
        }
        cout << "\n";
    }

    void show_elements(vector<int> &elements) {
        for (int i = 0; i < elements.size(); ++i) {
            cout << "p(" << i << "): ";
            show_polinom(elements[i]);
        }
    }

    void show_polynomial(vector<int> &gen_pol, vector<int> &backward_elements, const std::string &message) {
        cout << message << '\n';
        bool first = true;
        for (int i = 0; i < gen_pol.size(); ++i) {
            int alpha_num = backward_elements[gen_pol[i]];
            int degree = gen_pol.size() - 1 - i;
            if (!first) {
                cout << " + ";
            } else {
                first = false;
            }
            cout << "a^" << alpha_num << "*y^" << degree;
        }
        cout << '\n';
    }

    void show_straight_polynomial(vector<int> &gen_pol, const std::string &message) {
        cout << message << '\n';
        bool first = true;
        for (int i = 0; i < gen_pol.size(); ++i) {
            if (!first) {
                cout << " + ";
            } else {
                first = false;
            }
            cout << "a^" << gen_pol[i] << "*y^" << i;
        }
        cout << '\n';
    }

}

bool is_root(int alpha, vector<int> &check_polynomial, vector<int> &elements,
             vector<int> &backward_elements, bool contains_degree = true) {
    int cur_value = max_degree_eq_zero;
    if (contains_degree) {
        for (int i = 0; i < check_polynomial.size(); ++i) {
            int t = i * alpha;
            int k = check_polynomial[i];
            int cur_alpha = add_alpha_degrees(i * alpha % max_degree_eq_zero, check_polynomial[i]); // todo no need
            cur_value = backward_elements[elements[cur_alpha] xor elements[cur_value]];
        }
    } else {
        for (int i = 0; i < check_polynomial.size(); ++i) {
//            multiply_alpha_degrees((check_polynomial.size() - 1 - i), alpha)
            int cur_alpha = add_alpha_degrees(multiply_alpha_degrees((check_polynomial.size() - 1 - i), alpha),
                                              backward_elements[check_polynomial[i]]);
            cur_value = backward_elements[elements[cur_alpha] xor elements[cur_value]];
        }
    }
    return cur_value == max_degree_eq_zero;
}


std::pair<vector<int>, vector<int>> find_elements(int p) {
    vector<int> elements(max_degree_eq_zero + 1);
    vector<int> backward_elements(max_degree_eq_zero + 1);
    int cur = 1;
    int i = 0;
    elements[i] = cur;
    backward_elements[cur] = i;
    for (i = 1; i < max_degree_eq_zero; ++i) {
        cur <<= 1;
        if (cur >= max_degree_eq_zero + 1) {
            cur ^= p;
        }
        elements[i] = cur;
        backward_elements[cur] = i;
    }
    elements.back() = 0;
    backward_elements[0] = max_degree_eq_zero;
    return {elements, backward_elements};
}

int error() {
    throw bad_alloc();
}

vector<int>
create_gen_polynomial(set<int> &all_alphas, vector<int> &elements, vector<int> &backward_elements) {
    if (all_alphas.empty()) {
        return {};
    }
//    all_alphas = {1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17};
    vector<int> res(2);
    res.reserve(all_alphas.size());
    res[0] = 1;
    auto itr = all_alphas.begin();
    res[1] = elements[*itr];
    vector<int> pushed_alphas(1, *itr);
    ++itr;
    for (; itr != all_alphas.end(); ++itr) {
//        show::show_polynomial(res, backward_elements, "cur res:");
        res.push_back(0);
//        show::show_polynomial(res, backward_elements, "shifted res:");
        vector<int> alpha_vec(res.size());
        alpha_vec[0] = 0;
        int cur_alpha = *itr;
        for (int i = 0; i < res.size() - 1; ++i) {
            int new_degree = add_alpha_degrees(backward_elements[res[i]], cur_alpha);
            alpha_vec[i + 1] = elements[new_degree];
        }
//        show::show_polynomial(alpha_vec, backward_elements, "alpha res:");
        for (int i = 0; i < res.size(); ++i) {
            res[i] = res[i] xor alpha_vec[i];
        }
        pushed_alphas.push_back(cur_alpha);
//        show::show_polynomial(res, backward_elements, "after mult on " + std::to_string(cur_alpha) + " have poly:");
        for (int prev_alpha: pushed_alphas) {
            if (!is_root(prev_alpha, res, elements, backward_elements, false)) {
                error();
            }
        }
    }
    return res;
}

template<class T, class U>
bool contains(const T &container, const U &value) {
    return container.find(value) != container.end();
}

auto find_gen_polynomial(int delta, vector<int> &elements, vector<int> &backward_elements, ostream &output) {
    set<int> all_alphas;
    vector<pair<int, vector<int>>> generating_polynomial;
    set<int> cur_class_alphas;
    for (int alpha = 1; alpha <= 1 + delta - 2; ++alpha) {
        int cur_alpha_degree = alpha;
        cur_class_alphas.clear();
        for (int i = 0; i < max_degree_eq_zero; ++i) {
            cur_alpha_degree %= max_degree_eq_zero;
            if (contains(all_alphas, cur_alpha_degree) && i == 0) {
                break;
            }
            if (contains(all_alphas, cur_alpha_degree)) {
                throw std::out_of_range("all alpha contains not first current alpha");
            }

            auto[itr, was_inserted] = cur_class_alphas.insert(cur_alpha_degree);
            if (!was_inserted) {
                break;
            }
            cur_alpha_degree <<= 1;
        }
        if (!cur_class_alphas.empty()) {
            std::string class_log;
            for (int i: cur_class_alphas) {
                class_log += to_string(i) + ", ";
                all_alphas.insert(i);
            }
            if constexpr(verbose_debug) { cout << "Classes: " + class_log.substr(0, class_log.length() - 2) + ", \n"; }
        }
    }
    int g_degree = all_alphas.size();
    int k = max_degree_eq_zero - g_degree;
    output << k << '\n';
//    auto gen = create_gen_polynomial(cur_class_alphas, elements, backward_elements);
    if constexpr(verbose_debug) {
        cout << "All alphas:\n";
        for (auto i: all_alphas) {
            cout << i << " ";
        }
        cout << '\n';
    }
    vector<int> gen_pol;
    if (all_alphas.empty()) {
        gen_pol = {1};
    } else {
        gen_pol = create_gen_polynomial(all_alphas, elements, backward_elements);
    }
    if constexpr(verbose_debug) {
        show::show_polynomial(gen_pol, backward_elements, "Gen pol:");
    }
    for (auto itr = gen_pol.rbegin(); itr != gen_pol.rend(); ++itr) {
        if (*itr != 0 && *itr != 1) {
            error();
        }
        output << *itr << " ";
    }
    output << '\n';
    return std::make_pair(gen_pol, k);
}

auto divide(vector<int> dividend, const vector<int> &divisor) {
    vector<int> quotient;
    vector<int> remainder;
    for (int i = 0; i + divisor.size() <= dividend.size(); ++i) {
        quotient.push_back(dividend[i]);
        if (dividend[i]) {
            for (int j = 0; j < divisor.size(); ++j) {
                dividend[i + j] = dividend[i + j] xor divisor[j];
            }
        }
    }
    for (int i = dividend.size() - divisor.size() + 1; i < dividend.size(); ++i) {
        remainder.push_back(dividend[i]);
    }
    return std::make_pair(quotient, remainder);
}

void encode(const vector<int> &to_encode, vector<int> &encoded, const vector<int> &gen) {
    for (int &i: encoded) {
        i = 0;
    }
    for (int i = 0; i < to_encode.size(); ++i) {
        encoded[i] = to_encode[i];
    }
    auto[quotient, remainder] = divide(encoded, gen);
    for (int i = 0; i < to_encode.size(); ++i) {
        encoded[i] = to_encode[i];
    }
    for (int i = to_encode.size(), j = 0; i < encoded.size(); ++i, ++j) {
        encoded[i] = remainder[j];
    }
}

int calc_syndrome(const vector<int> &vec, int alpha_degree, vector<int> &elements,
                  vector<int> &backward_elements) {
    int cur_value = 0;
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i]) {
            int cur_alpha_degree = multiply_alpha_degrees(alpha_degree, i);
            cur_value = cur_value xor elements[cur_alpha_degree];
        }
    }
    return backward_elements[cur_value];
}

void decode(vector<int> &to_decode, vector<int> &syndromes, vector<int> &elements, vector<int> &backward_elements) {
    for (int i = 0; i < syndromes.size(); ++i) {
        syndromes[i] = calc_syndrome(to_decode, i + 1, elements,
                                     backward_elements); // contains degree of alphas
    }
    vector<vector<int>> Ls(syndromes.size() + 1); // one L -> [0, 1, 3] -> a^0*x^0 + a^1*x^1 + a^3*x^2
    Ls[0] = {0};
    vector<int> deltas(syndromes.size() + 1);
    deltas[0] = {0};
    int m = 0;
    for (int r = 1; r <= syndromes.size(); ++r) {
        int delta_r = syndromes[r - 1];
        const auto &lambda_r_minus_1 = Ls[r - 1];
        for (int j = 1; j < lambda_r_minus_1.size(); ++j) {
            int cur_alpha_degree = add_alpha_degrees(lambda_r_minus_1[j], syndromes[r - j - 1]);
            delta_r = backward_elements[(elements[delta_r] xor elements[cur_alpha_degree])];
        }
        deltas[r] = delta_r;
        if (delta_r != max_degree_eq_zero) { // not zero
            int delta_r_divided_delta_m = delta_r - deltas[m];
            while (delta_r_divided_delta_m < 0) {
                delta_r_divided_delta_m += max_degree_eq_zero;
            }
            delta_r_divided_delta_m %= max_degree_eq_zero;

            vector<int> L_m_minus_1(r - m, max_degree_eq_zero); // store alpha degrees
            for (auto t: (m == 0) ? Ls[0] : Ls[m - 1]) {
                int cur_degree = add_alpha_degrees(t, delta_r_divided_delta_m);
                L_m_minus_1.push_back(cur_degree);
            }
            Ls[r] = Ls[r - 1];
            for (int i = 0; i < L_m_minus_1.size(); i++) {
                if (Ls[r].size() <= i) {
                    Ls[r].push_back(L_m_minus_1[i]);
                } else {
                    Ls[r][i] = backward_elements[(elements[Ls[r][i]] xor elements[L_m_minus_1[i]])];
                }
            }
            if (Ls[r].size() > Ls[m].size()) {
                m = r;
            }

        } else {
            Ls[r] = Ls[r - 1];
        }
        if constexpr(verbose_debug) {
            show::show_straight_polynomial(Ls[r], "delta_r = " + to_string(delta_r) + " Ls[" + to_string(r) + "]:");
        }
    }
    vector<int> &check_polynomial = Ls.back();
    vector<int> errors;
    for (int i = 0; i < max_degree_eq_zero; ++i) {
        if (is_root(i, check_polynomial, elements, backward_elements, max_degree_eq_zero)) {
            int opposit_element = i == 0 ? 0 : (max_degree_eq_zero - i);
            errors.push_back(opposit_element);
            to_decode[opposit_element] = to_decode[opposit_element] xor 1;
        }
    }
}

double run_simulation(vector<int> &to_encode, vector<int> &encoded, vector<int> &gen_polynomial,
                      vector<int> &error_encoded, vector<int> &syndromes, vector<int> &elements,
                      vector<int> &backward_elements, double error_prob, int max_error, int iteration_number,
                      mt19937 &generator, std::uniform_int_distribution<std::mt19937::result_type> &int_distribution,
                      uniform_real_distribution<> &double_distribution) {
    long cur_error_number = 0;
    long cur_iteration = 0;
    while (cur_iteration < iteration_number && cur_error_number < max_error) {
        for (int &i: to_encode) {
            i = int_distribution(generator);
        }
        encode(to_encode, encoded, gen_polynomial);
        for (int i = 0; i < encoded.size(); ++i) {
            double random = double_distribution(generator);
            if (random <= error_prob) {
                error_encoded[i] = encoded[encoded.size() - 1 - i] xor 1;
            } else {
                error_encoded[i] = encoded[encoded.size() - 1 - i];
            }
        }
        decode(error_encoded, syndromes, elements, backward_elements);
        for (int i = 0; i < encoded.size(); ++i) {
            if (encoded[i] != error_encoded[encoded.size() - 1 - i]) {
                ++cur_error_number;
                break;
            }
        }
        ++cur_iteration;
    }
    return double(cur_error_number) / double(cur_iteration);
}

void solve(const std::string &input, const std::string &output) {
    ifstream fin(input);
    ofstream fout(output);
    fout << std::setprecision(10);
    int n, p, delta;
    fin >> n >> p >> delta;
    max_degree_eq_zero = n;
    if constexpr(verbose_debug) {
        show::show_polinom(p);
    }
    auto[elements, backward_elements] = find_elements(p);
    if constexpr(verbose_debug) {
        show::show_elements(elements);
    }

    auto[gen_polynomial, k] = find_gen_polynomial(delta, elements, backward_elements, fout);
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::uniform_real_distribution<> double_distr(0, 1);
    std::uniform_int_distribution<std::mt19937::result_type> int_distr(0, 1);

    vector<int> encoded(n);
    vector<int> decoded(n);
    vector<int> to_encode(k);
    vector<int> syndromes(delta - 1);
    string command;
    while (fin >> command) {
        if (command == "Encode") {
            for (int i = k - 1; i >= 0; --i) {
                int t;
                fin >> t;
                to_encode[i] = t;
            }
            encode(to_encode, encoded, gen_polynomial);
            for (auto itr = encoded.rbegin(); itr != encoded.rend(); ++itr) {
                fout << *itr << " ";
            }
            fout << "\n";
        } else if (command == "Decode") {
            for (int &decoded_el: decoded) {
                fin >> decoded_el;
            }
            decode(decoded, syndromes, elements, backward_elements);
            for (auto t: decoded) {
                fout << (int) t << " ";
            }
            fout << "\n";
        } else if (command == "Simulate") {
            double error_prob;
            int iteration_number;
            int max_error;
            fin >> error_prob >> iteration_number >> max_error;
            fout << std::fixed
                 << run_simulation(to_encode, encoded, gen_polynomial, decoded,
                                   syndromes, elements, backward_elements,
                                   error_prob, max_error, iteration_number,
                                   generator, int_distr, double_distr) << '\n';
        } else if (command == "end") {
            break;
        }
    }
    fin.close();
    fout.close();
}

int main() {
    solve("input.txt", "output.txt");
}