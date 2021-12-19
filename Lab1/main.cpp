#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <algorithm>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

using line_t = vector<bool>;
using matrix_t = vector<line_t>;


void set_false_vector(vector<bool> &v) {
    for (int i = 0; i < v.size(); ++i) {
        v[i] = false;
    }
}

pair<int, int> activation_start_end(line_t &line) {
    int start = -1;
    int end = -1;
    for (int i = 0; i < line.size(); ++i) {
        if (line[i]) {
            if (start == -1) {
                start = i;
            }
            end = i;
        }
    }
    return {start, end};
}

void subtract(matrix_t &matrix, int l1, int l2) {
    for (int i = 0; i < matrix[l1].size(); ++i) {
        matrix[l1][i] = matrix[l1][i] xor matrix[l2][i];
    }
}

void simplify_matrix(matrix_t &m, int n, int k) {
    for (int t = 0; t < k; ++t) {
        int min_start_pos = INT_MAX;
        int min_start_pos_index = -1;
        for (int i = t; i < k; ++i) {
            auto[start, end] = activation_start_end(m[i]);
            if (min_start_pos > start) {
                min_start_pos = start;
                min_start_pos_index = i;
            }
        }
        if (min_start_pos_index != t) {
            subtract(m, t, min_start_pos_index);
        }
        for (int i = t + 1; i < k; ++i) {
            if (m[i][t]) {
                subtract(m, i, t);
            }
        }
    }
    vector<bool> processed_lines(k, false);
    for (int t = 0; t < k; ++t) {
        int max_end_pos = -1;
        int max_end_pos_index = -1;
        for (int i = k - 1; i >= 0; --i) {
            if (processed_lines[i]) {
                continue;
            }
            auto[start, end] = activation_start_end(m[i]);
            if (max_end_pos < end) {
                max_end_pos = end;
                max_end_pos_index = i;
            }
        }
        processed_lines[max_end_pos_index] = true;
        for (int i = max_end_pos_index - 1; i >= 0; --i) {
            if (processed_lines[i]) {
                continue;
            }
            if (m[i][max_end_pos]) {
                subtract(m, i, max_end_pos_index);
            }
        }
    }
}

struct Node {
    std::vector<bool> activated_values;
    std::vector<int> input_edges;
    std::vector<double> input_edges_value;
    int prev_index = -1;
    double max_error = INT_MIN;
};

std::vector<bool> convert(int x, int num) {
    std::vector<bool> ret;
    for (int i = 0; i < num; ++i) {
        if (x & 1)
            ret.push_back(true);
        else
            ret.push_back(false);
        x >>= 1;
    }
    std::reverse(ret.begin(), ret.end());
    return ret;
}

std::vector<std::vector<Node>>
build_lattice(std::vector<std::vector<bool>> &matrix, int n, int m, std::ostream &ostream) {
    std::vector<std::vector<Node>> level_nodes(m + 1);
    vector<int> starts(n);
    vector<int> ends(n);
    for (int i = 0; i < n; ++i) {
        auto[start, end] = activation_start_end(matrix[i]);
        starts[i] = start;
        ends[i] = end;
    }
    vector<vector<int>> level_activated_vertexes(m + 1);
    for (int level = 0; level <= m; ++level) {
        std::vector<int> &activated_vertexes = level_activated_vertexes[level];
        for (int i = 0; i < n; ++i) {
            int start = starts[i];
            int end = ends[i];
            if (start <= level - 1 && end >= level) {
                activated_vertexes.push_back(i);
            }
        }
        int vertexes_size = 1 << activated_vertexes.size();
        ostream << vertexes_size << " ";
        if (level == 0) {
            level_nodes[level].push_back(Node());
            continue;
        }
        for (int i = 0; i < vertexes_size; ++i) {
            level_nodes[level].push_back(Node());
            Node &cur_node = level_nodes[level].back();
            cur_node.activated_values = convert(i, activated_vertexes.size());
        }
        std::vector<Node> &nodes = level_nodes[level];
        std::vector<int> &prev_activated_vertexes = level_activated_vertexes[level - 1];
        for (int i = 0; i < nodes.size(); ++i) {
            Node &node = nodes[i];
            for (int j = 0; j < level_nodes[level - 1].size(); ++j) {
                Node &prev_node = level_nodes[level - 1][j];
                bool need_to_connect = true;
                int prev_ind = 0;
                int ind = 0;
                while (prev_ind < prev_activated_vertexes.size() && ind < activated_vertexes.size()) {
                    if (activated_vertexes[ind] == prev_activated_vertexes[prev_ind]) {
                        if (node.activated_values[ind] != prev_node.activated_values[prev_ind]) {
                            need_to_connect = false;
                            break;
                        }
                        ++ind;
                        ++prev_ind;
                    } else if (activated_vertexes[ind] < prev_activated_vertexes[prev_ind]) {
                        ++ind;
                    } else {
                        ++prev_ind;
                    }
                }
                ind = prev_ind = 0;
                if (need_to_connect) {
                    bool current_value = false;
                    for (int line_num = 0; line_num < n; ++line_num) {
                        while (ind < activated_vertexes.size() && activated_vertexes[ind] < line_num) {
                            ++ind;
                        }
                        while (prev_ind < prev_activated_vertexes.size() &&
                               prev_activated_vertexes[prev_ind] < line_num) {
                            ++prev_ind;
                        }
                        if ((ind < activated_vertexes.size() && activated_vertexes[ind] == line_num &&
                            node.activated_values[ind])
                            || (prev_ind < prev_activated_vertexes.size() &&
                               prev_activated_vertexes[prev_ind] == line_num && prev_node.activated_values[prev_ind])) {
                            current_value = matrix[line_num][level - 1] ? !current_value : current_value;
                        }
                    }
                    node.input_edges.push_back(j);
                    node.input_edges_value.push_back(current_value ? -1. : 1.);
                }
            }
        }
    }
    ostream << '\n';
    return level_nodes;
}

void clean_nodes(std::vector<std::vector<Node>> &level_nodes) {
    for (auto &v: level_nodes) {
        for (auto &k: v) {
            k.max_error = INT_MIN;
            k.prev_index = -1;
        }
    }
}

void decode(std::vector<std::vector<Node>> &level_nodes, std::vector<double> &y, int k, int n, std::vector<bool> &res) {
    clean_nodes(level_nodes);
    level_nodes[0][0].max_error = {0};
    for (int level = 1; level <= n; ++level) {
        for (int j = 0; j < level_nodes[level].size(); ++j) {
            Node &node = level_nodes[level][j];
            for (int i = 0; i < node.input_edges.size(); ++i) {
                double edge_value = node.input_edges_value[i];
                double error_rate = edge_value * y[level - 1];
                double prev_error = level_nodes[level - 1][node.input_edges[i]].max_error;
                double cur_error = prev_error + error_rate;
                if (node.prev_index == -1 || node.max_error < cur_error) {
                    node.prev_index = i;
                    node.max_error = cur_error;
                }
            }
        }
    }
    int node_index = 0;
    for (int level = n; level > 0; --level) {
        Node &node = level_nodes[level][node_index];
        res[level - 1] = node.input_edges_value[node.prev_index] == -1;
        node_index = node.input_edges[node.prev_index];
    }
}


void encode(const vector<bool> &v, const matrix_t &m, int n, int k, vector<bool> &res) {
    for (int j = 0; j < k; ++j) {
        if (v[j]) {
            for (int i = 0; i < n; ++i) {
                res[i] = res[i] xor m[j][i];
            }
        }
    }
}

double
run_simulation(const double noise_level, const long iteration_number, const long max_error,
               std::vector<std::vector<Node>> &level_nodes,
               int n, int k, const matrix_t &m, std::mt19937 &generator,
               std::uniform_int_distribution<std::mt19937::result_type> &uniform_distribution, double R, vector<bool> &origin,
               vector<bool> &origin_encoded,
               vector<double> &noised_encoded, vector<bool> &decoded) {
//    double sigma = std::sqrt(1. / (2. * R * std::pow(10., noise_level / 10.)));
    double sigma = std::sqrt(
    0.5 * pow(10, (-noise_level) / 10) / ((double) R));
    long cur_error_number = 0;
    long cur_iteration = 0;
    std::normal_distribution<double> gauss_dist(0.0f, sigma);
    while (cur_iteration < iteration_number && cur_error_number < max_error) {
        for (int i = 0; i < k; ++i) {
            origin[i] = uniform_distribution(generator) == 1;
        }
        set_false_vector(origin_encoded);
        encode(origin, m, n, k, origin_encoded);
        for (int i = 0; i < n; ++i) {
            noised_encoded[i] = (origin_encoded[i] ? -1. : 1.) + gauss_dist(generator);
        }
        decode(level_nodes, noised_encoded, k, n, decoded);
        if (decoded != origin_encoded) {
            ++cur_error_number;
        }
        ++cur_iteration;
    }
    return double(cur_error_number) / double(cur_iteration);
}

int main() {
//    unsigned int start_time =  clock();
//    ios_base::sync_with_stdio(NULL);
//    cin.tie(NULL);
//    cout.tie(NULL);
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    fout << std::setprecision(10);
    int n, k;
    fin >> n >> k;
    double R = double(k) / n;
    matrix_t init_matrix(k, line_t(n, false));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            int t;
            fin >> t;
            init_matrix[i][j] = t == 1;
        }
    }
    matrix_t matrix = init_matrix;
    simplify_matrix(matrix, n, k);
    std::vector<std::vector<Node>> level_nodes = build_lattice(matrix, k, n, fout);
    string command;
//    std::default_random_engine generator;
//    std::bernoulli_distribution bernoulliDistribution;
    std::random_device rd{};
    std::mt19937 generator{rd()};
    std::uniform_int_distribution<std::mt19937::result_type> uniform_distr(0, 1);

    vector<bool> encoded_vector(n);
    vector<bool> decoded_vector(n);
    vector<double> noised_encoded(n);
    vector<bool> origin(k);
    while (fin >> command) {
        if (command == "Encode") {
            for (int i = 0; i < k; ++i) {
                int t;
                fin >> t;
                origin[i] = t == 1;
            }
            set_false_vector(encoded_vector);
            encode(origin, init_matrix, n, k, encoded_vector);
            for (auto t: encoded_vector) {
                fout << (int) t << " ";
            }
            fout << "\n";
        } else if (command == "Decode") {
            for (int i = 0; i < n; ++i) {
                fin >> noised_encoded[i];
            }
            decode(level_nodes, noised_encoded, k, n, decoded_vector);
            for (auto t: decoded_vector) {
                fout << (int) t << " ";
            }
            fout << "\n";
        } else if (command == "Simulate") {
            int noise_level;
            int iteration_number;
            int max_error;
            fin >> noise_level >> iteration_number >> max_error;
            fout << std::fixed << (double) run_simulation(noise_level, iteration_number, max_error, level_nodes, n, k, init_matrix,
                                   generator, uniform_distr, R, origin, encoded_vector, noised_encoded,
                                   decoded_vector) << '\n';
        }
    }
//    unsigned int end_time = clock(); // конечное время
//    unsigned int search_time = end_time - start_time;
//    std::cout << "time - " << search_time;
}
