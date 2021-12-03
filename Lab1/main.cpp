#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <climits>
#include <algorithm>
#include <cmath>
#include <random>

using namespace std;

using line_t = vector<bool>;
using matrix_t = vector<line_t>;

const bool verbose_debug = false;

pair<int, int> activation_start_end(line_t &line) {
    int start = -1;
    int end = -1;
    for (int i = 0; i < line.size(); ++i) {
        if (line[i]) {
            start = i;
            break;
        }
    }
    for (int i = line.size() - 1; i >= 0; --i) {
        if (line[i]) {
            end = i;
            break;
        }
    }
    return {start, end};
}

void show_matrix(matrix_t &m, std::ostream &ostream = std::cout) {
    if constexpr(verbose_debug) {
        for (int i = 0; i < m.size(); ++i) {
            for (int j = 0; j < m[i].size(); ++j) {
                ostream << m[i][j];
            }
            auto[start, end] = activation_start_end(m[i]);
            ostream << " (" << i << ") {" << start << ", " << end << "}\n";
        }
        ostream << '\n';
    }
}

void subtract(matrix_t &matrix, int l1, int l2) {
    for (int i = 0; i < matrix[l1].size(); ++i) {
        matrix[l1][i] = ((!matrix[l1][i]) != (!matrix[l2][i]));
    }
}


void simplify_matrix(matrix_t &m, int n, int k) {
    for (int t = 0; t < k; ++t) {
        int min_start_pos = INT_MAX;
        int min_start_pos_index = -1;
//        vector<int> start_pos(k);
        for (int i = t; i < k; ++i) {
            auto[start, end] = activation_start_end(m[i]);
//            start_pos[i] = start;
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
//        vector<int> end_pos(k);
        for (int i = k - 1; i >= 0; --i) {
            if (processed_lines[i]) {
                continue;
            }
            auto[start, end] = activation_start_end(m[i]);
//            end_pos[i] = end;
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
    int num;
    std::vector<int> activated_vertexes;
    std::vector<bool> activated_values;
    std::vector<Node *> input_edges;
    std::vector<Node *> output_edges;
    std::vector<double> output_edges_value;
    std::vector<double> input_edges_value;
    std::vector<double> error_input_edge;

    std::string toString() {
        std::string res = "activated nums: ";
        for (auto t: activated_vertexes) {
            res += std::to_string(t) + " ";
        }
        res += " values: ";
        for (auto t: activated_values) {
            res += std::to_string(t);
        }
        res += ", edges: ";
        int alignment_size = res.size();
        for (int i = 0; i < input_edges.size(); ++i) {
            if (i != 0) {
                res += "\n\t" + string(alignment_size, ' ');
            }
            res += "input edge num(" + std::to_string(input_edges[i]->num) + ")(edge=" +
                   to_string(input_edges_value[i]) + ")=" +
                   std::to_string(error_input_edge[i]) + " ";
        }
        return res;
    }

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

std::vector<std::vector<Node *>>
build_lattice(std::vector<std::vector<bool>> &matrix, int n, int m, std::ostream &ostream) {
    static const bool debug = false;
    std::vector<std::vector<Node *>> level_nodes;
    vector<int> starts(n);
    vector<int> ends(n);
    for (int i = 0; i < n; ++i) {
        auto[start, end] = activation_start_end(matrix[i]);
        starts[i] = start;
        ends[i] = end;
    }
    for (int level = 0; level <= m; ++level) {
        std::vector<int> activated_vertexes;
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
            Node *first_node = new Node();
            level_nodes.push_back({first_node});
            continue;
        }
        std::vector<Node *> nodes;
        for (int i = 0; i < vertexes_size; ++i) {
            Node *cur_node = new Node();
            cur_node->activated_vertexes = activated_vertexes;
            cur_node->activated_values = convert(i, activated_vertexes.size());
            cur_node->num = i;
            nodes.push_back(cur_node);
        }

        for (Node *node: nodes) {
            for (Node *prev_node: level_nodes[level - 1]) {
                bool need_to_connect = true;
                for (int i = 0; i < prev_node->activated_vertexes.size(); ++i) {
                    int prev_active_vertex_num = prev_node->activated_vertexes[i];
                    bool prev_active_vertex_value = prev_node->activated_values[i];
                    auto it = std::find(activated_vertexes.begin(), activated_vertexes.end(), prev_active_vertex_num);
                    if (it != activated_vertexes.end()) {
                        int index = std::distance(activated_vertexes.begin(), it);
                        if (node->activated_values[index] != prev_active_vertex_value) {
                            need_to_connect = false;
                            goto checked_node;
                        }
                    }
                }
                checked_node:
                if (need_to_connect) {
                    bool current_value = false;
                    for (int line_num = 0; line_num < n; ++line_num) {
                        auto it = std::find(activated_vertexes.begin(), activated_vertexes.end(), line_num);
                        auto it_prev = std::find(prev_node->activated_vertexes.begin(),
                                                 prev_node->activated_vertexes.end(), line_num);
                        bool is_line_active = false;
                        if (it != activated_vertexes.end()) {
                            int index = std::distance(activated_vertexes.begin(), it);
                            is_line_active = node->activated_values[index];
                        } else if (it_prev != prev_node->activated_vertexes.end()) {
                            size_t index = std::distance(prev_node->activated_vertexes.begin(), it_prev);
                            is_line_active = prev_node->activated_values[index];
                        }
                        if (is_line_active) {
                            current_value = matrix[line_num][level - 1] ? !current_value : current_value;
                        }
                    }
                    node->input_edges.push_back(prev_node);
                    node->input_edges_value.push_back(current_value ? -1. : 1.);
                    prev_node->output_edges.push_back(node);
                    prev_node->output_edges_value.push_back(current_value ? -1. : 1.);
                }
            }
        }
        level_nodes.push_back(std::move(nodes));
    }
    ostream << '\n';
    return level_nodes;
}

void clean_nodes(std::vector<std::vector<Node *>> &level_nodes) {
    for (auto &v: level_nodes) {
        for (auto &k: v) {
            k->error_input_edge.clear();
        }
    }
}

std::vector<bool> decode(std::vector<std::vector<Node *>> &level_nodes, std::vector<double> &y, int k, int n) {
    clean_nodes(level_nodes);
    std::vector<double> c;
    level_nodes[0][0]->error_input_edge = {0};
    for (int level = 1; level <= n; ++level) {
        for (Node *node: level_nodes[level]) {
            int node_num = node->num;
            for (int i = 0; i < node->input_edges.size(); ++i) {
                double edge_value = node->input_edges_value[i];
                double error_rate = edge_value * y[level - 1];
                double prev_error = *std::max_element(std::begin(node->input_edges[i]->error_input_edge),
                                                      std::end(node->input_edges[i]->error_input_edge));
                node->error_input_edge.push_back(prev_error + error_rate);
            }
        }
    }
    Node *node = level_nodes[n][0];
    for (int level = n; level > 0; --level) {
        auto max_it = std::max_element(std::begin(node->error_input_edge),
                                       std::end(node->error_input_edge));
        int index = std::distance(node->error_input_edge.begin(), max_it);
        c.push_back(node->input_edges_value[index]);
        node = node->input_edges[index];
    }
    vector<bool> result;
    result.reserve(c.size());
    for (auto itr = c.rbegin(); itr != c.rend(); ++itr) {
        result.push_back(*itr == -1);
    }
    return result;
}

template<typename T>
void show_vector(vector<T> &vec) {
    for (auto t: vec) {
        cout << t << " ";
    }
}

vector<bool> encode(const vector<bool> &v, const matrix_t &m, int n, int k) {
    vector<bool> res(n, false);
    for (int j = 0; j < k; ++j) {
        if (v[j]) {
            for (int i = 0; i < n; ++i) {
                res[i] = ((!res[i]) != (!m[j][i]));
            }
        }
    }
    return res;
}


double
run_simulation(const double noise_level, const long iteration_number, const long max_error,
               std::vector<std::vector<Node *>> &level_nodes,
               int n, int k, const matrix_t &m, std::default_random_engine &generator, double R) {
    double sigma = std::sqrt(1. / (2. * R * std::pow(10., noise_level / 10.)));
    long cur_error_number = 0;
    long cur_iteration = 0;
    std::normal_distribution<double> gauss_dist(0.0f, sigma);
    while (cur_iteration < iteration_number && cur_error_number < max_error) {
        vector<bool> origin(k);
        for (int i = 0; i < k; ++i) {
            origin[i] = rand() % 2;
        }
        vector<bool> origin_encoded = encode(origin, m, n, k);
        vector<double> noised_origin(n);
        for (int i = 0; i < n; ++i) {
            noised_origin[i] = (origin_encoded[i] ? -1. : 1.) + gauss_dist(generator);
        }
        vector<bool> decoded = decode(level_nodes, noised_origin, k, n);
        for (int i = 0; i < n; ++i) {
            if (origin_encoded[i] != decoded[i]) {
                ++cur_error_number;
                break;
            }
        }
        ++cur_iteration;
    }
    return double(cur_error_number) / double(cur_iteration);
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int n, k;
    fin >> n >> k;
    double R = double(k) / n;
    matrix_t init_matrix(k, line_t(n, false));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < n; ++j) {
            int t;
            fin >> t;
            if (t == 1) {
                init_matrix[i][j] = true;
            }
        }
    }
    matrix_t matrix = init_matrix;
    show_matrix(matrix, fout);
    simplify_matrix(matrix, n, k);
    show_matrix(matrix, fout);
    std::vector<std::vector<Node *>> level_nodes = build_lattice(matrix, k, n, fout);
    string command;
    std::default_random_engine generator;
    while (fin >> command) {
        if (command == "Encode") {
            vector<bool> to_encode(k);
            for (int i = 0; i < k; ++i) {
                int t;
                fin >> t;
                to_encode[i] = t == 1;
            }
            vector<bool> multiplication = encode(to_encode, init_matrix, n, k);
            for (auto t: multiplication) {
                fout << (int) t << " ";
            }
            fout << "\n";
        } else if (command == "Decode") {
            double d;
            vector<double> y(n, false);
            for (int i = 0; i < n; ++i) {
                fin >> y[i];
            }
            std::vector<bool> c = decode(level_nodes, y, k, n);
            for (auto t: c) {
                fout << (int) t << " ";
            }
            fout << "\n";
        } else if (command == "Simulate") {
            double noise_level;
            int iteration_number;
            int max_error;
            fin >> noise_level >> iteration_number >> max_error;
            double res = run_simulation(noise_level, iteration_number, max_error, level_nodes, n, k, init_matrix,
                                        generator, R);
            fout << res << '\n';
        } else {
            fout << "unknown command";
        }
    }
#ifdef DEBUG
    fout.close();
    ifstream t("output.txt");
    t.seekg(0, std::ios::end);
    size_t size = t.tellg();
    std::string buffer(size, ' ');
    t.seekg(0);
    t.read(&buffer[0], size);
    std::string expected = "1 2 4 8 4 8 4 2 1 \n"
                           "1 1 1 1 1 1 1 1 \n"
                           "0 0 1 1 1 1 0 0 \n"
                           "0.0264131\n"
                           "0.00830565\n";
    if (buffer != expected) {
        std::cout << "Output not equal to origin. Expected:\n" << expected << "Got:\n" << buffer;
        return 1;
    }
#endif
    return 0;
}
