#include <iostream>
#include <vector>
#include <cassert>
#include "PolarEncoder.hpp"
#include "PolarDecoder.hpp"

void testSelector() {
    std::vector<int> vec = {7, 10, 4, 3, 20, 15};
    assert(util::select<int>(vec, 0, vec.size(), 6) == 20);
    assert(util::select<int>(vec, 0, vec.size(), 3) == 7);
    std::vector<int> vec2 = {7, 10, 4, 3, 20, 15};
    std::vector<int> vec3 = {34};
    assert(util::select(vec2, 0, vec2.size(), 3) == 7);
    assert(util::select(vec2, 0, vec2.size(), 4) == 10);
    assert(util::select(vec2, 0, vec2.size(), 1) == 3);
    assert(util::select(vec2, 0, vec2.size(), 2) == 4);
    assert(util::select(vec3, 0, vec3.size(), 1) == 34);
}

void testPolarCodes(const std::vector<uint8_t> &information_seq, int n, int L) {
    PolarEncoder my_polar_encoder;
    PolarDecoder my_polar_decoder;
    int frozen_bits_num = (1 << n) - information_seq.size();
    std::vector<uint8_t> frozen_bits(frozen_bits_num, 0);
    std::vector<uint8_t> initial_word = information_seq;
    initial_word.insert(initial_word.end(), frozen_bits.begin(), frozen_bits.end());
    std::vector<uint8_t> coded_word = initial_word;
    my_polar_encoder.encode(coded_word);
    std::vector<int> frozen_nums;
    for (int i = information_seq.size(); i < initial_word.size(); ++i){
        frozen_nums.push_back(i);
    }
    std::vector<double> W_y_1;
    for (uint8_t t : coded_word){
        W_y_1.push_back(t);
    }
    std::vector<uint8_t> decoded_word = my_polar_decoder.SCL_decoder(L, n, W_y_1, frozen_nums);
    assert(information_seq == decoded_word);
}

void testDecoder(){
    testPolarCodes({1, 1, 0, 1}, 3, 5);
    testPolarCodes({1, 1, 0, 1}, 3, 10);
    testPolarCodes({1, 1, 0, 1}, 3, 3);
    testPolarCodes({1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0}, 4, 7);
    testPolarCodes({1, 0}, 4, 7);
    testPolarCodes({1, 1, 0, 1, 0, 1}, 3, 2);
    testPolarCodes({1, 1, 0, 0, 0, 0}, 3, 2);
    testPolarCodes({1, 1}, 1, 2);
    testPolarCodes({1, 1, 0, 1}, 2, 2);
}

int main(int argc, char *argv[]) {
    testDecoder();
    testSelector();
}