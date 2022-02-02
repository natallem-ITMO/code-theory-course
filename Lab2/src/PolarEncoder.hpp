//
// Created by natallem on 10/10/21.
//

#ifndef POLARCODING_POLARENCODER_HPP
#define POLARCODING_POLARENCODER_HPP

#include <vector>
#include <cstdint>

class PolarEncoder {

public:
    void encode(std::vector<uint8_t> &u) {
        for (int iteration = 1; iteration < u.size(); iteration <<= 1){
            for (int j = 0; j < iteration; ++j) {
                for (int i = 0; i < u.size(); i += 2 * iteration) {
                    u[i + j] = (uint8_t)((u[i + j] + u[i + j + iteration]) % 2);
                }
            }
        }
        return;
       /* if (end - start == 1) {
            return;
        }
        int half_elements_number = (end - start) / 2;
        std::vector<uint8_t> copy;

        for (int i = start; i < end; i += 2) {
            copy.push_back(u[i + 1]);
            u[i] ^= u[i + 1];
        }
        for (int i = 0; i < half_elements_number; ++i) {
            u[start + i] = u[start + i * 2];
        }
        for (int i = 0; i < half_elements_number; ++i) {
            u[start + half_elements_number + i] = copy[i];
        }
        encode(u, start, end - half_elements_number);
        encode(u, start + half_elements_number, end);*/
    }

};


#endif //POLARCODING_POLARENCODER_HPP
