/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2016 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#pragma once

#include <random>
#include <stdexcept>
#include <vector>

inline std::vector<double> generateGaussNoise(const std::vector<double>& data,
                                              const std::vector<double>& error,
                                              std::mt19937& rng) {
    if (data.size() != error.size())
        throw std::invalid_argument("data and error vectors must have equal size");

    std::vector<double> noisy_data(data.size());
    for (std::size_t i = 0; i < data.size(); ++i) {
        if (error[i] < 0.0)
            throw std::invalid_argument("Gaussian standard deviation must be nonnegative");
        if (error[i] == 0.0)
            noisy_data[i] = data[i];
        else {
            std::normal_distribution<> distribution(data[i], error[i]);
            noisy_data[i] = distribution(rng);
        }
    }
    return noisy_data;
}
