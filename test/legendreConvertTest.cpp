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

#include "../legendre_convert/gaussian_noise.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

TEST(LegendreConvert, ZeroDeviationDoesNotAdvanceGenerator) {
  std::mt19937 rng(1234);
  const std::mt19937 untouched = rng;
  const std::vector<double> data{1.0, -2.0};

  EXPECT_EQ(generateGaussNoise(data, {0.0, 0.0}, rng), data);
  EXPECT_EQ(rng, untouched);
}

TEST(LegendreConvert, PositiveDeviationUsesDistribution) {
  std::mt19937 rng(1234);
  const std::mt19937 untouched = rng;

  const auto noisy = generateGaussNoise({1.0}, {0.5}, rng);

  ASSERT_EQ(noisy.size(), 1U);
  EXPECT_TRUE(std::isfinite(noisy[0]));
  EXPECT_NE(rng, untouched);
}
