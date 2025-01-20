/*
 * This file is part of the Nano-Shell Simulation Project.
 * 
 * Copyright (C) 2025 Alessandro Veltri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cmath>

int main() {
    const int numPoints = 400;
    const double rangeMin = 0.0;
    const double rangeMax = 2.0;
    const double peakValue = 1.0;
    const double sigma = 0.05; // Controls the density around the peak

    // Calculate the step size
    const double stepSize = (rangeMax - rangeMin) / (numPoints - 1);

    for (int i = 0; i < numPoints; ++i) {
        double x = rangeMin + i * stepSize;
        double distance = std::abs(x - peakValue);
        double density = std::exp(-0.5 * std::pow(distance / sigma, 2));

        // Adjust step size based on density
        double adjustedStepSize = stepSize / density;

        // Calculate the new position
        x = rangeMin + i * adjustedStepSize;

        std::cout <<  x << " " << 1. << "\n";
    }

    return 0;
}
