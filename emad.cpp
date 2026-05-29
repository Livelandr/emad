//
// Copyright (C) 2026 Livelandr <livelandr@proton.me>
//
// This file is part of EMAD.
//
// EMAD is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// EMAD is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//

#include <iostream>
#include <fstream>
#include <vector>

double SMA(std::vector<double>& data, size_t period, size_t offset = 0) {
	double val = 0;

	for (size_t i = 0; i < period; i++) {
		val += data[data.size() - 1 - i - offset];
	}

	return val / period;
}

double WMA(std::vector<double>& data, size_t period, size_t offset = 0) {
	double val = 0;

	for (size_t i = 0; i < period; i++) {
		val += data[data.size() - 1 - i - offset] * (period - i);
	}

	return val /((period*(period+1.0))/2.0);
}

double HMA(std::vector<double>& data, size_t period, size_t offset = 0) {
	std::vector<double> HMA_Raw; 
	HMA_Raw.reserve(hma_period);

	if (data.size() < period + offset + hma_period) return 0.0;

	for (size_t i = 0; i < std::floor(std::sqrt(period)); i++) {
		double hma = 2 * WMA(data, period/2, offset + i) - WMA(data, period, offset + i);
		HMA_Raw.push_back(hma);
	}

	double val = 0.0;
	for (size_t i = 0; i < hma_period; i++) {
		val += HMA_Raw[i] * (hma_period - i);
	}
	return val / ((hma_period * (hma_period + 1.0)) / 2.0);
}

double EMAD(std::vector<double>& data, size_t period1, size_t period2, size_t delta, size_t offset = 0) {
	if (delta == 0) delta = 1;

	double diff = HMA(data, period1, offset) - HMA(data, period2, offset);
	double prevDiff = HMA(data, period1, delta + offset) - HMA(data, period2, delta + offset);

	double derivative = (diff - prevDiff) / delta;

	if (std::abs(derivative) < 1e-9) {
		return (diff >= 0) ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
	}

	return -diff / derivative;
}