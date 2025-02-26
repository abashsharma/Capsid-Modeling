#pragma once

#include <string_view>
#include "Capsid.h"

void optimize(capsid::Harmonics& h, std::string_view method);

double std_bending_energy(const capsid::Harmonics& h, double K);

double capsid_energy_function(const capsid::Harmonics& h);
