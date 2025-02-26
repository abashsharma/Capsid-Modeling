#pragma once

#include <filesystem>
#include "Capsid.h"

namespace fs = std::filesystem;

namespace capsid
{

void SaveRadii(const Harmonics& h, const fs::path& ofile, bool append=false);

const char* SupportedFormats() noexcept;

} // namespace capsid
