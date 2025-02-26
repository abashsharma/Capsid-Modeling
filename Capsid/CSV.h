#pragma once

#include <filesystem>
#include <fstream>

// Write data to CSV files
class CSVWriter
{
public:
    inline CSVWriter(const std::filesystem::path& ofile)
        : m_ofile(ofile)
    {
        if (!m_ofile.is_open())
            throw std::exception();
    }

    inline void write(auto t) 
    {
        m_ofile << t << ',';
    }

    // Write out anything with an iterator.
    template<typename T>
    inline void writeline(T begin, T end)
    {
        std::for_each(begin, end, [&] (auto v) {
            write(v);
        });
        m_ofile << '\n';
    }

private:
    std::ofstream m_ofile;
};
