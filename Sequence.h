#pragma once


// Generate pairs of [C0, Ct] sequences for multithreaded sampling
// Currently it just generates sequences in a line, but it can be altered
//  with a different sampling method in the future.
class Sequence 
{ 
public:
    using iterator_category = std::input_iterator_tag;
    using value_type        = std::pair<double, double>;
    using reference         = const value_type&;
    using pointer           = value_type const*;
    using difference_type   = ptrdiff_t;

public:
    Sequence(double C0, double Ct, double CStep)
        : maxC0(C0), maxCt(Ct), cStep(CStep), mCur(CStep, CStep)
    {
    }

    reference operator*() const
    {
        return mCur;
    }

    pointer operator->() const
    {
        return &mCur;
    }

    explicit operator bool() const
    {
        return mCur.first < maxCt;
    }

    Sequence& operator++()
    {
        if (mCur.second >= maxCt)
        {
            mCur.first += cStep;
            mCur.second = cStep;
        }
        else
        {
            mCur.second += cStep;
        }
        return *this;
    }

    Sequence operator++(int)
    {
        const Sequence tmp(*this);
        ++*this;
        return tmp;
    }

private:
    double maxC0, maxCt, cStep;

    value_type mCur;

};

