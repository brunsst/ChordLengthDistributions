#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <iostream>
#include <string.h>
#include <cstdint>
#include <vector>

namespace aux
{
    //2D-Conversions
    void idx2xy(const int64_t &pos, const int shape[2], int &outx, int &outy);
    int64_t xy2idx(const int &x,const int &y, const int shape[2]);

    //3D-Conversions
    void idx2xyz(const int64_t &pos, const int shape[3], int &outx, int &outy, int &outz);
    int64_t xyz2idx(const int &x, const int &y, const int &z, const int shape[3]);

    /*String-Manipulation
    *********************************************************/
    std::string zfill_int2string(int inint, const unsigned int &zfill);

    /*Numpy-like
    *********************************************************/
    std::vector<double> linspace(double startval, double endval, uint64_t bins);
    void rescale(std::vector<float> &imgdata, float lb, float ub);
    void rescale(std::vector<float> &imgdata, float lb, float ub, float maxval, float minval);
    void rescale(std::vector<float> &imgdata, float lb, float ub, float minval, float maxval, bool cutoffoob);
    std::vector<size_t> argsort_ascending(const std::vector<uint64_t> &v);
    std::vector<size_t> argsort_descending(const std::vector<uint64_t> &v);

    /*
    *********************************************************/

    // Get the median of an unordered set of numbers of arbitrary
    // type without modifying the underlying dataset.
    template <typename It>
    typename std::iterator_traits<It>::value_type median(It begin, It end)
    {
        using T = typename std::iterator_traits<It>::value_type;
        std::vector<T> data(begin, end);
        std::nth_element(data.begin(), data.begin() + data.size() / 2, data.end());
        return data[data.size() / 2];
    }
}

#endif // AUXILIARY_H



