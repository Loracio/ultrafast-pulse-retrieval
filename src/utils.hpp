#ifndef UTILS_INCLUDED
#define UTILS_INCLUDED

/**
 * @file utils.hpp
 * @author Víctor Loras Herrero
 * @brief Various auxiliary generic functions used in other modules.
 * 
 */

#include <cmath>
#include <numeric>
#include <vector>
#include <algorithm>


/**
 * \brief Computes the mean value of the element-wise division of two vectors.
 *
 * This function computes the mean value of the element-wise division of two vectors,
 * `x` and `y`. It calculates the sum of the element-wise product of `x` and `y`,
 * and divides it by the sum of the elements in `y`.
 *
 * \tparam T The type of elements in the vectors.
 * \param x The first vector.
 * \param y The second vector.
 * \return The mean value of the element-wise division of `x` and `y`.
 */
template <typename T>
T mean(const std::vector<T> &x, const std::vector<T> &y)
{
    T sumProduct = std::inner_product(x.begin(), x.end(), y.begin(), 0);
    T sumY = std::accumulate(y.begin(), y.end(), 0);

    return sumProduct / sumY;
}


/**
 * \brief Computes the standard deviation of a dataset.
 *
 * This function calculates the standard deviation of a dataset represented by two vectors,
 * `x` and `y`. It uses the mean value of `x` and `y` to calculate the standard deviation.
 *
 * \tparam T The type of elements in the vectors.
 * \param x The vector representing the x-values of the dataset.
 * \param y The vector representing the y-values of the dataset.
 * \return The standard deviation of the dataset.
 */
template <typename T>
T stdDev(const std::vector<T>& x, const std::vector<T>& y)
{
    T meanValue = mean(x, y);

    T sumNum = 0;
    T sumY = 0;

    for (size_t i = 0; i < x.size(); i++)
    {
        T diff = x[i] - meanValue;
        sumNum += diff * diff * y[i];
        sumY += y[i];
    }

    return std::sqrt(sumNum / sumY);
}


/**
 * \brief Finds the maximum element in a vector.
 *
 * This function finds and returns the maximum element in the given vector.
 *
 * \tparam T The type of elements in the vector.
 * \param vec The vector to search for the maximum element.
 * \return The maximum element in the vector.
 */
template <typename T>
T findMax(const std::vector<T> &vec)
{
    return *std::max_element(vec.begin(), vec.end());
}

/**
 * \brief Finds the index where the value in the vector exceeds a specified threshold starting from the beggining of the vector.
 *
 * This function searches the given vector for the leftmost index where the value exceeds
 * the specified threshold.
 *
 * \tparam T The type of elements in the vector.
 * \param x The vector to search for the left index.
 * \param value The threshold value to compare with.
 * \return The left index where the value exceeds the threshold, or -1 if not found.
 */
template <typename T>
int leftIndexValue(const std::vector<T> &x, double value)
{
    int N = x.size();

    for (int i = 0; i < N; i++)
    {
        if (x[i] - value > 0)
        {
            return i;
        }
    }

    return -1;
}

/**
 * \brief Finds the first index where the value in the vector exceeds a specified threshold starting from the end of the vector.
 *
 * This function searches the given vector for the rightmost index where the value exceeds
 * the specified threshold.
 *
 * \tparam T The type of elements in the vector.
 * \param x The vector to search for the right index.
 * \param value The threshold value to compare with.
 * \return The right index where the value exceeds the threshold, or -1 if not found.
 */
template <typename T>
int rightIndexValue(const std::vector<T> &x, double value)
{
    int N = x.size();

    for (int i = 0; i < N; i++)
    {
        if (x[N - i - 1] - value > 0)
        {
            return N - i - 1;
        }
    }

    return -1;
}

/**
 * \brief Calculates the Full Width at Half Maximum (FWHM) of a waveform.
 *
 * This function calculates the Full Width at Half Maximum (FWHM) of a waveform
 * given a vector of values. The FWHM is the width of the waveform at half of its maximum value.
 *
 * \tparam T The type of elements in the waveform vector.
 * \tparam Numeric The type of the delta_t parameter.
 * \param x The vector of waveform values.
 * \param delta_t The time interval between waveform samples.
 * \return The FWHM of the waveform.
 */
template <typename T, typename Numeric>
T FWHM(const std::vector<T> &x, Numeric delta_t)
{
    T halfMax = findMax(x) / 2;

    T leftIndex = leftIndexValue(x, halfMax);
    T rightIndex = rightIndexValue(x, halfMax);

    return (rightIndex - leftIndex) * delta_t;
}

/**
 * Converts a vector of frequencies to angular frequency units.
 *
 * @param frequencies The vector of frequencies to convert.
 * @return The vector of frequencies in angular frequency units.
 */
std::vector<double> toAngularFrequency(const std::vector<double>& frequencies)
{
    std::vector<double> angularFrequencies(frequencies.size());

    for (size_t i = 0; i < frequencies.size(); i++)
    {
        angularFrequencies[i] = 2 * M_PI * frequencies[i];
    }

    return angularFrequencies;
}

#endif // UTILS_INCLUDED