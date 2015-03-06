/*-----------------------------------------------------------------------------
@file

Declaration of filter coefficients, correlator coefficients and other tables

------------------------------------------------------------------------------*/



#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H


#include <array>


// Coefficients for the root raised cosine filter.
// The parameters are:
// Sampling frequency 38400 sps
// Symbol rate: 4800
// Filter length: 8 symbols
// Oversampling: 8x

extern const std::array<double, 65> coeffsRrc;



#endif