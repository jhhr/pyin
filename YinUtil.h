/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    pYIN - A fundamental frequency estimator for monophonic audio
    Centre for Digital Music, Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef _YINUTIL_H_
#define _YINUTIL_H_

#include <cmath>

#include <iostream>
#include <vector>
#include <exception>

#include "vamp-sdk/FFT.h"

using std::vector;

class YinUtil
{
public:
    YinUtil(int yinBufferSize);
    ~YinUtil();
    
    double sumSquare(const double *in, int startInd, int endInd);
    void difference(const double *in, double *yinBuffer);
    void fastDifference(const double *in, double *yinBuffer);
    void slowDifference(const double *in, double *yinBuffer);
    void cumulativeDifference(double *yinBuffer);
    int absoluteThreshold(const double *yinBuffer, double thresh);
    vector<double> yinProb(const double *yinBuffer, int prior,
                           int minTau = 0, int maxTau = 0);

    /**
     * Candidate probabilities by the pYIN paper's rule, softened.
     * yinProb() gives every dip below each threshold that threshold's
     * probability, which lets the deeper dip at twice a voice's period
     * (or at two voices' common period) outweigh the voice's own. Here,
     * as YIN chooses, order counts: for each threshold, the dips below
     * it share its probability by their order in period, the k-th
     * getting a weight exp(-lambda * k). A dip that is not first still
     * gets some, so that a voice whose strong second harmonic makes a
     * first dip at half its period keeps some at its own. When no dip
     * is below a threshold, the deepest dip gets 0.01 of its
     * probability, as in yinProb().
     *
     * Dips closer together than mergeRadius (relative to the period of
     * the first of them) count as one, at the deepest of them: ripples
     * at the bottom of one valley are one candidate, and their first,
     * which the order would favour, is sharp.
     *
     * prior holds the probability of each of 100 thresholds,
     * 0.01 to 1.00.
     */
    vector<double> yinProbOrdered(const double *yinBuffer,
                                  const vector<double> &prior,
                                  double lambda, double mergeRadius);

    /**
     * The probabilities of the 100 thresholds 0.01 to 1.00 under a
     * Beta(a, b) distribution, evaluated as pYIN's own tables are (the
     * density at each threshold, normalised).
     */
    static vector<double> betaThresholdPrior(double a, double b);

    double parabolicInterpolation(const double *yinBuffer, int tau);

private:
    const int m_yinBufferSize;
    Vamp::FFTReal m_fft;
};

#endif
