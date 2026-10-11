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

#ifndef _HARMONICSALIENCE_H_
#define _HARMONICSALIENCE_H_

#include "vamp-sdk/FFT.h"

#include <utility>
#include <vector>

/**
 * F0 candidates from a harmonic salience function of the spectrum.
 *
 * YIN's difference function of two simultaneous voices (a harmony,
 * a doubled part, the echo of the last note under the next) dips
 * deepest at the period the two have in common, which is far below
 * either of them. A sum of the spectrum's peaks at the harmonics of
 * each trial f0 instead peaks at each voice's own f0: the common
 * period collects only the few partials that fall on its harmonics,
 * with the small weights of high harmonic numbers.
 *
 * Salience of a trial f0, on a grid of 10-cent bins:
 *   S(f0) = sum_h w_h (P(h f0) - negWeight P((h - 1/2) f0)),
 *   w_h = alpha^(h-1)
 * where P is the peaks' compressed magnitude (a^beta) spread over a
 * cos^2 tolerance window. The negative lobes between the harmonics
 * (as in SWIPE) count against an octave-up trial f0, whose (h - 1/2)
 * multiples are the true f0's odd harmonics.
 *
 * The host's blocks overlap, so the class keeps the stream of
 * samples and takes a window centred on each frame's time (a given
 * offset from the start of its block), which may reach past the
 * block. A frame whose window reaches past the samples given so far
 * waits: the caller pushes each block, then takes the frames that
 * are ready, and at the end takes the rest (zeros past the end of
 * the stream).
 */
class HarmonicSalience
{
public:
    /**
     * The analysis settings. The plugin uses these defaults, tuned in
     * Tony's pitch lab on singing with harmonies, echoes and the damage
     * a separation from a mix does.
     */
    struct Params {
        int windowSize = 2048;    // samples, Hann
        int fftSize = 8192;       // zero-padded
        double minPeakHz = 50;
        double maxPeakHz = 6000;
        double peakFloorDb = -40; // peaks further below the frame's highest are ignored
        double beta = 0.5;        // magnitude compression
        int nHarmonics = 10;
        double alpha = 0.85;      // harmonic weight decay
        double negWeight = 0.75;  // weight of the lobes at (h - 1/2) f0
        double tolBins = 5;       // half-width of the tolerance window, in 10-cent bins
        double minF0 = 60;
        double maxF0 = 1600;
        int maxCandidates = 5;
        int refineHarmonics = 5;  // harmonics fitted to refine a candidate's frequency
    };

    struct Peak {
        double freq;
        double amp;
    };

    struct Candidate {
        double freq;      // Hz
        double salience;  // > 0
    };

    /**
     * centreOffset is the number of samples from the start of a block
     * to the time of its frame, where its window is centred.
     */
    HarmonicSalience(float sampleRate, int blockSize, int stepSize,
                     int centreOffset, const Params &params);
    ~HarmonicSalience();

    void reset();

    /**
     * Takes one block as the host gives it, each block overlapping
     * the last by blockSize - stepSize samples.
     */
    void push(const float *block);

    /** Whether the next frame's window lies within the samples pushed. */
    bool ready() const;

    /**
     * Computes the salience of the next frame (as zeros past the
     * samples pushed, if it is not ready) and returns its candidates,
     * highest salience first.
     */
    std::vector<Candidate> next();

    /**
     * The highest salience of the last frame over the sum of its
     * peaks' compressed magnitudes: high for a harmonic sound, low
     * for noise or for many unrelated partials.
     */
    double harmonicity() const {
        return m_peakSum > 0 ? m_maxSalience / m_peakSum : 0.0;
    }

    /** The spectral peaks of the last frame, in order of frequency. */
    const std::vector<Peak> &peaks() const { return m_peaks; }

    /**
     * Lowers by penalty (0 to 1) the salience of each candidate of the
     * last frame that looks like the common fundamental of two voices
     * rather than a voice of its own, and sorts the candidates again.
     * See the static version.
     */
    void penaliseCommonFundamentals(std::vector<Candidate> &candidates,
                                    double penalty) const;

    /**
     * Two voices at a ratio of small integers p:q have a common
     * fundamental an octave (a fifth, 3:2) or a twelfth (a fourth or
     * a sixth) below the lower voice, whose harmonics p and q are the
     * two voices' own fundamentals: it collects both voices' partials
     * and can be the highest candidate. A single voice's candidates at
     * twice and three times its f0 are much weaker than its own (the
     * negative lobes see to that), so a candidate with two others at
     * distinct multiples of it (2 to 5), each with at least 0.7 of its
     * salience, is taken for a common fundamental, unless it has a
     * partial of its own: a peak among peaks, within 50 cents or
     * resolutionHz of it, no more than 20 dB below the highest. A
     * candidate taken for one loses penalty of its salience.
     */
    static void penaliseCommonFundamentals(std::vector<Candidate> &candidates,
                                           const std::vector<Peak> &peaks,
                                           double resolutionHz,
                                           double penalty);

    /**
     * The candidates as (frequency, probability) pairs for the pitch
     * HMM. They share voicedMass, the probability that the frame is
     * voiced, by their salience relative to the highest, to the fourth
     * power, which leaves a little to a second voice and next to
     * nothing to the sub-harmonics. Each frequency is refined by the
     * frame's peaks.
     */
    std::vector<std::pair<double, double> >
    candidateProbabilities(const std::vector<Candidate> &candidates,
                           double voicedMass) const;

    /**
     * Refines a candidate f0 by the least-squares fit of h f0 to the
     * frequencies of the spectral peaks within the tolerance of its
     * first nHarm harmonics. Returns f0 unchanged if none is found.
     */
    double refine(double f0, int nHarm) const;

private:
    float m_sampleRate;
    int m_blockSize;
    int m_stepSize;
    int m_centreOffset;
    Params m_p;

    std::vector<double> m_stream; // the most recent samples
    long long m_streamStart;      // absolute index of m_stream[0]
    long long m_streamEnd;        // absolute index one past the last sample
    long long m_pushed;           // blocks pushed
    long long m_frame;            // the next frame to compute

    std::vector<double> m_window;
    std::vector<double> m_fftIn;
    std::vector<double> m_fftOut;
    std::vector<double> m_mag;
    Vamp::FFTReal m_fft;

    double m_gridRef;   // Hz of grid bin 0
    int m_nBins;
    std::vector<double> m_sal;
    double m_maxSalience;
    double m_peakSum;
    std::vector<Peak> m_peaks;

    void append(const float *x, int n);
    long long centreOf(long long frame) const;
    void findPeaks();
    void accumulate();
};

#endif
