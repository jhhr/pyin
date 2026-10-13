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

#ifndef _LEVELGATE_H_
#define _LEVELGATE_H_

#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>

/**
 * A gate on voicing by level, relative to the singing's own loudness
 * rather than to a fixed level. Echo and reverb tails are quieter than
 * the voice they repeat, but they are as periodic; a quiet verse is
 * not quieter than itself.
 *
 * The reference is the 90th percentile of the levels of the sung
 * frames among the last windowFrames frames, the level the singing has
 * rather than its peaks. Only the sung frames count: the silence, the
 * breaths and the noise between phrases would otherwise pull the
 * reference down over a rest, and let the next phrase's echo through.
 * Where none of the window's frames is sung there is no reference, and
 * the gate is open. Only the past counts, so that a caller can emit its
 * results as it goes. A frame more than floorDb below the reference has
 * its gate faded from 1 to 0 over the next 6 dB.
 */
class LevelGate
{
public:
    LevelGate(double floorDb, int windowFrames) :
        m_floorDb(floorDb),
        m_windowFrames(std::max(1, windowFrames)),
        m_hist(nBins, 0),
        m_sung(0)
    { }

    void reset() {
        std::fill(m_hist.begin(), m_hist.end(), 0);
        m_bins.clear();
        m_sung = 0;
    }

    /**
     * Adds a frame of this RMS level, to the reference if it is sung,
     * and returns its gate, from 0 (not voiced) to 1 (as it was).
     */
    double process(float rms, bool sung) {
        double db = 20 * std::log10(std::max(double(rms), 1e-10));
        int bin = int(std::floor((db - bottomDb) * binsPerDb));
        bin = std::max(0, std::min(nBins - 1, bin));
        if (sung) {
            ++m_hist[bin];
            ++m_sung;
        }
        m_bins.push_back(sung ? bin : -1);
        while (int(m_bins.size()) > m_windowFrames) {
            if (m_bins.front() >= 0) {
                --m_hist[m_bins.front()];
                --m_sung;
            }
            m_bins.pop_front();
        }
        int n = m_sung;
        if (n == 0) return 1.0;
        int target = std::max(1, int(std::ceil(percentile * n)));
        int acc = 0, refBin = nBins - 1;
        for (int b = 0; b < nBins; ++b) {
            acc += m_hist[b];
            if (acc >= target) {
                refBin = b;
                break;
            }
        }
        double refDb = bottomDb + (refBin + 0.5) / binsPerDb;
        double floorDb = refDb - m_floorDb;
        double g = (db - floorDb + rampDb) / rampDb;
        return std::max(0.0, std::min(1.0, g));
    }

private:
    // The levels' histogram: quarter-dB bins from -120 dB to +10 dB
    static const int nBins = 520;
    static constexpr double binsPerDb = 4;
    static constexpr double bottomDb = -120;
    static constexpr float percentile = 0.9f;
    static constexpr double rampDb = 6;

    double m_floorDb;
    int m_windowFrames;
    std::vector<int> m_hist;  // the sung frames' levels
    std::deque<int> m_bins;   // each frame's bin, oldest first; -1 if not sung
    int m_sung;               // the sung frames in the window
};

#endif
