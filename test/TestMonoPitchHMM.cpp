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

#include "MonoPitchHMM.h"

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <vector>

using std::vector;
using std::pair;

// A candidate as PYinVamp gives MonoPitchHMM one: (MIDI pitch, probability)
static pair<double, double> candidate(double freq, double prob)
{
    return pair<double, double>(12 * std::log2(freq / 440) + 69, prob);
}

BOOST_AUTO_TEST_SUITE(TestMonoPitchHMM)

BOOST_AUTO_TEST_CASE(voicedFloorWithoutCandidates)
{
    vector<pair<double, double> > none;

    MonoPitchHMM plain(0);
    vector<double> obs = plain.calculateObsProb(none);
    for (int i = 0; i < plain.m_nPitch; ++i) {
        BOOST_CHECK_EQUAL(obs[i], 0.0);
    }

    MonoPitchHMM floored(0, 0.2);
    obs = floored.calculateObsProb(none);
    for (int i = 0; i < floored.m_nPitch; ++i) {
        BOOST_CHECK_CLOSE(obs[i], 0.2 * obs[i + floored.m_nPitch], 1e-9);
        BOOST_CHECK_GT(obs[i], 0.0);
    }
}

BOOST_AUTO_TEST_CASE(voicedFloorKeepsCandidates)
{
    // A candidate's state keeps its observation; the others get the
    // floor, a fifth of an unvoiced state's
    vector<pair<double, double> > one { candidate(220, 0.8) };
    MonoPitchHMM plain(0);
    MonoPitchHMM floored(0, 0.2);
    vector<double> p = plain.calculateObsProb(one);
    vector<double> f = floored.calculateObsProb(one);
    int n = plain.m_nPitch;
    int withCandidate = 0;
    for (int i = 0; i < n; ++i) {
        if (p[i] > 0) {
            ++withCandidate;
            BOOST_CHECK_EQUAL(f[i], p[i]);
        } else {
            BOOST_CHECK_CLOSE(f[i], 0.2 * f[i + n], 1e-9);
        }
        BOOST_CHECK_EQUAL(f[i + n], p[i + n]);
    }
    BOOST_CHECK_EQUAL(withCandidate, 1);
}

BOOST_AUTO_TEST_CASE(voicedFloorKeepsStateFrequency)
{
    // A voiced state the path reached over the floor, with no candidate
    // within a semitone, keeps its own frequency
    vector<pair<double, double> > far { candidate(440, 0.8) };
    MonoPitchHMM plain(0);
    MonoPitchHMM floored(0, 0.2);
    int state = 0;
    while (plain.m_freqs[state] < 220) ++state;
    BOOST_CHECK_CLOSE(plain.nearestFreq(state, far), 440.0, 1e-4);
    BOOST_CHECK_CLOSE(floored.nearestFreq(state, far), plain.m_freqs[state], 1e-4);

    // a candidate near it still gives its own frequency
    vector<pair<double, double> > near { candidate(plain.m_freqs[state] * 1.01, 0.8) };
    BOOST_CHECK_CLOSE(floored.nearestFreq(state, near), plain.m_freqs[state] * 1.01, 1e-4);
}

BOOST_AUTO_TEST_CASE(unvoicedHoldDividesUnvoicedMoves)
{
    MonoPitchHMM plain(0);
    MonoPitchHMM held(0, 0, 3);
    BOOST_REQUIRE_EQUAL(plain.m_transProb.size(), held.m_transProb.size());
    size_t n = size_t(plain.m_nPitch);
    int moves = 0;
    for (size_t t = 0; t < plain.m_transProb.size(); ++t) {
        BOOST_REQUIRE_EQUAL(plain.m_from[t], held.m_from[t]);
        BOOST_REQUIRE_EQUAL(plain.m_to[t], held.m_to[t]);
        size_t from = plain.m_from[t], to = plain.m_to[t];
        if (from >= n && to >= n && from != to) {
            // an unvoiced state moving to another pitch
            BOOST_CHECK_CLOSE(held.m_transProb[t], plain.m_transProb[t] / 4, 1e-9);
            ++moves;
        } else {
            // staying, voicing, unvoicing and every voiced transition
            BOOST_CHECK_EQUAL(held.m_transProb[t], plain.m_transProb[t]);
        }
    }
    BOOST_CHECK_GT(moves, 0);
}

BOOST_AUTO_TEST_SUITE_END()
