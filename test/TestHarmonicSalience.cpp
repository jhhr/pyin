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

#include "HarmonicSalience.h"

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdlib>
#include <vector>

using std::vector;
using std::pair;

typedef HarmonicSalience::Candidate Candidate;
typedef HarmonicSalience::Peak Peak;

static const int rate = 44100;
static const int blockSize = 2048;
static const int stepSize = 256;
static const double pi = 3.14159265358979323846;

static double cents(double a, double b)
{
    return 1200.0 * std::log2(a / b);
}

// A harmonic tone of f0 with partials of amplitude amp / k up to 8 kHz
static void addTone(vector<float> &audio, double f0, double amp)
{
    for (int k = 1; k * f0 < 8000; ++k) {
        double phase = 0.7 * k; // not all in phase, as in a voice
        for (size_t i = 0; i < audio.size(); ++i) {
            audio[i] += float(amp / k * std::sin(2 * pi * k * f0 * i / rate + phase));
        }
    }
}

// The candidates of the frame in the middle of the audio, which the
// salience is given block by block as a host gives them
static vector<Candidate> middleFrame(const vector<float> &audio,
                                     HarmonicSalience &hs)
{
    int frames = int(audio.size() - blockSize) / stepSize + 1;
    vector<Candidate> result;
    int taken = 0;
    for (int frame = 0; frame < frames; ++frame) {
        hs.push(&audio[frame * stepSize]);
        while (hs.ready()) {
            vector<Candidate> c = hs.next();
            if (taken++ == frames / 2) return c;
        }
    }
    BOOST_FAIL("no middle frame");
    return result;
}

static bool hasCandidateNear(const vector<Candidate> &c, double freq, size_t first,
                             double tolCents)
{
    for (size_t i = 0; i < c.size() && i < first; ++i) {
        if (std::fabs(cents(c[i].freq, freq)) < tolCents) return true;
    }
    return false;
}

BOOST_AUTO_TEST_SUITE(TestHarmonicSalience)

BOOST_AUTO_TEST_CASE(toneGivesItsF0)
{
    vector<float> audio(rate, 0.f);
    addTone(audio, 220.0, 0.3);
    HarmonicSalience hs(rate, blockSize, stepSize, blockSize / 4,
                        HarmonicSalience::Params());
    vector<Candidate> c = middleFrame(audio, hs);
    BOOST_REQUIRE(!c.empty());
    BOOST_CHECK_SMALL(cents(c[0].freq, 220.0), 5.0);
    // and refined by the partials, closer still
    BOOST_CHECK_SMALL(cents(hs.refine(c[0].freq, 5), 220.0), 1.0);
}

BOOST_AUTO_TEST_CASE(toneWithoutItsFundamentalGivesF0)
{
    // A low voice whose fundamental the separation (or a small speaker)
    // took away: its even partials alone are an octave higher tone's,
    // whose negative lobes, between them, fall on the odd ones
    vector<float> audio(rate, 0.f);
    addTone(audio, 220.0, 0.3);
    for (size_t i = 0; i < audio.size(); ++i) {
        audio[i] -= float(0.3 * std::sin(2 * pi * 220.0 * i / rate + 0.7));
    }
    HarmonicSalience hs(rate, blockSize, stepSize, blockSize / 4,
                        HarmonicSalience::Params());
    vector<Candidate> c = middleFrame(audio, hs);
    BOOST_REQUIRE(!c.empty());
    BOOST_CHECK_SMALL(cents(c[0].freq, 220.0), 5.0);
}

BOOST_AUTO_TEST_CASE(thirdGivesBothVoices)
{
    // A just major third, 4:5: YIN's deepest dip is at the period the
    // two have in common, two octaves below the lower
    vector<float> audio(rate, 0.f);
    addTone(audio, 330.0, 0.3);
    addTone(audio, 412.5, 0.3);
    HarmonicSalience hs(rate, blockSize, stepSize, blockSize / 4,
                        HarmonicSalience::Params());
    vector<Candidate> c = middleFrame(audio, hs);
    BOOST_CHECK(hasCandidateNear(c, 330.0, 3, 10));
    BOOST_CHECK(hasCandidateNear(c, 412.5, 3, 10));
    BOOST_CHECK(!hasCandidateNear(c, 82.5, 3, 50));
}

BOOST_AUTO_TEST_CASE(silenceGivesNothing)
{
    vector<float> audio(rate, 0.f);
    HarmonicSalience hs(rate, blockSize, stepSize, blockSize / 4,
                        HarmonicSalience::Params());
    vector<Candidate> c = middleFrame(audio, hs);
    BOOST_CHECK(c.empty());
    BOOST_CHECK_EQUAL(hs.harmonicity(), 0.0);
}

BOOST_AUTO_TEST_CASE(harmonicityHighForToneLowForNoise)
{
    // PYinVamp gives a frame no voiced probability for harmonicity
    // below 0.05 and its most at 0.12
    vector<float> tone(rate, 0.f);
    addTone(tone, 196.0, 0.3);
    HarmonicSalience hs(rate, blockSize, stepSize, blockSize / 4,
                        HarmonicSalience::Params());
    middleFrame(tone, hs);
    BOOST_CHECK_GT(hs.harmonicity(), 0.12);

    vector<float> noise(rate);
    srand(1);
    for (size_t i = 0; i < noise.size(); ++i) {
        noise[i] = 0.3f * (float(rand()) / float(RAND_MAX) - 0.5f);
    }
    HarmonicSalience hn(rate, blockSize, stepSize, blockSize / 4,
                        HarmonicSalience::Params());
    middleFrame(noise, hn);
    BOOST_CHECK_LT(hn.harmonicity(), 0.05);
}

BOOST_AUTO_TEST_CASE(windowCentredOnTheFrame)
{
    // A tone that starts at sample 22050. Frame f's window is centred
    // a quarter block after its block's start, at 256 f + 512, and
    // spans 256 f - 512 to 256 f + 1536: frame 80's ends before the
    // tone, frame 81's takes in its first samples, and frame 89's
    // starts after the tone has
    vector<float> audio(rate, 0.f);
    vector<float> tone(rate / 2, 0.f);
    addTone(tone, 220.0, 0.3);
    std::copy(tone.begin(), tone.end(), audio.begin() + rate / 2);
    HarmonicSalience hs(rate, blockSize, stepSize, blockSize / 4,
                        HarmonicSalience::Params());
    int frames = int(audio.size() - blockSize) / stepSize + 1;
    int taken = 0;
    for (int frame = 0; frame < frames; ++frame) {
        hs.push(&audio[frame * stepSize]);
        while (hs.ready()) {
            vector<Candidate> c = hs.next();
            if (taken == 80) BOOST_CHECK(c.empty());
            if (taken == 81) BOOST_CHECK(!c.empty());
            if (taken == 89) BOOST_CHECK(hasCandidateNear(c, 220.0, 1, 5));
            ++taken;
        }
    }
    BOOST_CHECK_GT(taken, 89);
}

BOOST_AUTO_TEST_CASE(framesWaitForTheirWindow)
{
    // Blocks shorter than the window: a frame is ready only once the
    // blocks after it have brought the end of its window, and the rest
    // come at the end, so no frame is lost
    const int shortBlock = 1024;
    vector<float> audio(rate, 0.f);
    addTone(audio, 220.0, 0.3);
    HarmonicSalience hs(rate, shortBlock, stepSize, shortBlock / 4,
                        HarmonicSalience::Params());
    int blocks = int(audio.size() - shortBlock) / stepSize + 1;
    int taken = 0;
    for (int b = 0; b < blocks; ++b) {
        hs.push(&audio[b * stepSize]);
        while (hs.ready()) {
            // frame taken: its window, 1024 past its centre, which is
            // 256 into its block, is within the blocks pushed
            BOOST_CHECK_LE(taken * stepSize + 256 + 1024, b * stepSize + shortBlock);
            vector<Candidate> c = hs.next();
            if (taken > 4 && taken < blocks - 8) {
                BOOST_CHECK(hasCandidateNear(c, 220.0, 1, 5));
            }
            ++taken;
        }
    }
    BOOST_CHECK_LT(taken, blocks);
    while (taken < blocks) {
        hs.next();
        ++taken;
    }
    BOOST_CHECK(!hs.ready());
}

BOOST_AUTO_TEST_CASE(candidateProbabilitiesFollowSalience)
{
    // Before any frame there are no peaks to refine the frequencies by
    HarmonicSalience hs(rate, blockSize, stepSize, blockSize / 4,
                        HarmonicSalience::Params());
    vector<Candidate> c { { 200.0, 1.0 }, { 300.0, 0.5 }, { 400.0, 0.25 } };
    vector<pair<double, double> > p = hs.candidateProbabilities(c, 0.6);
    BOOST_REQUIRE_EQUAL(p.size(), 3u);
    double w[3] = { 1.0, std::pow(0.5, 4), std::pow(0.25, 4) };
    double wsum = w[0] + w[1] + w[2];
    double sum = 0;
    for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_EQUAL(p[i].first, c[i].freq);
        BOOST_CHECK_CLOSE(p[i].second, 0.6 * w[i] / wsum, 1e-9);
        sum += p[i].second;
    }
    BOOST_CHECK_CLOSE(sum, 0.6, 1e-9);

    // no voiced probability, or no candidates: an unvoiced frame
    BOOST_CHECK(hs.candidateProbabilities(c, 0.0).empty());
    BOOST_CHECK(hs.candidateProbabilities(vector<Candidate>(), 0.6).empty());
}

// The spectral peaks of two voices at 330 and 495 Hz, a fifth apart:
// their harmonics, and none at their common fundamental, 165 Hz
static vector<Peak> fifthPeaks()
{
    vector<Peak> peaks;
    double freqs[] = { 330, 495, 660, 990, 1320, 1485 };
    for (double f : freqs) peaks.push_back(Peak { f, 1.0 });
    return peaks;
}

BOOST_AUTO_TEST_CASE(commonFundamentalPenalised)
{
    vector<Candidate> c { { 165.0, 1.0 }, { 495.0, 0.9 }, { 330.0, 0.8 } };
    HarmonicSalience::penaliseCommonFundamentals(c, fifthPeaks(), 21.5, 0.8);
    BOOST_REQUIRE_EQUAL(c.size(), 3u);
    BOOST_CHECK_EQUAL(c[0].freq, 495.0);
    BOOST_CHECK_EQUAL(c[1].freq, 330.0);
    BOOST_CHECK_EQUAL(c[2].freq, 165.0);
    BOOST_CHECK_CLOSE(c[2].salience, 0.2, 1e-9);
}

BOOST_AUTO_TEST_CASE(commonFundamentalWithAPartialSpared)
{
    // A low voice with its own fundamental (here 20 dB below its
    // strongest partial, the most that still counts)
    vector<Peak> peaks = fifthPeaks();
    peaks.insert(peaks.begin(), Peak { 165.0, 0.1 });
    vector<Candidate> c { { 165.0, 1.0 }, { 495.0, 0.9 }, { 330.0, 0.8 } };
    HarmonicSalience::penaliseCommonFundamentals(c, peaks, 21.5, 0.8);
    BOOST_CHECK_EQUAL(c[0].freq, 165.0);
    BOOST_CHECK_EQUAL(c[0].salience, 1.0);
}

BOOST_AUTO_TEST_CASE(singleVoiceNotPenalised)
{
    // One voice's candidates at twice and three times its f0 are weak,
    // even without a partial of its own at f0
    vector<Peak> peaks { { 400.0, 1.0 }, { 600.0, 1.0 }, { 800.0, 1.0 } };
    vector<Candidate> c { { 200.0, 1.0 }, { 400.0, 0.3 }, { 600.0, 0.3 } };
    HarmonicSalience::penaliseCommonFundamentals(c, peaks, 21.5, 0.8);
    BOOST_CHECK_EQUAL(c[0].freq, 200.0);
    BOOST_CHECK_EQUAL(c[0].salience, 1.0);
    BOOST_CHECK_EQUAL(c[1].salience, 0.3);
    BOOST_CHECK_EQUAL(c[2].salience, 0.3);
}

BOOST_AUTO_TEST_SUITE_END()
