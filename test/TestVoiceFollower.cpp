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

// The voice follower on candidates made up frame by frame: which
// voice its second path takes, where its gate lets that through, and
// that the first path's entries all come out, in order, as the input
// goes on.

#include "VoiceFollower.h"

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <vector>

using std::vector;
using std::pair;

namespace {

const int lag = 100;
const double lowerHz = 220.0;
const double upperHz = 220.0 * std::pow(2.0, 4 / 12.0); // a major third up

double midi(double hz) { return 12 * std::log2(hz / 440.0) + 69; }

double cents(double a, double b) { return 1200 * std::log2(a / b); }

// One frame: its voices (Hz), each with the probability it gets and a
// likelihood of 0.9 of being a voice, and the first path's frequency
struct Frame {
    vector<pair<double, double> > voices; // Hz, probability
    float first;
};

// Feeds the frames as PYinVamp does, the first path's entries coming
// with the fixed lag (frame 0 none, frame n - lag twice), and takes the
// entries as they come out. Each entry is paired with the number of
// frames given when it came out (-1: after the end of the input)
struct Run {
    vector<VoiceFollower::Entry> entries;
    vector<int> givenWhenOut;
};

Run run(VoiceFollower::Voice voice, const vector<Frame> &frames)
{
    VoiceFollower vf(voice, lag, 0.2, 3);
    Run r;
    auto take = [&](int given) {
        for (const VoiceFollower::Entry &e : vf.takeEntries()) {
            r.entries.push_back(e);
            r.givenWhenOut.push_back(given);
        }
    };
    auto first = [&](int frame) {
        VoiceFollower::Entry e;
        e.frame = frame;
        e.freq = frames[frame].first;
        e.timestamp = Vamp::RealTime::frame2RealTime(frame * 256, 44100);
        vf.addFirst(e);
    };
    for (size_t i = 0; i < frames.size(); ++i) {
        vector<pair<double, double> > pitchProb;
        vector<double> likelihoods;
        vector<VoiceFollower::GateVoice> gateVoices;
        for (const auto &v : frames[i].voices) {
            pitchProb.push_back(pair<double, double>(midi(v.first), v.second));
            likelihoods.push_back(0.9);
            VoiceFollower::GateVoice gv;
            gv.pitch = float(midi(v.first));
            gv.likelihood = 0.9f;
            gateVoices.push_back(gv);
        }
        vf.addFrame(pitchProb, likelihoods, gateVoices, 1.0);
        int n = vf.frameCount();
        if (n >= lag + 1) first(n - lag);
        take(n);
    }
    int n = vf.frameCount();
    for (int k = std::max(0, n - lag); k < n; ++k) first(k);
    vf.finish();
    take(-1);
    return r;
}

// Frames a to b: the given voices, the first path on firstHz
void fill(vector<Frame> &frames, int a, int b,
          const vector<pair<double, double> > &voices, double firstHz)
{
    for (int i = a; i < b; ++i) {
        frames[i].voices = voices;
        frames[i].first = float(firstHz);
    }
}

const vector<pair<double, double> > both { { lowerHz, 0.3 }, { upperHz, 0.6 } };
const vector<pair<double, double> > upperOnly { { upperHz, 0.9 } };

// The entries' frequencies by frame (the last entry for a frame twice
// given)
vector<float> byFrame(const Run &r, int n)
{
    vector<float> f(n, 0.f);
    for (const auto &e : r.entries) f[e.frame] = e.freq;
    return f;
}

}

BOOST_AUTO_TEST_SUITE(TestVoiceFollower)

BOOST_AUTO_TEST_CASE(oneVoiceLeavesTheFirstPath)
{
    vector<Frame> frames(1200);
    fill(frames, 0, 1200, upperOnly, upperHz);
    for (VoiceFollower::Voice v : { VoiceFollower::Lower, VoiceFollower::Higher,
                                    VoiceFollower::Middle }) {
        Run r = run(v, frames);
        for (const auto &e : r.entries) {
            BOOST_CHECK_EQUAL(e.freq, frames[e.frame].first);
        }
    }
}

BOOST_AUTO_TEST_CASE(lowerAndHigherInsideAPassage)
{
    // Two voices a third apart all through, the upper the stronger: the
    // first path on the voice the follower is not to take
    vector<Frame> frames(1200);
    fill(frames, 0, 1200, both, upperHz);
    vector<float> lower = byFrame(run(VoiceFollower::Lower, frames), 1200);
    fill(frames, 0, 1200, both, lowerHz);
    vector<float> higher = byFrame(run(VoiceFollower::Higher, frames), 1200);
    for (int k = 200; k < 1000; ++k) {
        BOOST_TEST_CONTEXT("frame " << k) {
            BOOST_CHECK_LT(std::fabs(cents(lower[k], lowerHz)), 50);
            BOOST_CHECK_LT(std::fabs(cents(higher[k], upperHz)), 50);
        }
    }
}

BOOST_AUTO_TEST_CASE(voicingIsTheFirstPaths)
{
    // The first path unvoiced for a moment inside the passage, and the
    // second path voiced there: the entries stay unvoiced
    vector<Frame> frames(1200);
    fill(frames, 0, 1200, both, upperHz);
    for (int k = 500; k < 520; ++k) frames[k].first = -float(upperHz);
    Run r = run(VoiceFollower::Lower, frames);
    for (const auto &e : r.entries) {
        if (e.frame >= 500 && e.frame < 520) {
            BOOST_CHECK_EQUAL(e.freq, frames[e.frame].first);
        } else if (e.frame >= 200 && e.frame < 1000) {
            BOOST_CHECK_GT(e.freq, 0.f);
        }
    }
}

BOOST_AUTO_TEST_CASE(gateReachesAboutASecondPastAPassage)
{
    // Two voices for 600 frames, then the upper alone, the first path a
    // little off it, so that where the gate lets the second path's
    // pitch through, the entry shows it. A centred window of 344 frames
    // holds 80% two-voice frames up to frame 496, and the passage is
    // widened by 172 frames: to frame 668
    const double offHz = upperHz * 1.01;
    vector<Frame> frames(1400);
    fill(frames, 0, 600, both, offHz);
    fill(frames, 600, 1400, upperOnly, offHz);
    vector<float> f = byFrame(run(VoiceFollower::Lower, frames), 1400);
    for (int k = 200; k < 590; ++k) {
        BOOST_TEST_CONTEXT("frame " << k) {
            BOOST_CHECK_LT(std::fabs(cents(f[k], lowerHz)), 50);
        }
    }
    for (int k = 610; k <= 668; ++k) {
        BOOST_TEST_CONTEXT("frame " << k) {
            BOOST_CHECK_LT(std::fabs(cents(f[k], upperHz)), 5);
        }
    }
    for (int k = 669; k < 1400; ++k) {
        BOOST_TEST_CONTEXT("frame " << k) {
            BOOST_CHECK_EQUAL(f[k], float(offHz));
        }
    }
}

BOOST_AUTO_TEST_CASE(gateWantsFourFifthsOfTheVoicedFrames)
{
    // The second voice in five frames of every six: a passage. In three
    // of every four: none
    for (int every : { 6, 4 }) {
        vector<Frame> frames(1200);
        fill(frames, 0, 1200, both, upperHz);
        for (int k = 0; k < 1200; k += every) frames[k].voices = upperOnly;
        vector<float> f = byFrame(run(VoiceFollower::Lower, frames), 1200);
        int lower = 0;
        for (int k = 300; k < 900; ++k) {
            if (std::fabs(cents(f[k], lowerHz)) < 50) ++lower;
        }
        BOOST_TEST_CONTEXT("one frame in " << every << " with one voice") {
            if (every == 6) BOOST_CHECK_GE(lower, 600 * 9 / 10);
            else BOOST_CHECK_EQUAL(lower, 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(gateWantsAFifthVoiced)
{
    // Short notes, one frame in six voiced: the frames that are voiced
    // all have two voices, but a fifth of the window is not voiced
    vector<Frame> frames(1200);
    for (int k = 0; k < 1200; ++k) {
        frames[k].first = -float(upperHz);
        if (k % 6 == 0) {
            frames[k].voices = both;
            frames[k].first = float(upperHz);
        }
    }
    Run r = run(VoiceFollower::Lower, frames);
    for (const auto &e : r.entries) {
        BOOST_CHECK_EQUAL(e.freq, frames[e.frame].first);
    }
}

BOOST_AUTO_TEST_CASE(entriesComeOutAsTheInputGoesOn)
{
    // Every entry the first path gives comes out once, in order, with
    // its timestamp: during the input, no more than 443 frames after its
    // own frame, but for those the end of the input lets out
    const int n = 3000;
    vector<Frame> frames(n);
    fill(frames, 0, 1000, upperOnly, upperHz);
    fill(frames, 1000, 2000, both, upperHz);
    fill(frames, 2000, n, upperOnly, upperHz);
    Run r = run(VoiceFollower::Lower, frames);

    // frames 1 to n - lag, then n - lag again and the rest of the lag
    BOOST_REQUIRE_EQUAL(int(r.entries.size()), n);
    int during = 0;
    for (size_t i = 0; i < r.entries.size(); ++i) {
        const auto &e = r.entries[i];
        int expected = int(i) < n - lag ? int(i) + 1 : int(i);
        BOOST_CHECK_EQUAL(e.frame, expected);
        BOOST_CHECK(e.timestamp ==
                    Vamp::RealTime::frame2RealTime(e.frame * 256, 44100));
        if (r.givenWhenOut[i] >= 0) {
            ++during;
            BOOST_CHECK_LE(r.givenWhenOut[i] - 1 - e.frame, 443);
        }
    }
    BOOST_CHECK_GE(during, n - 500);
}

BOOST_AUTO_TEST_SUITE_END()
