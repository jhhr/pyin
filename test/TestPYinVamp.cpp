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

// The plugin's parameters for voices with harmonies, doubles and
// echoes, end to end on synthetic voices that pYIN 1.2 gets an octave
// or two wrong.

#include "PYinVamp.h"

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::string;
using std::vector;

static const int rate = 44100;
static const double pi = 3.14159265358979323846;

typedef map<string, float> Params;

// Tony's settings for a song's analysis
static Params robust()
{
    Params p;
    p["lowampsuppression"] = 0.2f;
    p["onsetsensitivity"] = 0.7f;
    p["prunethresh"] = 0.1f;
    p["salmode"] = 1;
    p["candidaterule"] = 1;
    p["commonf0"] = 0.8f;
    p["voicehold"] = 3;
    p["voicedfloor"] = 0.2f;
    p["levelfloor"] = 13;
    return p;
}

// A harmonic tone of f0 with partials of amplitude amp / k up to 8 kHz,
// and with half, partials between them at (k + 1/2) f0, as of a voice
// whose period doubles (vocal fry, a rough low note)
static vector<float> tone(double f0, double amp, double seconds, double half = 0)
{
    vector<float> audio(size_t(seconds * rate), 0.f);
    for (int k = 1; k * f0 < 8000; ++k) {
        for (size_t i = 0; i < audio.size(); ++i) {
            audio[i] += float(amp / k * std::sin(2 * pi * k * f0 * i / rate + 0.7 * k));
        }
    }
    for (int k = 0; half > 0 && (k + 0.5) * f0 < 8000; ++k) {
        for (size_t i = 0; i < audio.size(); ++i) {
            audio[i] += float(half * amp / (k + 0.5) *
                              std::sin(2 * pi * (k + 0.5) * f0 * i / rate + 0.3 * k));
        }
    }
    return audio;
}

static vector<float> mix(const vector<float> &a, const vector<float> &b)
{
    vector<float> m(a);
    for (size_t i = 0; i < m.size() && i < b.size(); ++i) m[i] += b[i];
    return m;
}

static double seconds(const Vamp::RealTime &t)
{
    return t.sec + t.nsec / 1e9;
}

typedef map<string, Vamp::Plugin::FeatureList> Outputs;

// Every output's features, by the output's identifier, as a host gets
// them: blocks of the preferred size from the first frame, the last
// ones padded with zeros
static Outputs run(const vector<float> &audio, const Params &params)
{
    PYinVamp plugin(rate);
    for (const auto &p : params) plugin.setParameter(p.first, p.second);
    size_t blockSize = plugin.getPreferredBlockSize();
    size_t stepSize = plugin.getPreferredStepSize();
    BOOST_REQUIRE(plugin.initialise(1, stepSize, blockSize));

    Vamp::Plugin::OutputList list = plugin.getOutputDescriptors();
    Outputs outputs;
    auto collect = [&](Vamp::Plugin::FeatureSet &fs) {
        for (auto &o : fs) {
            Vamp::Plugin::FeatureList &to = outputs[list.at(o.first).identifier];
            to.insert(to.end(), o.second.begin(), o.second.end());
        }
    };
    vector<float> block(blockSize);
    for (size_t start = 0; start < audio.size(); start += stepSize) {
        for (size_t i = 0; i < blockSize; ++i) {
            block[i] = start + i < audio.size() ? audio[start + i] : 0.f;
        }
        const float *buffers[1] = { &block[0] };
        Vamp::Plugin::FeatureSet fs = plugin.process
            (buffers, Vamp::RealTime::frame2RealTime(long(start), rate));
        collect(fs);
    }
    Vamp::Plugin::FeatureSet rest = plugin.getRemainingFeatures();
    collect(rest);
    return outputs;
}

// The smoothed pitch track as (time, f0): the voiced frames, or every
// frame with outputunvoiced 2 (the unvoiced ones negative)
static vector<pair<double, double> > pitchTrack(const vector<float> &audio,
                                                const Params &params)
{
    Outputs outputs = run(audio, params);
    vector<pair<double, double> > track;
    for (const auto &f : outputs["smoothedpitchtrack"]) {
        if (!f.values.empty()) {
            track.push_back(pair<double, double>(seconds(f.timestamp), f.values[0]));
        }
    }
    return track;
}

// The probability the pitch HMM is given, on average over the frames
// from 0.25 to 1.25 s, for candidates within 50 cents of freq
static double meanProbabilityNear(const vector<float> &audio, const Params &params,
                                  double freq)
{
    Outputs outputs = run(audio, params);
    const Vamp::Plugin::FeatureList &candidates = outputs["f0candidates"];
    const Vamp::Plugin::FeatureList &probs = outputs["f0probs"];
    BOOST_REQUIRE_EQUAL(candidates.size(), probs.size());
    double sum = 0;
    int frames = 0;
    for (size_t i = 0; i < candidates.size(); ++i) {
        double t = seconds(candidates[i].timestamp);
        if (t < 0.25 || t > 1.25) continue;
        ++frames;
        BOOST_REQUIRE_EQUAL(candidates[i].values.size(), probs[i].values.size());
        for (size_t j = 0; j < candidates[i].values.size(); ++j) {
            if (std::fabs(1200 * std::log2(candidates[i].values[j] / freq)) < 50) {
                sum += probs[i].values[j];
            }
        }
    }
    BOOST_REQUIRE(frames > 0);
    return sum / frames;
}

struct Count {
    int voiced = 0;   // voiced frames in the span
    int onVoice = 0;  // ... within 50 cents of one of the voices
};

static Count count(const vector<pair<double, double> > &track, double from, double to,
                   const vector<double> &voices)
{
    Count c;
    for (const auto &p : track) {
        if (p.first < from || p.first > to || p.second <= 0) continue;
        ++c.voiced;
        for (double v : voices) {
            if (std::fabs(1200 * std::log2(p.second / v)) < 50) {
                ++c.onVoice;
                break;
            }
        }
    }
    return c;
}

// The frames of a second of audio: every one from 0.25 to 1.25 s
static const int framesInSpan = int(1.0 * rate / 256);

BOOST_AUTO_TEST_SUITE(TestPYinVamp)

BOOST_AUTO_TEST_CASE(parametersDefaultToPyin12)
{
    PYinVamp plugin(rate);
    const char *ids[] = { "candidaterule", "salmode", "commonf0",
                          "levelfloor", "voicedfloor", "voicehold",
                          "voicefollow" };
    Vamp::Plugin::ParameterList list = plugin.getParameterDescriptors();
    for (const char *id : ids) {
        BOOST_TEST_CONTEXT(id) {
            bool found = false;
            for (const auto &d : list) {
                if (d.identifier != id) continue;
                found = true;
                BOOST_CHECK_EQUAL(d.defaultValue, 0.f);
                BOOST_CHECK(d.minValue <= 0.f && d.maxValue > 0.f);
            }
            BOOST_CHECK(found);
            BOOST_CHECK_EQUAL(plugin.getParameter(id), 0.f);
            plugin.setParameter(id, 1.f);
            BOOST_CHECK_EQUAL(plugin.getParameter(id), 1.f);
        }
    }
}

BOOST_AUTO_TEST_CASE(zeroIsThePlainAnalysis)
{
    // Setting the parameters to 0 changes nothing at all
    vector<float> audio = mix(tone(330, 0.25, 1.5), tone(495, 0.25, 1.5));
    Params zero;
    const char *ids[] = { "candidaterule", "salmode", "commonf0",
                          "levelfloor", "voicedfloor", "voicehold",
                          "voicefollow" };
    for (const char *id : ids) zero[id] = 0.f;
    vector<pair<double, double> > plain = pitchTrack(audio, Params());
    vector<pair<double, double> > zeroed = pitchTrack(audio, zero);
    BOOST_REQUIRE(!plain.empty());
    BOOST_CHECK(plain == zeroed);
}

BOOST_AUTO_TEST_CASE(everyFrameStillHasAnEstimate)
{
    // The salience's frames come out as YIN's do, one for each block
    // (and in fixed-lag mode one twice, see docs/forks.md in Tony)
    vector<float> audio = mix(tone(330, 0.25, 1.5), tone(495, 0.25, 1.5));
    Params plain;
    plain["outputunvoiced"] = 2;
    Params params = robust();
    params["outputunvoiced"] = 2;
    vector<pair<double, double> > a = pitchTrack(audio, plain);
    vector<pair<double, double> > b = pitchTrack(audio, params);
    BOOST_REQUIRE_GE(a.size(), (audio.size() + 255) / 256);
    BOOST_REQUIRE_EQUAL(a.size(), b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        BOOST_CHECK_EQUAL(a[i].first, b[i].first);
    }
}

BOOST_AUTO_TEST_CASE(fifthTracksAVoiceNotTheCommonFundamental)
{
    // E4 and B4 a just fifth apart, 2:3: pYIN 1.2 tracks E3, their
    // common fundamental, and E2
    vector<float> audio = mix(tone(330, 0.25, 1.5), tone(495, 0.25, 1.5));
    Count c = count(pitchTrack(audio, robust()), 0.25, 1.25, { 330, 495 });
    BOOST_CHECK_GE(c.voiced, framesInSpan * 9 / 10);
    BOOST_CHECK_EQUAL(c.onVoice, c.voiced);
}

BOOST_AUTO_TEST_CASE(commonFundamentalLosesItsProbability)
{
    // In the same fifth the salience's candidate at E3, the common
    // fundamental, has a quarter of each frame's probability; the
    // penalty, which finds E4 and B4 at twice and three times it and no
    // partial at E3 itself, takes nearly all of it away
    vector<float> audio = mix(tone(330, 0.25, 1.5), tone(495, 0.25, 1.5));
    Params without = robust();
    without["commonf0"] = 0;
    double penalised = meanProbabilityNear(audio, robust(), 165);
    double unpenalised = meanProbabilityNear(audio, without, 165);
    BOOST_CHECK_GT(unpenalised, 0.1);
    BOOST_CHECK_LT(penalised, unpenalised / 50);
}

BOOST_AUTO_TEST_CASE(thirdTracksAVoiceNotTheCommonPeriod)
{
    // E4 and G#4 a just major third apart, 4:5: pYIN 1.2 tracks E2, two
    // octaves below the lower voice
    vector<float> audio = mix(tone(330, 0.25, 1.5), tone(412.5, 0.25, 1.5));
    Count c = count(pitchTrack(audio, robust()), 0.25, 1.25, { 330, 412.5 });
    BOOST_CHECK_GE(c.voiced, framesInSpan * 9 / 10);
    BOOST_CHECK_EQUAL(c.onVoice, c.voiced);
}

BOOST_AUTO_TEST_CASE(periodDoubledToneKeepsItsPitch)
{
    // D3 with partials between its harmonics 20 dB down: its waveform
    // repeats every two periods, and pYIN 1.2 tracks D2
    vector<float> audio = tone(147, 0.3, 1.5, 0.1);
    Count c = count(pitchTrack(audio, robust()), 0.25, 1.25, { 147 });
    BOOST_CHECK_GE(c.voiced, framesInSpan * 9 / 10);
    BOOST_CHECK_EQUAL(c.onVoice, c.voiced);
}

BOOST_AUTO_TEST_CASE(echoTailIsNotVoiced)
{
    // A note and then its echo 20 dB down: pYIN 1.2 tracks the echo as
    // if it were sung
    vector<float> audio = tone(220, 0.3, 0.8);
    vector<float> echo = tone(220, 0.03, 0.6);
    audio.insert(audio.end(), echo.begin(), echo.end());
    vector<pair<double, double> > track = pitchTrack(audio, robust());
    Count note = count(track, 0.1, 0.75, { 220 });
    BOOST_CHECK_GE(note.voiced, int(0.65 * rate / 256) * 9 / 10);
    BOOST_CHECK_EQUAL(note.onVoice, note.voiced);
    BOOST_CHECK_EQUAL(count(track, 0.9, 1.35, { 220 }).voiced, 0);
}

// A harmonic tone of f0 at one level for the first half and another
// for the second
static vector<float> toneInTwoHalves(double f0, double seconds,
                                     double firstAmp, double secondAmp)
{
    vector<float> audio = tone(f0, 1.0, seconds);
    for (size_t i = 0; i < audio.size(); ++i) {
        audio[i] *= float(i < audio.size() / 2 ? firstAmp : secondAmp);
    }
    return audio;
}

// Two voices a just major third apart for four seconds, the upper the
// louder for the first two and the lower for the last two: the
// strongest voice changes at 2 s
static vector<float> thirdWithLevelsCrossing()
{
    return mix(toneInTwoHalves(220, 4.0, 0.12, 0.3),
               toneInTwoHalves(275, 4.0, 0.3, 0.12));
}

static Params following(int voice)
{
    Params p = robust();
    p["voicefollow"] = float(voice);
    return p;
}

BOOST_AUTO_TEST_CASE(voiceFollowZeroIsTheStrongest)
{
    // 0 is the analysis without it, pitch track and notes alike
    vector<float> audio = thirdWithLevelsCrossing();
    Outputs without = run(audio, robust());
    Outputs zero = run(audio, following(0));
    for (const char *output : { "smoothedpitchtrack", "notes" }) {
        BOOST_TEST_CONTEXT(output) {
            const Vamp::Plugin::FeatureList &a = without[output], &b = zero[output];
            BOOST_REQUIRE(!a.empty());
            BOOST_REQUIRE_EQUAL(a.size(), b.size());
            for (size_t i = 0; i < a.size(); ++i) {
                BOOST_CHECK(a[i].timestamp == b[i].timestamp);
                BOOST_CHECK(a[i].values == b[i].values);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(lowerAndHigherVoiceAcrossTheCrossing)
{
    // The strongest changes voice where the levels cross; the lower and
    // the higher hold theirs through it, and are voiced where it is
    vector<float> audio = thirdWithLevelsCrossing();
    vector<pair<double, double> > strongest = pitchTrack(audio, following(0));
    vector<pair<double, double> > lower = pitchTrack(audio, following(1));
    vector<pair<double, double> > higher = pitchTrack(audio, following(2));

    Count before = count(strongest, 0.5, 1.5, { 275 });
    Count after = count(strongest, 2.5, 3.5, { 220 });
    BOOST_CHECK_GE(before.onVoice, framesInSpan * 9 / 10);
    BOOST_CHECK_GE(after.onVoice, framesInSpan * 9 / 10);

    Count low = count(lower, 0.5, 3.5, { 220 });
    Count high = count(higher, 0.5, 3.5, { 275 });
    BOOST_CHECK_GE(low.onVoice, 3 * framesInSpan * 9 / 10);
    BOOST_CHECK_GE(high.onVoice, 3 * framesInSpan * 9 / 10);

    BOOST_REQUIRE_EQUAL(lower.size(), strongest.size());
    BOOST_REQUIRE_EQUAL(higher.size(), strongest.size());
    for (size_t i = 0; i < strongest.size(); ++i) {
        BOOST_CHECK_EQUAL(lower[i].first, strongest[i].first);
        BOOST_CHECK_EQUAL(higher[i].first, strongest[i].first);
    }
}

BOOST_AUTO_TEST_CASE(soloWithAnEchoIsUnchanged)
{
    // A note, the next a third up, and the first one's echo 12 dB down
    // under the start of the second: below it, so "the lower voice" if
    // anything were, but not long enough together for a passage of two
    // voices. Every voice to follow gives the strongest's track
    vector<float> audio = tone(220, 0.3, 1.0);
    vector<float> second = tone(262, 0.3, 1.0);
    audio.insert(audio.end(), second.begin(), second.end());
    vector<float> echo = tone(220, 0.075, 0.4);
    for (size_t i = 0; i < echo.size(); ++i) audio[rate + i] += echo[i];
    vector<pair<double, double> > strongest = pitchTrack(audio, following(0));
    BOOST_REQUIRE(!strongest.empty());
    for (int voice = 1; voice <= 3; ++voice) {
        BOOST_TEST_CONTEXT("voicefollow " << voice) {
            BOOST_CHECK(pitchTrack(audio, following(voice)) == strongest);
        }
    }
}

BOOST_AUTO_TEST_CASE(trackComesOutDuringTheAnalysis)
{
    // With a voice to follow, the track still comes out of process()
    // as the input goes on, a little later than without one: at most
    // 443 frames after its own (2.6 s), against 100
    vector<float> audio = mix(tone(220, 0.2, 8.0), tone(275, 0.3, 8.0));
    PYinVamp plugin(rate);
    for (const auto &p : following(1)) plugin.setParameter(p.first, p.second);
    size_t blockSize = plugin.getPreferredBlockSize();
    size_t stepSize = plugin.getPreferredStepSize();
    BOOST_REQUIRE(plugin.initialise(1, stepSize, blockSize));
    int track = -1;
    Vamp::Plugin::OutputList list = plugin.getOutputDescriptors();
    for (size_t i = 0; i < list.size(); ++i) {
        if (list[i].identifier == "smoothedpitchtrack") track = int(i);
    }
    BOOST_REQUIRE(track >= 0);
    vector<float> block(blockSize);
    int during = 0;
    double worst = 0;
    for (size_t start = 0; start < audio.size(); start += stepSize) {
        for (size_t i = 0; i < blockSize; ++i) {
            block[i] = start + i < audio.size() ? audio[start + i] : 0.f;
        }
        const float *buffers[1] = { &block[0] };
        Vamp::Plugin::FeatureSet fs = plugin.process
            (buffers, Vamp::RealTime::frame2RealTime(long(start), rate));
        for (const auto &f : fs[track]) {
            ++during;
            worst = std::max(worst, double(start) / rate - seconds(f.timestamp));
        }
    }
    Vamp::Plugin::FeatureSet rest = plugin.getRemainingFeatures();
    BOOST_CHECK_GT(during, int(5.0 * rate / stepSize));
    BOOST_CHECK_LE(worst, 443.0 * stepSize / rate);
    BOOST_CHECK(!rest[track].empty());
}

BOOST_AUTO_TEST_CASE(shortInputsAndFullViterbi)
{
    // Less input than the fixed lag, or than the gate's window, and the
    // full Viterbi of fixedlag 0: the voicing is the strongest's
    vector<float> shortThird = mix(tone(220, 0.2, 0.3), tone(275, 0.3, 0.3));
    vector<float> third = mix(tone(220, 0.2, 3.0), tone(275, 0.3, 3.0));
    for (bool full : { false, true }) {
        for (const vector<float> *audio : { &shortThird, &third }) {
            BOOST_TEST_CONTEXT((full ? "full Viterbi, " : "fixed lag, ")
                               << audio->size() << " samples") {
                Params strongest = following(0), lower = following(1);
                strongest["outputunvoiced"] = lower["outputunvoiced"] = 2;
                if (full) strongest["fixedlag"] = lower["fixedlag"] = 0;
                vector<pair<double, double> > a = pitchTrack(*audio, strongest);
                vector<pair<double, double> > b = pitchTrack(*audio, lower);
                BOOST_REQUIRE_EQUAL(a.size(), b.size());
                for (size_t i = 0; i < a.size(); ++i) {
                    BOOST_CHECK_EQUAL(a[i].first, b[i].first);
                    BOOST_CHECK_EQUAL(a[i].second > 0, b[i].second > 0);
                }
                if (audio == &third) {
                    Count low = count(b, 0.5, 2.5, { 220 });
                    BOOST_CHECK_GE(low.onVoice, 2 * framesInSpan * 9 / 10);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
