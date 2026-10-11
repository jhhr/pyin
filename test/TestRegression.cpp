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

// The regression test of test/regression/regression.sh, run in the
// test executable so that it needs no host. The plugin is fed
// testdata/bob_02.wav as sonic-annotator feeds it (blocks of the
// preferred size from the first frame, the last ones padded with
// zeros), with its default parameters, and its smoothed pitch track and
// notes, with fixed-lag smoothing off and on, are compared with the
// expected CSV files of test/regression. Those hold six significant
// digits, so values are compared to within a part in 10^5, which no
// change to the analysis gets under.

#include "PYinVamp.h"

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Where testdata/ and test/ are. The Makefiles run the tests from the
// source directory; a build elsewhere says where it is
#ifndef PYIN_SOURCE_DIR
#define PYIN_SOURCE_DIR "."
#endif

using std::string;
using std::vector;

typedef vector<vector<double> > Rows;

static uint32_t le32(const unsigned char *p)
{
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8) |
        (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}

static uint16_t le16(const unsigned char *p)
{
    return uint16_t(p[0] | (p[1] << 8));
}

// A RIFF WAVE file of 16-bit PCM or 32-bit float, its channels mixed
// to their mean
static bool readWav(string path, vector<float> &mono, float &rate)
{
    std::ifstream f(path.c_str(), std::ios::binary);
    if (!f) return false;
    vector<unsigned char> data((std::istreambuf_iterator<char>(f)),
                               std::istreambuf_iterator<char>());
    if (data.size() < 12 || memcmp(&data[0], "RIFF", 4) ||
        memcmp(&data[8], "WAVE", 4)) {
        return false;
    }
    int format = 0, channels = 0, bits = 0;
    size_t pos = 12;
    while (pos + 8 <= data.size()) {
        uint32_t size = le32(&data[pos + 4]);
        const unsigned char *body = &data[pos + 8];
        if (pos + 8 + size > data.size()) size = uint32_t(data.size() - pos - 8);
        if (!memcmp(&data[pos], "fmt ", 4) && size >= 16) {
            format = le16(body);
            channels = le16(body + 2);
            rate = float(le32(body + 4));
            bits = le16(body + 14);
        } else if (!memcmp(&data[pos], "data", 4)) {
            if (channels < 1) return false;
            if (format == 1 && bits == 16) {
                size_t n = size / (2 * channels);
                mono.assign(n, 0.f);
                for (size_t i = 0; i < n; ++i) {
                    float sum = 0.f;
                    for (int c = 0; c < channels; ++c) {
                        sum += int16_t(le16(body + 2 * (i * channels + c))) / 32768.f;
                    }
                    mono[i] = sum / float(channels);
                }
                return true;
            }
            if (format == 3 && bits == 32) {
                size_t n = size / (4 * channels);
                mono.assign(n, 0.f);
                for (size_t i = 0; i < n; ++i) {
                    float sum = 0.f;
                    for (int c = 0; c < channels; ++c) {
                        uint32_t u = le32(body + 4 * (i * channels + c));
                        float v;
                        memcpy(&v, &u, 4);
                        sum += v;
                    }
                    mono[i] = sum / float(channels);
                }
                return true;
            }
            return false;
        }
        pos += 8 + size + (size & 1);
    }
    return false;
}

static bool readCsv(string path, Rows &rows)
{
    std::ifstream f(path.c_str());
    if (!f) return false;
    string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        vector<double> row;
        std::istringstream ls(line);
        string cell;
        while (std::getline(ls, cell, ',')) {
            row.push_back(atof(cell.c_str()));
        }
        rows.push_back(row);
    }
    return true;
}

static double seconds(const Vamp::RealTime &t)
{
    return t.sec + t.nsec / 1e9;
}

// The plugin's smoothed pitch track (time, f0) and notes (time,
// duration, f0), as sonic-annotator writes them to CSV
static void runPlugin(const vector<float> &audio, float rate, bool fixedLag,
                      Rows &pitch, Rows &notes)
{
    PYinVamp plugin(rate);
    plugin.setParameter("fixedlag", fixedLag ? 1.f : 0.f);
    size_t blockSize = plugin.getPreferredBlockSize();
    size_t stepSize = plugin.getPreferredStepSize();
    BOOST_REQUIRE(plugin.initialise(1, stepSize, blockSize));

    std::map<string, int> outputs;
    Vamp::Plugin::OutputList list = plugin.getOutputDescriptors();
    for (size_t i = 0; i < list.size(); ++i) {
        outputs[list[i].identifier] = int(i);
    }
    BOOST_REQUIRE(outputs.count("smoothedpitchtrack") && outputs.count("notes"));
    int pitchOutput = outputs["smoothedpitchtrack"];
    int notesOutput = outputs["notes"];

    Vamp::Plugin::FeatureSet fs;
    vector<float> block(blockSize);
    for (size_t start = 0; start < audio.size(); start += stepSize) {
        for (size_t i = 0; i < blockSize; ++i) {
            block[i] = start + i < audio.size() ? audio[start + i] : 0.f;
        }
        const float *buffers[1] = { &block[0] };
        Vamp::Plugin::FeatureSet out = plugin.process
            (buffers, Vamp::RealTime::frame2RealTime(long(start), int(lrintf(rate))));
        for (auto &p : out) {
            fs[p.first].insert(fs[p.first].end(), p.second.begin(), p.second.end());
        }
    }
    Vamp::Plugin::FeatureSet rest = plugin.getRemainingFeatures();
    for (auto &p : rest) {
        fs[p.first].insert(fs[p.first].end(), p.second.begin(), p.second.end());
    }

    for (const auto &f : fs[pitchOutput]) {
        if (f.values.empty()) continue;
        pitch.push_back({ seconds(f.timestamp), double(f.values[0]) });
    }
    for (const auto &f : fs[notesOutput]) {
        if (f.values.empty()) continue;
        notes.push_back({ seconds(f.timestamp), seconds(f.duration), double(f.values[0]) });
    }
}

static void compare(string name, const Rows &obtained, const Rows &expected)
{
    BOOST_TEST_CONTEXT(name) {
        BOOST_REQUIRE_EQUAL(obtained.size(), expected.size());
        int reported = 0;
        for (size_t i = 0; i < expected.size(); ++i) {
            BOOST_REQUIRE_EQUAL(obtained[i].size(), expected[i].size());
            for (size_t j = 0; j < expected[i].size(); ++j) {
                double e = expected[i][j], o = obtained[i][j];
                if (std::fabs(o - e) > 1e-5 * std::fabs(e) + 1e-6) {
                    // The first few differences say enough
                    if (++reported <= 5) {
                        BOOST_ERROR("row " << i + 1 << " column " << j + 1
                                    << ": expected " << e << ", obtained " << o);
                    }
                }
            }
        }
        BOOST_CHECK_MESSAGE(reported == 0, reported << " values differ");
    }
}

static void regression(bool fixedLag)
{
    string root = PYIN_SOURCE_DIR;
    vector<float> audio;
    float rate = 0.f;
    if (!readWav(root + "/testdata/bob_02.wav", audio, rate)) {
        // As regression.sh: a source package without testdata/
        BOOST_TEST_MESSAGE("No " << root << "/testdata/bob_02.wav, skipping regression test");
        return;
    }
    Rows pitch, notes;
    runPlugin(audio, rate, fixedLag, pitch, notes);

    string suffix = fixedLag ? "-fixedlag" : "";
    string dir = root + "/test/regression/";
    Rows expectedPitch, expectedNotes;
    BOOST_REQUIRE(readCsv(dir + "expected-smoothedpitchtrack" + suffix + ".csv", expectedPitch));
    BOOST_REQUIRE(readCsv(dir + "expected-notes" + suffix + ".csv", expectedNotes));
    compare("smoothedpitchtrack" + suffix, pitch, expectedPitch);
    compare("notes" + suffix, notes, expectedNotes);
}

BOOST_AUTO_TEST_SUITE(TestRegression)

BOOST_AUTO_TEST_CASE(smoothedPitchTrackAndNotes)
{
    regression(false);
}

BOOST_AUTO_TEST_CASE(smoothedPitchTrackAndNotesFixedLag)
{
    regression(true);
}

BOOST_AUTO_TEST_SUITE_END()
