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

#include "VoiceFollower.h"

#include <algorithm>
#include <cmath>

using std::vector;
using std::pair;

namespace {

// The fill-in: the frames each side whose candidates a frame takes,
// which are also those whose voices count for its rank among them, and
// how much each frame of distance fades them
const size_t fillFrames = 10;
const float fillDecay = 0.85f;

// Another voice below or above a candidate is 1 to 19 semitones from it
const double rankGap = 1.0;
const double rankRange = 19.0;

// The gate's voices of a frame: those at all likely, at least a
// semitone apart, up to eight; of them, those likely enough count
const double keepLikelihood = 0.1;
const size_t keepVoices = 8;
const float voiceGap = 1.f;
const float gateLikelihood = 0.85f;
const float gateLow = 1.f;
const float gateHigh = 9.5f;

// The passage window, centred on its frame; the share of its voiced
// frames that must have two voices, and of all of its frames that must
// be voiced; and how far a passage is widened each way
const int gateWindow = 344;
const float gateShare = 0.8f;
const float gateMinVoiced = 0.2f;
const int gateWiden = 172;

// Inside a passage, where the second path is unvoiced for a moment: its
// nearest voiced pitch within this many frames
const int secondHold = 20;

}

VoiceFollower::VoiceFollower(Voice voice, int fixedLag,
                             double voicedFloor, double unvoicedHold) :
    m_voice(voice),
    m_voicedFloor(voicedFloor),
    m_rankNext(0),
    m_hmm(fixedLag, voicedFloor, unvoicedHold),
    m_fed(0),
    m_next(1),
    m_frames(0),
    m_finished(false),
    m_firstKnown(0),
    m_secondKnown(0),
    m_voicedSum(1, 0),
    m_twoSum(1, 0),
    m_windowSum(1, 0)
{
}

HarmonicSalience::VoiceRule
VoiceFollower::preferenceRule(Voice voice)
{
    // The values as they were tuned, some in single precision
    HarmonicSalience::VoiceRule rule;
    rule.ownPartialResolution = 0.5;
    rule.explainLow = 0.3f;
    rule.explainHigh = 0.6f;
    rule.explainAccept = 0.5f;
    if (voice == Higher) {
        // A voice above another is also any of the lower one's partials
        // that it does not explain away
        rule.relLow = 0.01f;
        rule.relHigh = 0.1f;
        rule.harmonicPenalty = 1;
        rule.highHarmonicPenalty = 1;
        rule.harmonicOfAny = true;
        rule.explain = true;
        rule.explainPower = 2;
        rule.explainWhole = true;
    } else {
        rule.relLow = 0.002f;
        rule.relHigh = 0.02f;
        rule.harmonicPenalty = 0.8f;
        rule.highHarmonicPenalty = 0.8f;
    }
    return rule;
}

HarmonicSalience::VoiceRule
VoiceFollower::gateRule()
{
    HarmonicSalience::VoiceRule rule;
    rule.relLow = 0.002f;
    rule.relHigh = 0.02f;
    rule.ownPartialResolution = 0.5;
    return rule;
}

void
VoiceFollower::addFrame(const vector<pair<double, double> > &pitchProb,
                        const vector<double> &likelihoods,
                        const vector<GateVoice> &gateVoices,
                        double levelGate)
{
    ++m_frames;
    m_first.push_back(0.f);

    RankFrame rf;
    rf.pitchProb = pitchProb;
    rf.q = likelihoods;
    rf.q.resize(pitchProb.size(), 0.0);
    rf.voices = gateVoices;
    rf.gate = levelGate;
    m_rankFrames.push_back(rf);

    // A frame's fill-in and rank need the frames after it
    while (m_rankNext < m_rankFrames.size() &&
           m_rankFrames.size() - 1 - m_rankNext >= fillFrames) {
        feedRankFrame();
    }
}

void
VoiceFollower::addFirst(const Entry &entry)
{
    m_first[entry.frame] = entry.freq;
    m_firstKnown = std::max(m_firstKnown, entry.frame + 1);
    m_entries.push_back(entry);
}

void
VoiceFollower::feedRankFrame()
{
    const size_t c = m_rankNext;
    const RankFrame &rf = m_rankFrames[c];
    vector<pair<double, double> > pp = rf.pitchProb;
    vector<double> q = rf.q;
    vector<GateVoice> gv = rf.voices;

    double psum = 0;
    for (size_t i = 0; i < pp.size(); ++i) psum += pp[i].second;

    // Fill-in: the candidates of the frames around, fading with their
    // distance, join this frame's and share its voiced mass, so that a
    // voice missing from a frame or two is still there
    if (psum > 0) {
        size_t lo = c >= fillFrames ? c - fillFrames : 0;
        size_t hi = std::min(m_rankFrames.size() - 1, c + fillFrames);
        for (size_t s = lo; s <= hi; ++s) {
            if (s == c) continue;
            const RankFrame &o = m_rankFrames[s];
            double osum = 0;
            for (const auto &x : o.pitchProb) osum += x.second;
            if (osum <= 0) continue;
            double w = std::pow(double(fillDecay), double(s > c ? s - c : c - s));
            for (size_t j = 0; j < o.pitchProb.size(); ++j) {
                double pj = o.pitchProb[j].second / osum * psum * w;
                double qj = o.q[j] * w;
                bool merged = false;
                for (size_t k = 0; k < pp.size(); ++k) {
                    if (std::fabs(pp[k].first - o.pitchProb[j].first) < 0.35) {
                        pp[k].second = std::max(pp[k].second, pj);
                        q[k] = std::max(q[k], qj);
                        merged = true;
                        break;
                    }
                }
                if (!merged) {
                    pp.push_back(pair<double, double>(o.pitchProb[j].first, pj));
                    q.push_back(qj);
                }
            }
            // and the gate's voices the same way
            for (const GateVoice &ov : o.voices) {
                float qj = float(ov.likelihood * w);
                bool merged = false;
                for (GateVoice &v : gv) {
                    if (std::fabs(v.pitch - ov.pitch) < 0.35f) {
                        v.likelihood = std::max(v.likelihood, qj);
                        merged = true;
                        break;
                    }
                }
                if (!merged) gv.push_back(GateVoice { ov.pitch, qj });
            }
        }
        double nsum = 0;
        for (const auto &x : pp) nsum += x.second;
        for (auto &x : pp) x.second *= psum / nsum;
    }

    // The gate's voices of this frame: the distinct ones, most likely
    // first, of which it keeps those likely enough to count. A frame
    // with no voiced mass has none
    {
        vector<int> order(gv.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = int(i);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            return gv[a].likelihood > gv[b].likelihood;
        });
        vector<int> picked;
        for (int i : order) {
            if (gv[i].likelihood < keepLikelihood) continue;
            bool distinct = true;
            for (int j : picked) {
                if (std::fabs(gv[i].pitch - gv[j].pitch) < voiceGap) distinct = false;
            }
            if (distinct) picked.push_back(i);
            if (picked.size() == keepVoices) break;
        }
        m_voiceStart.push_back(m_voices.size());
        if (psum > 0) {
            for (int i : picked) {
                if (gv[i].likelihood >= gateLikelihood) {
                    m_voices.push_back(gv[i].pitch);
                }
            }
        }
    }

    const int n = int(pp.size());

    // Presence: a candidate that is clearly a voice is as probable as
    // the strongest, so that of two voices at different levels neither
    // is much more likely
    if (n > 0 && psum > 0) {
        double pmax = 0;
        for (int i = 0; i < n; ++i) pmax = std::max(pmax, pp[i].second);
        double nsum = 0;
        for (int i = 0; i < n; ++i) {
            pp[i].second = std::max(pp[i].second, q[i] * pmax);
            nsum += pp[i].second;
        }
        for (int i = 0; i < n; ++i) pp[i].second *= psum / nsum;
    }

    if (n > 0 && psum > 0) {
        // For each candidate, the probability that no voice sings below
        // it and that none sings above it, by the voices of the frames
        // around as they came (not filled in), averaged over those
        // frames: a voice missing for a frame or two does not make the
        // next one up the lowest
        size_t lo = c >= fillFrames ? c - fillFrames : 0;
        size_t hi = std::min(m_rankFrames.size() - 1, c + fillFrames);
        vector<double> noBelow(n), noAbove(n);
        for (int i = 0; i < n; ++i) {
            double pitch = pp[i].first;
            double sumBelow = 0, sumAbove = 0;
            int frames = 0;
            for (size_t s = lo; s <= hi; ++s) {
                const RankFrame &o = m_rankFrames[s];
                double nb = 1, na = 1;
                for (size_t j = 0; j < o.pitchProb.size(); ++j) {
                    if (s == c && int(j) == i) continue;
                    double d = o.pitchProb[j].first - pitch;
                    double qj = o.q[j] * o.gate;
                    if (d <= -rankGap && d >= -rankRange) nb *= 1 - qj;
                    else if (d >= rankGap && d <= rankRange) na *= 1 - qj;
                }
                sumBelow += nb;
                sumAbove += na;
                ++frames;
            }
            noBelow[i] = sumBelow / frames;
            noAbove[i] = sumAbove / frames;
        }

        // How likely some candidate is a middle voice, with voices both
        // below and above it: where none is, the middle preference has
        // nothing to choose and leaves the frame alone
        double middlePresent = 0;
        if (m_voice == Middle) {
            for (int i = 0; i < n; ++i) {
                middlePresent = std::max(middlePresent,
                                         q[i] * (1 - noBelow[i]) * (1 - noAbove[i]));
            }
        }

        double nsum = 0;
        for (int i = 0; i < n; ++i) {
            double factor = 1;
            if (m_voice == Lower) {
                factor *= q[i] * noBelow[i];
            } else if (m_voice == Higher) {
                factor *= q[i] * noAbove[i];
            } else {
                double g = q[i] * (1 - noBelow[i]) * (1 - noAbove[i]);
                factor *= 1 - (middlePresent - g);
            }
            pp[i].second *= factor;
            nsum += pp[i].second;
        }
        // The frame keeps its voiced mass: what the voices not preferred
        // lose goes to the one preferred
        if (nsum > 0) {
            double target = nsum + (psum - nsum);
            for (int i = 0; i < n; ++i) pp[i].second *= target / nsum;
        }
    }

    feedSecondPath(pp, rf.pitchProb, rf.gate);

    ++m_rankNext;
    while (m_rankNext > fillFrames) {
        m_rankFrames.pop_front();
        --m_rankNext;
    }
}

void
VoiceFollower::feedSecondPath(const vector<pair<double, double> > &pitchProb,
                              const vector<pair<double, double> > &own,
                              double gate)
{
    // the floor bridges gaps in loud singing, not into a tail, as on the
    // first path
    m_hmm.m_voicedFloor = m_voicedFloor * gate;
    vector<double> obs = m_hmm.calculateObsProb(pitchProb);
    if (m_fed == 0) {
        m_hmm.initialise(obs);
    } else {
        m_hmm.process(obs);
    }
    ++m_fed;
    m_second.push_back(0.f);
    // The path's pitch is looked up among the frame's own candidates:
    // those filled in from the frames around are the voice where it was
    // then, tens of cents away in a vibrato. With none of its own near
    // it, the frame keeps the path's own pitch (MonoPitchHMM::nearestFreq())
    m_pending.push_back(own);
    int lag = m_hmm.m_fixedLag;
    if (lag > 0 && int(m_pending.size()) == lag + 1) {
        // as on the first path, the first frame has no estimate
        m_pending.pop_front();
        vector<int> path = m_hmm.track();
        m_second[m_next++] = m_hmm.nearestFreq(path[0], m_pending[0]);
        m_secondKnown = m_next;
    }
}

void
VoiceFollower::finish()
{
    while (m_rankNext < m_rankFrames.size()) {
        feedRankFrame();
    }

    // the second path's last lag (or, without one, the whole of it)
    if (m_fed > 0) {
        vector<int> path = m_hmm.track();
        int n = m_frames;
        int from = n - int(path.size());
        for (size_t k = 0; k < path.size(); ++k) {
            int frame = from + int(k);
            if (frame < 0 || frame >= int(m_second.size())) continue;
            const auto &candidates = m_pending[m_pending.size() - path.size() + k];
            m_second[frame] = m_hmm.nearestFreq(path[k], candidates);
        }
    }

    m_finished = true;
}

bool
VoiceFollower::twoVoices(int k) const
{
    if (m_first[k] <= 0 || k >= int(m_voiceStart.size())) return false;
    double p1 = 69 + 12 * std::log2(m_first[k] / 440.0);
    size_t a = m_voiceStart[k];
    size_t b = k + 1 < int(m_voiceStart.size()) ? m_voiceStart[k + 1] : m_voices.size();
    for (size_t j = a; j < b; ++j) {
        double d = m_voices[j] - p1;
        if (std::fabs(d) < gateLow || std::fabs(d) > gateHigh) continue;
        return true;
    }
    return false;
}

float
VoiceFollower::secondPitch(int k, int n) const
{
    float f2 = k < int(m_second.size()) ? m_second[k] : 0.f;
    if (f2 > 0) return f2;
    if (f2 == 0) return m_first[k]; // no estimate (the first frame)

    // Where the second path is unvoiced for a moment: its nearest voiced
    // pitch within secondHold frames inside the first path's voiced run,
    // earlier first
    bool left = true, right = true;
    for (int d = 1; d <= secondHold && (left || right); ++d) {
        int j = k - d;
        if (left) {
            if (j < 0 || m_first[j] <= 0) left = false;
            else if (m_second[j] > 0) return m_second[j];
        }
        j = k + d;
        if (right) {
            if (j >= n || m_first[j] <= 0) right = false;
            else if (m_second[j] > 0) return m_second[j];
        }
    }
    // else the pitch its unvoiced state holds, unless that has wandered
    // more than an octave from the first path's
    if (std::fabs(12 * std::log2(-f2 / m_first[k])) > 12) return m_first[k];
    return -f2;
}

vector<VoiceFollower::Entry>
VoiceFollower::takeEntries()
{
    const int n = m_frames;
    const bool end = m_finished;
    const int W = gateWindow;
    const int R = gateWiden;

    // Whether each frame has two voices, as far as its first path and
    // its voices are known
    int known = end ? n : std::min(m_firstKnown, int(m_voiceStart.size()));
    while (int(m_twoSum.size()) - 1 < known) {
        int k = int(m_twoSum.size()) - 1;
        m_voicedSum.push_back(m_voicedSum.back() + (m_first[k] > 0));
        m_twoSum.push_back(m_twoSum.back() + twoVoices(k));
    }
    int twoKnown = int(m_twoSum.size()) - 1;

    // Whether each frame's window holds a passage, as far as the frames
    // it reaches are known
    while (int(m_windowSum.size()) - 1 < n) {
        int i = int(m_windowSum.size()) - 1;
        if (!end && i + W - W / 2 > twoKnown) break;
        int a = std::max(0, i - W / 2), b = std::min(n, i + W - W / 2);
        int v1 = m_voicedSum[b] - m_voicedSum[a], v2 = m_twoSum[b] - m_twoSum[a];
        bool two = v1 > 0 && v1 >= gateMinVoiced * (b - a) && v2 >= gateShare * v1;
        m_windowSum.push_back(m_windowSum.back() + two);
    }
    int windowKnown = int(m_windowSum.size()) - 1;

    vector<Entry> out;
    while (!m_entries.empty()) {
        Entry e = m_entries.front();
        int k = e.frame;
        // inside a passage if a window within the widening holds one
        if (!end && k + R + 1 > windowKnown) break;
        int a = std::max(0, k - R), b = std::min(n, k + R + 1);
        bool gate = m_windowSum[b] - m_windowSum[a] > 0;
        if (e.freq > 0 && gate) {
            // there, where the first path is voiced, the second path's
            // pitch, once known around it
            if (!end && (k + secondHold >= m_secondKnown ||
                         k + secondHold >= m_firstKnown)) break;
            e.freq = secondPitch(k, n);
        }
        out.push_back(e);
        m_entries.pop_front();
    }
    return out;
}
