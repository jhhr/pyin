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

#ifndef _VOICEFOLLOWER_H_
#define _VOICEFOLLOWER_H_

#include "HarmonicSalience.h"
#include "MonoPitchHMM.h"

#include <vamp-sdk/RealTime.h>

#include <deque>
#include <utility>
#include <vector>

/**
 * Makes the pitch track follow one voice where two or three sing: the
 * lower, the higher or the middle one. The pitch HMM alone follows the
 * strongest, and so changes voice wherever the voices' levels cross,
 * note to note.
 *
 * A second pitch HMM decodes the same frames with their candidates
 * reshaped to prefer the chosen voice:
 *  - fill-in: a frame takes the candidates of the ten frames each side
 *    of it, fading by 0.85 a frame, as a voice the salience misses for a
 *    frame or two is what makes the track change voice most often;
 *  - presence: every candidate likely to be a voice of its own
 *    (HarmonicSalience::voiceLikelihoods()) is made as probable as the
 *    strongest, so that continuity counts for more than level;
 *  - preference: each candidate keeps its probability as far as no
 *    voice sings below it (Lower) or above it (Higher) over those
 *    frames, or as far as it is the middle one of three where one is
 *    likely (Middle).
 *
 * The track takes that path's pitch only inside passages where two or
 * more voices sing, as a gate finds them: a frame has two voices where
 * the first path is voiced and another voice, 1 to 9.5 semitones from
 * its pitch, is at least 0.85 likely; a passage is where that holds in
 * 80% of the first path's voiced frames over two seconds, at least a
 * fifth of them voiced, widened by a second each way. Elsewhere the
 * track is the first path's, and its voicing is the first path's
 * everywhere: in a solo, the lowest "voice" would be an echo or an
 * instrument's bleed.
 *
 * A track entry comes out once the gate and the second path are known
 * around it, at most 443 frames (2.6 s at a hop of 256) after its own
 * frame, so that a host still gets the track as the analysis runs.
 */
class VoiceFollower
{
public:
    enum Voice { Lower = 1, Higher = 2, Middle = 3 };

    /// A voice for the gate: its MIDI pitch, and its likelihood by
    /// gateRule()
    struct GateVoice {
        float pitch;
        float likelihood;
    };

    /// A track entry: its frame, frequency (not positive where
    /// unvoiced) and timestamp
    struct Entry {
        int frame;
        float freq;
        Vamp::RealTime timestamp;
    };

    /// fixedLag, voicedFloor and unvoicedHold as the first path's
    /// MonoPitchHMM has them
    VoiceFollower(Voice voice, int fixedLag,
                  double voicedFloor, double unvoicedHold);

    Voice voice() const { return m_voice; }

    /// How voiceLikelihoods() is to find the voices the preference
    /// chooses among. Higher has to tell a voice from the partials of
    /// the voices below it, which a voice above them looks like
    static HarmonicSalience::VoiceRule preferenceRule(Voice voice);

    /// How it is to find the voices for the gate
    static HarmonicSalience::VoiceRule gateRule();

    /**
     * The next frame, given before the first path gets it: its
     * candidates (MIDI pitch, probability) as the first path gets them,
     * each one's likelihood of being a voice by preferenceRule(), the
     * frame's voices for the gate, and the factor the first path's voiced
     * floor is scaled by (the level gate's, 1 with none), which scales the
     * floor here as well and counts against the frame's voices in the
     * preference.
     */
    void addFrame(const std::vector<std::pair<double, double> > &pitchProb,
                  const std::vector<double> &likelihoods,
                  const std::vector<GateVoice> &gateVoices,
                  double levelGate);

    /// The frames given so far
    int frameCount() const { return m_frames; }

    /// The first path's entry for a frame. In fixed-lag mode the first
    /// frame has none and frame n - lag two, as the first path gives them
    void addFirst(const Entry &entry);

    /// The end of the input: the frames still waiting for those after
    /// them, and the second path's last lag
    void finish();

    /// The entries that can come out now (after finish(), all the
    /// rest), in order
    std::vector<Entry> takeEntries();

private:
    // A frame waiting for the candidates of the frames after it
    struct RankFrame {
        std::vector<std::pair<double, double> > pitchProb;
        std::vector<double> q;
        std::vector<GateVoice> voices;
        double gate;
    };

    void feedRankFrame();
    // pitchProb, the frame's candidates as the second path weighs them
    // (filled in from the frames around); own, its candidates as they came
    void feedSecondPath(const std::vector<std::pair<double, double> > &pitchProb,
                        const std::vector<std::pair<double, double> > &own,
                        double gate);
    bool twoVoices(int frame) const;
    float secondPitch(int frame, int n) const;

    Voice m_voice;
    double m_voicedFloor;

    std::deque<RankFrame> m_rankFrames;
    size_t m_rankNext;

    MonoPitchHMM m_hmm;
    std::deque<std::vector<std::pair<double, double> > > m_pending;
    int m_fed;
    int m_next;

    int m_frames;
    bool m_finished;

    // Per frame: the first path's frequency and the second path's (0
    // where there is no estimate, as for the first frame in fixed-lag
    // mode), and how far each is known
    std::vector<float> m_first;
    std::vector<float> m_second;
    int m_firstKnown;
    int m_secondKnown;

    // Per frame whose candidates are filled in: the pitches of the
    // voices the gate counts (likely enough), frame k's from
    // m_voiceStart[k] to m_voiceStart[k + 1]
    std::vector<float> m_voices;
    std::vector<size_t> m_voiceStart;

    // The gate so far: prefix sums of the voiced frames and of those
    // with two voices, and of the frames whose window holds a passage
    std::vector<int> m_voicedSum;
    std::vector<int> m_twoSum;
    std::vector<int> m_windowSum;

    std::deque<Entry> m_entries;
};

#endif
