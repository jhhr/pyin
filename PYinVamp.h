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

#ifndef _PYINVAMP_H_
#define _PYINVAMP_H_

#include <vamp-sdk/Plugin.h>

#include "Yin.h"
#include "MonoPitchHMM.h"
#include "HarmonicSalience.h"
#include "LevelGate.h"
#include "VoiceFollower.h"

#include <memory>

class PYinVamp : public Vamp::Plugin
{
public:
    PYinVamp(float inputSampleRate);
    virtual ~PYinVamp();

    std::string getIdentifier() const;
    std::string getName() const;
    std::string getDescription() const;
    std::string getMaker() const;
    int getPluginVersion() const;
    std::string getCopyright() const;

    InputDomain getInputDomain() const;
    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(std::string identifier) const;
    void setParameter(std::string identifier, float value);

    ProgramList getPrograms() const;
    std::string getCurrentProgram() const;
    void selectProgram(std::string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    size_t m_channels;
    size_t m_stepSize;
    size_t m_blockSize;
    float m_fmin;
    float m_fmax;
    Yin m_yin;
    
    mutable int m_oF0Candidates;
    mutable int m_oF0Probs;
    mutable int m_oVoicedProb;
    mutable int m_oCandidateSalience;
    mutable int m_oSmoothedPitchTrack;
    mutable int m_oNotes;

    float m_threshDistr;
    float m_fixedLag;
    float m_outputUnvoiced;
    float m_preciseTime;
    float m_lowAmp;
    float m_onsetSensitivity;
    float m_pruneThresh;

    MonoPitchHMM m_pitchHmm;

    deque<vector<pair<double, double> > > m_pitchProb;
    deque<Vamp::RealTime> m_timestamp;
    vector<float> m_level;
    vector<float> m_pitchTrack;

    // The analysis of voices with harmonies, doubles and echoes. Each
    // parameter's default, 0, is pYIN 1.2 exactly.

    // Candidates from a harmonic salience of the spectrum (see
    // HarmonicSalience.h) in place of YIN's. YIN still runs: its
    // probability that the frame is voiced is what they share
    float m_salienceCandidates;

    // How YIN's dips get their probabilities: 0, every dip (pYIN 1.2);
    // 1, the dips in order of period (see YinUtil::yinProbOrdered())
    float m_candidateRule;

    // With salience candidates, how much of its salience a candidate
    // loses that looks like the common fundamental of two voices (see
    // HarmonicSalience::penaliseCommonFundamentals())
    float m_commonF0;

    // Frames more than this many dB below the level of the last ten
    // seconds' singing are faded out of voicing (see LevelGate.h)
    float m_levelFloor;

    // MonoPitchHMM's voiced floor (scaled by the level gate, so that it
    // bridges gaps in loud singing and not into a quiet tail) and its
    // unvoiced hold
    float m_voicedFloor;
    float m_voiceHold;

    // With salience candidates, which voice the track follows where two
    // or three sing: 0 the strongest, as the pitch HMM alone does; 1 the
    // lower, 2 the higher, 3 the middle one (see VoiceFollower.h)
    float m_voiceFollow;

    std::unique_ptr<HarmonicSalience> m_salience;
    std::unique_ptr<LevelGate> m_levelGate;
    std::unique_ptr<VoiceFollower> m_follower;

    // With a voice to follow, the frame salienceCandidates() last gave:
    // each candidate's likelihood of being a voice by the follower's
    // rule, and its voices for the follower's gate
    vector<double> m_lastLikelihoods;
    vector<VoiceFollower::GateVoice> m_lastGateVoices;

    // Frames whose salience window reaches past the blocks given so far
    struct WaitingFrame {
        Yin::YinOutput yo;
        float rms;
        Vamp::RealTime timestamp;
    };
    deque<WaitingFrame> m_waiting;

    // Feeds a frame's candidates to the pitch HMM and adds its outputs
    void addFrame(const Yin::YinOutput &yo, float rms,
                  Vamp::RealTime timestamp, FeatureSet &fs);

    // The salience candidates of the next frame, with their share of
    // the probability that it is voiced
    vector<pair<double, double> > salienceCandidates
        (const vector<pair<double, double> > &yinFreqProb, float rms);

    // Adds the follower's track entries to the pitch track and its output
    void addEntries(const vector<VoiceFollower::Entry> &entries,
                    FeatureSet &fs);

    void addNoteFeatures(FeatureSet &fs);
};

#endif
