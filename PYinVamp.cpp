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

#include "PYinVamp.h"
#include "MonoNote.h"
#include "MonoPitchHMM.h"

#include <vector>
#include <algorithm>

#include <cstdio>
#include <cmath>
#include <complex>

using std::string;
using std::vector;
using Vamp::RealTime;


PYinVamp::PYinVamp(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_stepSize(256),
    m_blockSize(2048),
    m_fmin(40),
    m_fmax(1600),
    m_yin(2048, inputSampleRate, 0.0),
    m_oF0Candidates(0),
    m_oF0Probs(0),
    m_oVoicedProb(0),
    m_oCandidateSalience(0),
    m_oSmoothedPitchTrack(0),
    m_oNotes(0),
    m_threshDistr(2.0f),
    m_fixedLag(1.0f),
    m_outputUnvoiced(0.0f),
    m_preciseTime(0.0f),
    m_lowAmp(0.1f),
    m_onsetSensitivity(0.7f),
    m_pruneThresh(0.1f),
    m_pitchHmm(0),
    m_pitchProb(0),
    m_timestamp(0),
    m_level(0),
    m_pitchTrack(0),
    m_salienceCandidates(0.0f),
    m_candidateRule(0.0f),
    m_commonF0(0.0f),
    m_levelFloor(0.0f),
    m_voicedFloor(0.0f),
    m_voiceHold(0.0f),
    m_voiceFollow(0.0f)
{
}

PYinVamp::~PYinVamp()
{
}

string
PYinVamp::getIdentifier() const
{
    return "pyin";
}

string
PYinVamp::getName() const
{
    return "pYin";
}

string
PYinVamp::getDescription() const
{
    return "Monophonic pitch and note tracking based on a probabilistic Yin extension.";
}

string
PYinVamp::getMaker() const
{
    return "Matthias Mauch";
}

int
PYinVamp::getPluginVersion() const
{
    // Increment this each time you release a version that behaves
    // differently from the previous one
    return 3;
}

string
PYinVamp::getCopyright() const
{
    return "GPL";
}

PYinVamp::InputDomain
PYinVamp::getInputDomain() const
{
    return TimeDomain;
}

size_t
PYinVamp::getPreferredBlockSize() const
{
    return 2048;
}

size_t 
PYinVamp::getPreferredStepSize() const
{
    return 256;
}

size_t
PYinVamp::getMinChannelCount() const
{
    return 1;
}

size_t
PYinVamp::getMaxChannelCount() const
{
    return 1;
}

PYinVamp::ParameterList
PYinVamp::getParameterDescriptors() const
{
    ParameterList list;
    
    ParameterDescriptor d;

    d.identifier = "threshdistr";
    d.name = "Yin threshold distribution";
    d.description = ".";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 7.0f;
    d.defaultValue = 2.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    d.valueNames.push_back("Uniform");
    d.valueNames.push_back("Beta (mean 0.10)");
    d.valueNames.push_back("Beta (mean 0.15)");
    d.valueNames.push_back("Beta (mean 0.20)");
    d.valueNames.push_back("Beta (mean 0.30)");
    d.valueNames.push_back("Single Value 0.10");
    d.valueNames.push_back("Single Value 0.15");
    d.valueNames.push_back("Single Value 0.20");
    list.push_back(d);

    d.valueNames.clear();

    d.identifier = "fixedlag";
    d.name = "Fixed-lag smoothing";
    d.description = "Use fixed lag smoothing, not full Viterbi smoothing.";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 1.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    list.push_back(d);

    d.identifier = "outputunvoiced";
    d.valueNames.clear();
    d.name = "Output estimates classified as unvoiced?";
    d.description = ".";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 2.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    d.valueNames.push_back("No");
    d.valueNames.push_back("Yes");
    d.valueNames.push_back("Yes, as negative frequencies");
    list.push_back(d);

    d.identifier = "precisetime";
    d.valueNames.clear();
    d.name = "Use non-standard precise YIN timing (slow).";
    d.description = ".";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    list.push_back(d);

    d.identifier = "lowampsuppression";
    d.valueNames.clear();
    d.name = "Suppress low amplitude pitch estimates.";
    d.description = ".";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.1f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "onsetsensitivity";
    d.valueNames.clear();
    d.name = "Onset sensitivity";
    d.description = "Adds additional note onsets when RMS increases.";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.7f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "prunethresh";
    d.valueNames.clear();
    d.name = "Duration pruning threshold.";
    d.description = "Prune notes that are shorter than this value.";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 0.2f;
    d.defaultValue = 0.1f;
    d.isQuantized = false;
    list.push_back(d);

    // For voices with harmonies, doubles and echoes, as in vocals
    // separated from a mix. Each default is pYIN 1.2's behaviour.

    d.identifier = "candidaterule";
    d.valueNames.clear();
    d.name = "Candidate probabilities";
    d.description = "How the Yin candidates share the probability that a frame is voiced: by every dip of the difference function, or by the dips in order of period, mostly to the first, as Yin chooses, which avoids most octave-down errors.";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    d.valueNames.push_back("Every dip");
    d.valueNames.push_back("Dips in order");
    list.push_back(d);

    d.identifier = "salmode";
    d.valueNames.clear();
    d.name = "Pitch candidates";
    d.description = "Where the pitch track's candidates come from: the Yin difference function, or a harmonic salience of the spectrum, which finds each of two voices singing together where Yin finds the period they have in common, an octave or two below both.";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    d.valueNames.push_back("Yin");
    d.valueNames.push_back("Harmonic salience");
    list.push_back(d);

    d.identifier = "commonf0";
    d.valueNames.clear();
    d.name = "Common fundamental penalty";
    d.description = "With harmonic salience candidates: the share of its salience a candidate loses where two voices, a fifth or a fourth apart, make it look like their common fundamental and it has no partial of its own (0: none).";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "levelfloor";
    d.valueNames.clear();
    d.name = "Level floor";
    d.description = "Frames more than this many dB below the level of the last ten seconds' singing are not voiced, which removes most echo and reverb tails (0: none).";
    d.unit = "dB";
    d.minValue = 0.0f;
    d.maxValue = 60.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "voicedfloor";
    d.valueNames.clear();
    d.name = "Voiced gap bridging";
    d.description = "Lets the pitch track cross frames that have no candidate near it rather than break: the least likelihood of a voiced pitch, relative to that of no pitch (0: none).";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 1.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "voicehold";
    d.valueNames.clear();
    d.name = "Voice hold across gaps";
    d.description = "Makes the pitch track come back after an unvoiced gap near the pitch it left, rather than on another voice: how much more a move costs during the gap than staying (0: no more).";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 10.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "voicefollow";
    d.valueNames.clear();
    d.name = "Voice to follow";
    d.description = "Where two or three voices sing at once, which one the pitch track follows: the strongest, which makes the track change voice wherever their levels cross; the lower; the higher; or the middle one of three. Inside such passages a second pitch track decoded with that preference gives the pitch; the voicing, and the track outside them, are as with the strongest. Needs the salience candidates.";
    d.unit = "";
    d.minValue = 0.0f;
    d.maxValue = 3.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    d.valueNames.push_back("Strongest");
    d.valueNames.push_back("Lower");
    d.valueNames.push_back("Higher");
    d.valueNames.push_back("Middle");
    list.push_back(d);

    return list;
}

float
PYinVamp::getParameter(string identifier) const
{
    if (identifier == "threshdistr") {
            return m_threshDistr;
    }
    if (identifier == "fixedlag") {
            return m_fixedLag;
    }
    if (identifier == "outputunvoiced") {
            return m_outputUnvoiced;
    }
    if (identifier == "precisetime") {
            return m_preciseTime;
    }
    if (identifier == "lowampsuppression") {
            return m_lowAmp;
    }
    if (identifier == "onsetsensitivity") {
            return m_onsetSensitivity;
    }
    if (identifier == "prunethresh") {
            return m_pruneThresh;
    }
    if (identifier == "candidaterule") {
            return m_candidateRule;
    }
    if (identifier == "salmode") {
            return m_salienceCandidates;
    }
    if (identifier == "commonf0") {
            return m_commonF0;
    }
    if (identifier == "levelfloor") {
            return m_levelFloor;
    }
    if (identifier == "voicedfloor") {
            return m_voicedFloor;
    }
    if (identifier == "voicehold") {
            return m_voiceHold;
    }
    if (identifier == "voicefollow") {
            return m_voiceFollow;
    }
    return 0.f;
}

void
PYinVamp::setParameter(string identifier, float value) 
{
    if (identifier == "threshdistr")
    {
        m_threshDistr = value;
    }
    if (identifier == "fixedlag")
    {
        m_fixedLag = value;
    }
    if (identifier == "outputunvoiced")
    {
        m_outputUnvoiced = value;
    }
    if (identifier == "precisetime")
    {
        m_preciseTime = value;
    }
    if (identifier == "lowampsuppression")
    {
        m_lowAmp = value;
    }
    if (identifier == "onsetsensitivity")
    {
        m_onsetSensitivity = value;
    }
    if (identifier == "prunethresh")
    {
        m_pruneThresh = value;
    }
    if (identifier == "candidaterule")
    {
        m_candidateRule = value;
    }
    if (identifier == "salmode")
    {
        m_salienceCandidates = value;
    }
    if (identifier == "commonf0")
    {
        m_commonF0 = value;
    }
    if (identifier == "levelfloor")
    {
        m_levelFloor = value;
    }
    if (identifier == "voicedfloor")
    {
        m_voicedFloor = value;
    }
    if (identifier == "voicehold")
    {
        m_voiceHold = value;
    }
    if (identifier == "voicefollow")
    {
        m_voiceFollow = value;
    }
}

PYinVamp::ProgramList
PYinVamp::getPrograms() const
{
    ProgramList list;
    return list;
}

string
PYinVamp::getCurrentProgram() const
{
    return ""; // no programs
}

void
PYinVamp::selectProgram(string)
{
}

PYinVamp::OutputList
PYinVamp::getOutputDescriptors() const
{
    OutputList outputs;

    OutputDescriptor d;
    
    int outputNumber = 0;

    d.identifier = "f0candidates";
    d.name = "F0 Candidates";
    d.description = "Estimated fundamental frequency candidates.";
    d.unit = "Hz";
    d.hasFixedBinCount = false;
    d.hasKnownExtents = true;
    d.minValue = m_fmin;
    d.maxValue = 500;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oF0Candidates = outputNumber++;

    d.identifier = "f0probs";
    d.name = "Candidate Probabilities";
    d.description = "Probabilities of estimated fundamental frequency candidates.";
    d.unit = "";
    d.hasFixedBinCount = false;
    d.hasKnownExtents = true;
    d.minValue = 0;
    d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oF0Probs = outputNumber++;
    
    d.identifier = "voicedprob";
    d.name = "Voiced Probability";
    d.description = "Probability that the signal is voiced according to Probabilistic Yin.";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = true;
    d.minValue = 0;
    d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oVoicedProb = outputNumber++;

    d.identifier = "candidatesalience";
    d.name = "Candidate Salience";
    d.description = "Candidate Salience";
    d.hasFixedBinCount = true;
    d.binCount = m_blockSize / 2;
    d.hasKnownExtents = true;
    d.minValue = 0;
    d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oCandidateSalience = outputNumber++;
    
    d.identifier = "smoothedpitchtrack";
    d.name = "Smoothed Pitch Track";
    d.description = "Frame-by-frame pitch estimate after smoothing";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = false;
    outputs.push_back(d);
    m_oSmoothedPitchTrack = outputNumber++;

    d.identifier = "notes";
    d.name = "Notes";
    d.description = "Derived fixed-pitch note frequencies";
    d.unit = "Hz";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.sampleRate = (m_inputSampleRate / m_stepSize);
    d.hasDuration = true;
    outputs.push_back(d);
    m_oNotes = outputNumber++;

    return outputs;
}

bool
PYinVamp::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;
    
    reset();

    return true;
}

void
PYinVamp::reset()
{    
    m_yin.setThresholdDistr(m_threshDistr);
    m_yin.setFrameSize(m_blockSize);
    m_yin.setFast(!m_preciseTime);
    m_yin.setCandidateRule(m_candidateRule > 0.5f ? 1 : 0);

    int lag = m_fixedLag > 0.5f ? 100 : 0;
    m_pitchHmm = MonoPitchHMM(lag, m_voicedFloor, m_voiceHold);

    if (m_salienceCandidates > 0.5f) {
        // centred where the frame's timestamp is
        int offset = m_preciseTime == 1.0 ? m_blockSize/2 : m_blockSize/4;
        m_salience.reset(new HarmonicSalience
                         (m_inputSampleRate, int(m_blockSize), int(m_stepSize),
                          offset, HarmonicSalience::Params()));
    } else {
        m_salience.reset();
    }
    if (m_levelFloor > 0) {
        m_levelGate.reset(new LevelGate
                          (m_levelFloor, int(10.0f * m_inputSampleRate / m_stepSize)));
    } else {
        m_levelGate.reset();
    }
    int follow = int(lrintf(m_voiceFollow));
    if (m_salience && follow >= VoiceFollower::Lower && follow <= VoiceFollower::Middle) {
        m_follower.reset(new VoiceFollower(VoiceFollower::Voice(follow), lag,
                                           m_voicedFloor, m_voiceHold));
    } else {
        m_follower.reset();
    }
    m_lastLikelihoods.clear();
    m_lastGateVoices.clear();
    m_waiting.clear();
    
    m_pitchProb.clear();
    m_timestamp.clear();
    m_level.clear();
    m_pitchTrack.clear();
}

PYinVamp::FeatureSet
PYinVamp::process(const float *const *inputBuffers, RealTime timestamp)
{
    int offset = m_preciseTime == 1.0 ? m_blockSize/2 : m_blockSize/4;
    timestamp = timestamp + Vamp::RealTime::frame2RealTime(offset, 
        lrintf(m_inputSampleRate));

    FeatureSet fs;
    
    float rms = 0;
    
    double *dInputBuffers = new double[m_blockSize];
    for (size_t i = 0; i < m_blockSize; ++i) {
        dInputBuffers[i] = inputBuffers[0][i];
        rms += inputBuffers[0][i] * inputBuffers[0][i];
    }
    rms /= m_blockSize;
    rms = sqrt(rms);
    
    Yin::YinOutput yo = m_yin.processProbabilisticYin(dInputBuffers);
    delete [] dInputBuffers;

    m_level.push_back(yo.rms);

    if (!m_salience) {
        addFrame(yo, rms, timestamp, fs);
        return fs;
    }

    // A frame's salience window can reach past its block, so the
    // frame waits until the blocks that complete it have come
    m_salience->push(inputBuffers[0]);
    WaitingFrame waiting;
    waiting.yo = yo;
    waiting.rms = rms;
    waiting.timestamp = timestamp;
    m_waiting.push_back(waiting);
    while (!m_waiting.empty() && m_salience->ready()) {
        WaitingFrame &w = m_waiting.front();
        w.yo.freqProb = salienceCandidates(w.yo.freqProb, w.rms);
        addFrame(w.yo, w.rms, w.timestamp, fs);
        m_waiting.pop_front();
    }
    return fs;
}

void
PYinVamp::addFrame(const Yin::YinOutput &yo, float rms, RealTime timestamp,
                   FeatureSet &fs)
{
    bool isLowAmplitude = (rms < m_lowAmp);

    vector<pair<double, double> > tempPitchProb;
    for (size_t iCandidate = 0; iCandidate < yo.freqProb.size(); ++iCandidate)
    {
        double tempPitch = 12 * std::log(yo.freqProb[iCandidate].first/440)/std::log(2.) + 69;
        if (!isLowAmplitude)
        {
            tempPitchProb.push_back(pair<double, double>
                (tempPitch, yo.freqProb[iCandidate].second));
        } else {
            float factor = ((rms+0.01*m_lowAmp)/(1.01*m_lowAmp));
            tempPitchProb.push_back(pair<double, double>
                (tempPitch, yo.freqProb[iCandidate].second*factor));
        }
    }

    double gate = 1;
    if (m_levelGate) {
        gate = m_levelGate->process(rms);
        if (gate < 1) {
            for (size_t i = 0; i < tempPitchProb.size(); ++i) {
                tempPitchProb[i].second *= gate;
            }
        }
        // the floor bridges gaps in loud singing, not into a tail
        m_pitchHmm.m_voicedFloor = m_voicedFloor * gate;
    }

    // The follower's second path takes the same candidates
    if (m_follower) {
        m_follower->addFrame(tempPitchProb, m_lastLikelihoods,
                             m_lastGateVoices, gate);
    }

    vector<double> tempObsProb = m_pitchHmm.calculateObsProb(tempPitchProb);
    if (m_timestamp.empty())
    {
        m_pitchHmm.initialise(tempObsProb);
    } else {
        m_pitchHmm.process(tempObsProb);
    }

    m_pitchProb.push_back(tempPitchProb);
    m_timestamp.push_back(timestamp);

    int lag = m_pitchHmm.m_fixedLag;

    if (m_fixedLag > 0.5f) // do fixed-lag smoothing instead of full Viterbi
    {
        if (int(m_timestamp.size()) == lag + 1)
        {
            m_timestamp.pop_front();
            m_pitchProb.pop_front();

            Feature f;
            f.hasTimestamp = true;
            vector<int> rawPitchPath = m_pitchHmm.track();
            float freq = m_pitchHmm.nearestFreq(rawPitchPath[0], 
                                                m_pitchProb[0]);
            if (m_follower) {
                // the entry waits for the follower's gate
                VoiceFollower::Entry e;
                e.frame = m_follower->frameCount() - lag;
                e.freq = freq;
                e.timestamp = m_timestamp[0];
                m_follower->addFirst(e);
            } else {
                m_pitchTrack.push_back(freq);
                f.timestamp = m_timestamp[0];
                f.values.clear();

                // different output modes
                if (freq < 0 && (m_outputUnvoiced==0))
                {

                } else {
                    if (m_outputUnvoiced == 1)
                    {
                        f.values.push_back(fabs(freq));
                    } else {
                        f.values.push_back(freq);
                    }
                    fs[m_oSmoothedPitchTrack].push_back(f);
                }
            }
        }
    }

    if (m_follower) {
        addEntries(m_follower->takeEntries(), fs);
    }

    // F0 CANDIDATES
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = timestamp;
    for (size_t i = 0; i < yo.freqProb.size(); ++i)
    {
        f.values.push_back(yo.freqProb[i].first);
    }
    fs[m_oF0Candidates].push_back(f);
    
    // VOICEDPROB
    f.values.clear();
    float voicedProb = 0;
    for (size_t i = 0; i < yo.freqProb.size(); ++i)
    {
        f.values.push_back(yo.freqProb[i].second);
        voicedProb += yo.freqProb[i].second;
    }
    fs[m_oF0Probs].push_back(f);
    
    f.values.clear();
    f.values.push_back(voicedProb);
    fs[m_oVoicedProb].push_back(f);

    // SALIENCE -- maybe this should eventually disappear
    f.values.clear();
    float salienceSum = 0;
    for (size_t iBin = 0; iBin < yo.salience.size(); ++iBin)
    {
        f.values.push_back(yo.salience[iBin]);
        salienceSum += yo.salience[iBin];
    }
    fs[m_oCandidateSalience].push_back(f);
}

PYinVamp::FeatureSet
PYinVamp::getRemainingFeatures()
{
    FeatureSet fs;

    // The frames still waiting for their salience, as if the audio
    // went on in silence
    while (!m_waiting.empty()) {
        WaitingFrame &w = m_waiting.front();
        w.yo.freqProb = salienceCandidates(w.yo.freqProb, w.rms);
        addFrame(w.yo, w.rms, w.timestamp, fs);
        m_waiting.pop_front();
    }

    if (m_pitchProb.empty()) {
        return fs;
    }

    Feature f;
    f.hasTimestamp = true;
    f.hasDuration = false;

    // ================== P I T C H  T R A C K =================================

    // NB we do this even in fixed-lag mode, as we still have the last
    // lag's-worth of pitch probs to consume
    
    vector<int> rawPitchPath = m_pitchHmm.track();
    
    for (size_t iFrame = 0; iFrame < rawPitchPath.size(); ++iFrame)
    {
        float freq = m_pitchHmm.nearestFreq(rawPitchPath[iFrame], 
                                            m_pitchProb[iFrame]);
        if (m_follower) {
            VoiceFollower::Entry e;
            e.frame = m_follower->frameCount() - int(rawPitchPath.size()) + int(iFrame);
            e.freq = freq;
            e.timestamp = m_timestamp[iFrame];
            m_follower->addFirst(e);
            continue;
        }
        m_pitchTrack.push_back(freq); // for note processing below
        
        f.timestamp = m_timestamp[iFrame];
        f.values.clear();
        
        // different output modes
        if (freq < 0 && (m_outputUnvoiced==0)) continue;
        if (m_outputUnvoiced == 1)
        {
            f.values.push_back(fabs(freq));
        } else {
            f.values.push_back(freq);
        }
        fs[m_oSmoothedPitchTrack].push_back(f);
    }

    if (m_follower) {
        m_follower->finish();
        addEntries(m_follower->takeEntries(), fs);
    }

    addNoteFeatures(fs);

    return fs;
}

void
PYinVamp::addEntries(const vector<VoiceFollower::Entry> &entries,
                     FeatureSet &fs)
{
    Feature f;
    f.hasTimestamp = true;
    f.hasDuration = false;
    for (const VoiceFollower::Entry &e : entries) {
        m_pitchTrack.push_back(e.freq); // for the notes at the end
        if (e.freq < 0 && m_outputUnvoiced == 0) continue;
        f.timestamp = e.timestamp;
        f.values.clear();
        f.values.push_back(m_outputUnvoiced == 1 ? fabs(e.freq) : e.freq);
        fs[m_oSmoothedPitchTrack].push_back(f);
    }
}

void
PYinVamp::addNoteFeatures(FeatureSet &fs)
{
    std::vector<std::vector<std::pair<double, double> > > smoothedPitch;
    for (size_t iFrame = 0; iFrame < m_pitchTrack.size(); ++iFrame) {
        std::vector<std::pair<double, double> > temp;
        if (m_pitchTrack[iFrame] > 0)
        {
            double tempPitch = 12 * 
                std::log(m_pitchTrack[iFrame]/440)/std::log(2.) + 69;
            temp.push_back(std::pair<double,double>(tempPitch, .9));
        }
        smoothedPitch.push_back(temp);
    }

    // In fixed-lag mode, we use fixed-lag processing for the note
    // transitions here as well as for the pitch transitions in
    // process. The main reason we provide the fixed-lag option is so
    // that we can get pitch results incrementally from process; we
    // don't get that outcome here, but we do benefit from its bounded
    // memory usage, which can be quite a big deal. So if the caller
    // asked for it there, we use it here too. (It is a bit slower,
    // but not much.)
    
    MonoNote mn(m_fixedLag > 0.5f);
    vector<MonoNote::FrameOutput> mnOut = mn.process(smoothedPitch);

    std::cerr << "mnOut size: " << mnOut.size() << std::endl;
    std::cerr << "m_pitchTrack size: " << m_pitchTrack.size() << std::endl;
    
    // turning feature into a note feature
    Feature f;
    f.hasTimestamp = true;
    f.hasDuration = true;
    f.values.clear();
        
    int onsetFrame = 0;
    bool isVoiced = 0;
    bool oldIsVoiced = 0;
    size_t nFrame = m_pitchTrack.size();

    float minNoteFrames = (m_inputSampleRate*m_pruneThresh) / m_stepSize;
    
    // the body of the loop below should be in a function/method
    // but what does it actually do??
    // * takes the result of the note tracking HMM
    // * collects contiguously pitched pitches
    // * writes a note once it notices the voiced segment has ended
    // complications:
    // * it needs a lookahead of two frames for m_level (wtf was I thinking)
    // * it needs to know the timestamp (which can be guessed from the frame no)
    // * 
    int offset = m_preciseTime == 1.0 ? m_blockSize/2 : m_blockSize/4;
    RealTime timestampOffset = Vamp::RealTime::frame2RealTime(offset, 
        lrintf(m_inputSampleRate));

    std::vector<float> notePitchTrack; // collects pitches for 1 note at a time
    for (size_t iFrame = 0; iFrame < nFrame; ++iFrame)
    {
        isVoiced = mnOut[iFrame].noteState < 3 
            && smoothedPitch[iFrame].size() > 0 
            && (iFrame >= nFrame-2
                || ((m_level[iFrame]/m_level[iFrame+2]) > m_onsetSensitivity));
        if (isVoiced && iFrame != nFrame-1)
        {
            if (oldIsVoiced == 0) // beginning of a note
            {
                onsetFrame = iFrame;
            }
            float pitch = smoothedPitch[iFrame][0].first;
            notePitchTrack.push_back(pitch); // add to the note's pitch track
        } else { // not currently voiced
            if (oldIsVoiced == 1) // end of note
            {
                if (notePitchTrack.size() >= minNoteFrames)
                {
                    std::sort(notePitchTrack.begin(), notePitchTrack.end());
                    float medianPitch = notePitchTrack[notePitchTrack.size()/2];
                    float medianFreq = 
                        std::pow(2,(medianPitch - 69) / 12) * 440;
                    f.values.clear();
                    f.values.push_back(medianFreq);
                    RealTime start = RealTime::frame2RealTime(
                        onsetFrame * m_stepSize, lrintf(m_inputSampleRate)) + 
                        timestampOffset;
                    RealTime end   = RealTime::frame2RealTime(
                            iFrame * m_stepSize, lrintf(m_inputSampleRate)) + 
                        timestampOffset;
                    f.timestamp = start;
                    f.duration = end - start;
                    fs[m_oNotes].push_back(f);
                }
                notePitchTrack.clear();
            }
        }
        oldIsVoiced = isVoiced;
    }
}

vector<pair<double, double> >
PYinVamp::salienceCandidates(const vector<pair<double, double> > &yinFreqProb,
                             float rms)
{
    vector<HarmonicSalience::Candidate> candidates = m_salience->next();
    if (m_commonF0 > 0) {
        m_salience->penaliseCommonFundamentals(candidates, m_commonF0);
    }

    // A ramp in the frame's level from `from` to `to` times the
    // low-amplitude threshold (1 if there is none)
    auto levelRamp = [&](float from, float to) {
        if (m_lowAmp <= 0) return 1.0;
        double g = (rms - from * m_lowAmp) / ((to - from) * m_lowAmp);
        return std::max(0.0, std::min(1.0, g));
    };

    // YIN's candidates' probabilities sum to its probability that the
    // frame is voiced. The part on periods longer than the HMM's lowest
    // pitch never reaches the HMM, here as in pYIN 1.2: it leaves a
    // frame of hum, or of a quiet chord's common period, unvoiced.
    double voiced = 0;
    for (size_t i = 0; i < yinFreqProb.size(); ++i) {
        if (yinFreqProb[i].first > m_pitchHmm.m_minFreq) {
            voiced += yinFreqProb[i].second;
        }
    }

    // YIN's aperiodicity stays high where two or three voices sound at
    // once, though each of them is harmonic: the salience's
    // harmonicity gives such a frame a floor of voiced probability, a
    // half at most, if it is loud enough
    const float harmonicityLow = 0.05f, harmonicityHigh = 0.12f;
    double r = (m_salience->harmonicity() - harmonicityLow) /
        (harmonicityHigh - harmonicityLow);
    r = std::max(0.0, std::min(1.0, r));
    voiced = std::max(voiced, 0.5 * r * levelRamp(0.2f, 0.4f));

    // A very quiet frame (hum, breath, a separation's residue) is not
    // voiced at all
    voiced *= levelRamp(0.05f, 0.15f);

    if (m_follower) {
        m_lastLikelihoods = m_salience->voiceLikelihoods
            (candidates, VoiceFollower::preferenceRule(m_follower->voice()));
        vector<double> gateLikelihoods = m_salience->voiceLikelihoods
            (candidates, VoiceFollower::gateRule());
        m_lastGateVoices.clear();
        for (size_t i = 0; i < candidates.size(); ++i) {
            double f = m_salience->refine(candidates[i].freq,
                                          HarmonicSalience::Params().refineHarmonics);
            VoiceFollower::GateVoice v;
            v.pitch = float(12 * std::log2(f / 440.0) + 69);
            v.likelihood = float(gateLikelihoods[i]);
            m_lastGateVoices.push_back(v);
        }
    }

    // No harmonic peak, or no voiced probability left, makes an
    // unvoiced frame: YIN's candidates would bring back the long
    // periods the salience is there to avoid
    return m_salience->candidateProbabilities(candidates, voiced);
}
