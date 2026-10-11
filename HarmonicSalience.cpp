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

#include <algorithm>
#include <cmath>

using std::vector;
using std::pair;

static const int binsPerOctave = 120; // 10 cents
static const double pi = 3.14159265358979323846;

static double centsBetween(double a, double b)
{
    return 1200.0 * std::log(a / b) / std::log(2.0);
}

static bool bySalience(const HarmonicSalience::Candidate &a,
                       const HarmonicSalience::Candidate &b)
{
    return a.salience > b.salience;
}

HarmonicSalience::HarmonicSalience(float sampleRate, int blockSize, int stepSize,
                                   int centreOffset, const Params &params) :
    m_sampleRate(sampleRate),
    m_blockSize(blockSize),
    m_stepSize(stepSize),
    m_centreOffset(centreOffset),
    m_p(params),
    m_streamStart(0),
    m_streamEnd(0),
    m_pushed(0),
    m_frame(0),
    m_fft(params.fftSize),
    m_maxSalience(0),
    m_peakSum(0)
{
    int n = m_p.windowSize;
    m_window.resize(n);
    double wsum = 0;
    for (int i = 0; i < n; ++i) {
        m_window[i] = 0.5 - 0.5 * cos(2 * pi * i / n);
        wsum += m_window[i];
    }
    // so that a sinusoid of amplitude a peaks at about a
    for (int i = 0; i < n; ++i) m_window[i] *= 2.0 / wsum;

    m_fftIn.resize(m_p.fftSize, 0.0);
    m_fftOut.resize(m_p.fftSize + 2, 0.0);
    m_mag.resize(m_p.fftSize / 2 + 1, 0.0);

    // bin 0 a tolerance window below the lowest f0, so that a
    // candidate at the lowest f0 is an interior maximum
    m_gridRef = m_p.minF0 * pow(2.0, -(m_p.tolBins + 2) / binsPerOctave);
    m_nBins = int(ceil(binsPerOctave * log2(m_p.maxF0 / m_gridRef))) + int(m_p.tolBins) + 3;
    m_sal.resize(m_nBins, 0.0);
}

HarmonicSalience::~HarmonicSalience()
{
}

void
HarmonicSalience::reset()
{
    m_stream.clear();
    m_streamStart = 0;
    m_streamEnd = 0;
    m_pushed = 0;
    m_frame = 0;
    m_maxSalience = 0;
    m_peakSum = 0;
    m_peaks.clear();
}

void
HarmonicSalience::append(const float *x, int n)
{
    for (int i = 0; i < n; ++i) m_stream.push_back(x[i]);
    m_streamEnd += n;
    // keep what the windows of the frames not yet computed can reach
    long long from = centreOf(m_frame) - m_p.windowSize / 2;
    long long drop = from - m_streamStart;
    if (drop > (long long)m_stream.size() / 2 && drop > 4096) {
        drop = std::min<long long>(drop, m_stream.size());
        m_stream.erase(m_stream.begin(), m_stream.begin() + drop);
        m_streamStart += drop;
    }
}

long long
HarmonicSalience::centreOf(long long frame) const
{
    return frame * (long long)m_stepSize + m_centreOffset;
}

void
HarmonicSalience::push(const float *block)
{
    // The host's blocks overlap: only the last step's worth is new
    if (m_pushed == 0) append(block, m_blockSize);
    else if (m_stepSize <= m_blockSize) append(block + m_blockSize - m_stepSize, m_stepSize);
    else append(block, m_blockSize); // blocks with gaps between them: the gaps are lost
    ++m_pushed;
}

bool
HarmonicSalience::ready() const
{
    return m_frame < m_pushed &&
        centreOf(m_frame) + m_p.windowSize / 2 <= m_streamEnd;
}

vector<HarmonicSalience::Candidate>
HarmonicSalience::next()
{
    long long start = centreOf(m_frame) - m_p.windowSize / 2;
    ++m_frame;

    for (int i = 0; i < m_p.windowSize; ++i) {
        long long j = start + i;
        double x = (j >= m_streamStart && j < m_streamEnd) ? m_stream[j - m_streamStart] : 0.0;
        m_fftIn[i] = x * m_window[i];
    }
    m_fft.forward(m_fftIn.data(), m_fftOut.data());
    for (int k = 0; k <= m_p.fftSize / 2; ++k) {
        m_mag[k] = sqrt(m_fftOut[2*k] * m_fftOut[2*k] + m_fftOut[2*k+1] * m_fftOut[2*k+1]);
    }

    findPeaks();
    accumulate();

    // The candidates: the salience's maxima within the f0 range, each
    // placed between its grid bins by a parabola
    vector<Candidate> out;
    m_maxSalience = 0;
    for (int b = 1; b + 1 < m_nBins; ++b) {
        double s1 = m_sal[b];
        if (s1 <= 0 || !(s1 > m_sal[b-1] && s1 >= m_sal[b+1])) continue;
        double s0 = m_sal[b-1], s2 = m_sal[b+1];
        double den = s0 - 2 * s1 + s2;
        double d = den < 0 ? 0.5 * (s0 - s2) / den : 0.0;
        d = std::max(-0.5, std::min(0.5, d));
        double f = m_gridRef * pow(2.0, (b + d) / binsPerOctave);
        if (f < m_p.minF0 || f > m_p.maxF0) continue;
        m_maxSalience = std::max(m_maxSalience, s1);
        Candidate c;
        c.freq = f;
        c.salience = s1;
        out.push_back(c);
    }
    std::sort(out.begin(), out.end(), bySalience);
    if (int(out.size()) > m_p.maxCandidates) out.resize(m_p.maxCandidates);
    return out;
}

void
HarmonicSalience::findPeaks()
{
    m_peaks.clear();
    double binHz = m_sampleRate / m_p.fftSize;
    int lo = std::max(1, int(ceil(m_p.minPeakHz / binHz)));
    int hi = std::min(m_p.fftSize / 2 - 1, int(floor(m_p.maxPeakHz / binHz)));
    double top = 0;
    for (int b = lo; b <= hi; ++b) top = std::max(top, m_mag[b]);
    if (top < 1e-7) return; // digital silence or nearly
    double thr = top * pow(10.0, m_p.peakFloorDb / 20.0);
    for (int b = lo; b <= hi; ++b) {
        double m1 = m_mag[b];
        if (m1 <= thr || !(m1 > m_mag[b-1] && m1 >= m_mag[b+1])) continue;
        // a parabola through the log magnitudes places the peak
        double a0 = log(m_mag[b-1] + 1e-12), a1 = log(m1), a2 = log(m_mag[b+1] + 1e-12);
        double den = a0 - 2 * a1 + a2;
        double d = den < 0 ? 0.5 * (a0 - a2) / den : 0.0;
        d = std::max(-0.5, std::min(0.5, d));
        Peak p;
        p.freq = (b + d) * binHz;
        p.amp = exp(a1 - 0.25 * (a0 - a2) * d);
        m_peaks.push_back(p);
    }
}

void
HarmonicSalience::accumulate()
{
    std::fill(m_sal.begin(), m_sal.end(), 0.0);
    const double tol = m_p.tolBins;
    // cos^2 of the distance from pos across the window, stepping the
    // angle by a rotation rather than calling cos for every bin
    const double step = pi / 2 / tol;
    const double cstep = cos(step), sstep = sin(step);
    auto spread = [&](double pos, double weight) {
        int b0 = std::max(0, int(ceil(pos - tol)));
        int b1 = std::min(m_nBins - 1, int(floor(pos + tol)));
        if (b1 < b0) return;
        double a = step * (b0 - pos);
        double c = cos(a), s = sin(a);
        for (int b = b0; b <= b1; ++b) {
            m_sal[b] += weight * c * c;
            double cn = c * cstep - s * sstep;
            s = s * cstep + c * sstep;
            c = cn;
        }
    };
    m_peakSum = 0;
    for (const Peak &pk : m_peaks) {
        double a = pow(pk.amp, m_p.beta);
        m_peakSum += a;
        double base = binsPerOctave * log2(pk.freq / m_gridRef);
        double w = 1.0;
        for (int h = 1; h <= m_p.nHarmonics; ++h) {
            if (h > 1) w *= m_p.alpha;
            // this peak as the h-th harmonic of f0 = freq / h
            double pos = base - binsPerOctave * log2(double(h));
            // and as the (h - 1/2)-th multiple, between harmonics,
            // which is a higher f0
            double npos = base - binsPerOctave * log2(h - 0.5);
            if (npos < -tol) break; // every later h is lower still
            if (pos >= -tol && pos <= m_nBins - 1 + tol) spread(pos, w * a);
            if (h > 1 && m_p.negWeight > 0 && npos <= m_nBins - 1 + tol) {
                spread(npos, -m_p.negWeight * w * a);
            }
        }
    }
}

void
HarmonicSalience::penaliseCommonFundamentals(vector<Candidate> &candidates,
                                             double penalty) const
{
    penaliseCommonFundamentals(candidates, m_peaks,
                               1.0 * m_sampleRate / m_p.windowSize, penalty);
}

void
HarmonicSalience::penaliseCommonFundamentals(vector<Candidate> &sc,
                                             const vector<Peak> &peaks,
                                             double resolutionHz,
                                             double penalty)
{
    // The share of a candidate's salience each of the two at its
    // multiples must have, the highest multiple looked at, and how
    // far below the frame's highest peak a partial of its own counts
    static const float minShare = 0.7f;
    static const int maxMultiple = 5;
    static const double partialFloorDb = -20.0;

    int n = int(sc.size());
    if (n < 3) return;
    double topAmp = 0;
    for (const auto &pk : peaks) topAmp = std::max(topAmp, pk.amp);
    double minAmp = topAmp * std::pow(10.0, partialFloorDb / 20.0);
    // A partial within 50 cents, or the window's resolution, of freq:
    // two low voices a third apart make one spectral peak between them
    auto partialAt = [&](double freq) {
        double tolHz = std::max(freq * (std::pow(2.0, 50.0 / 1200) - 1), resolutionHz);
        for (const auto &pk : peaks) {
            if (pk.amp >= minAmp && std::fabs(pk.freq - freq) < tolHz) return true;
        }
        return false;
    };
    vector<bool> common(n, false);
    for (int i = 0; i < n; ++i) {
        int multiples = 0;
        for (int m = 2; m <= maxMultiple && multiples < 2; ++m) {
            for (int j = 0; j < n; ++j) {
                if (j == i || sc[j].salience < minShare * sc[i].salience) continue;
                if (std::fabs(centsBetween(sc[j].freq, m * sc[i].freq)) < 50) {
                    ++multiples;
                    break;
                }
            }
        }
        if (multiples < 2) continue;
        // A low voice whose formants lift its third and fourth partials
        // has strong candidates at those multiples too, but a partial
        // of its own, which a common fundamental has not
        if (partialAt(sc[i].freq)) continue;
        common[i] = true;
    }
    bool any = false;
    for (int i = 0; i < n; ++i) {
        if (common[i]) {
            sc[i].salience *= 1.0 - penalty;
            any = true;
        }
    }
    if (!any) return;
    std::stable_sort(sc.begin(), sc.end(), bySalience);
    while (!sc.empty() && sc.back().salience <= 0) sc.pop_back();
}

vector<pair<double, double> >
HarmonicSalience::candidateProbabilities(const vector<Candidate> &candidates,
                                         double voicedMass) const
{
    vector<pair<double, double> > out;
    if (candidates.empty() || voicedMass <= 0) return out;
    // candidates are sorted, highest salience first
    double smax = candidates[0].salience;
    double psum = 0;
    for (const auto &c : candidates) {
        double p = std::pow(c.salience / smax, 4.0);
        out.push_back(pair<double, double>(refine(c.freq, m_p.refineHarmonics), p));
        psum += p;
    }
    for (auto &o : out) o.second *= voicedMass / psum;
    return out;
}

double
HarmonicSalience::refine(double f0, int nHarm) const
{
    if (m_peaks.empty() || f0 <= 0) return f0;
    double tolRatio = pow(2.0, m_p.tolBins / binsPerOctave);
    double num = 0, den = 0;
    for (int h = 1; h <= nHarm; ++h) {
        double target = h * f0;
        double lo = target / tolRatio, hi = target * tolRatio;
        // the strongest peak within the tolerance of this harmonic
        auto it = std::lower_bound(m_peaks.begin(), m_peaks.end(), lo,
                                   [](const Peak &p, double f) { return p.freq < f; });
        double bestAmp = 0, bestFreq = 0;
        for (; it != m_peaks.end() && it->freq <= hi; ++it) {
            if (it->amp > bestAmp) { bestAmp = it->amp; bestFreq = it->freq; }
        }
        if (bestAmp <= 0) continue;
        double w = bestAmp * bestAmp;
        num += w * h * bestFreq;
        den += w * h * h;
    }
    if (den <= 0) return f0;
    double f = num / den;
    // a fit that moves the candidate further than its tolerance has
    // found another voice's partials
    if (fabs(1200 * log2(f / f0)) > 10 * m_p.tolBins) return f0;
    return f;
}
