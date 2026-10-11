/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "YinUtil.h"

#include <boost/test/unit_test.hpp>

#include <utility>
#include <vector>

using std::pair;
using std::vector;

BOOST_AUTO_TEST_SUITE(TestYinUtil)

BOOST_AUTO_TEST_CASE(difference1)
{
    // difference always needs to be non-negative
    double in[] = { 1, -1, 1, -1, 1, -1, 1, -1 };
    double yinBuffer[4];

    YinUtil yinUtil(4);

    yinUtil.slowDifference(in, yinBuffer);

    BOOST_CHECK(yinBuffer[0] >= 0);
    BOOST_CHECK(yinBuffer[1] >= 0);
    BOOST_CHECK(yinBuffer[2] >= 0);
    BOOST_CHECK(yinBuffer[3] >= 0);

    yinUtil.fastDifference(in, yinBuffer);

    BOOST_CHECK(yinBuffer[0] >= 0);
    BOOST_CHECK(yinBuffer[1] >= 0);
    BOOST_CHECK(yinBuffer[2] >= 0);
    BOOST_CHECK(yinBuffer[3] >= 0);
}

// BOOST_AUTO_TEST_CASE(difference2)
// {
//     // difference for fast sinusoid
//     double in[] = { 1, -1, 1, -1, 1, -1, 1, -1 };
//     double yinBuffer[4];
//     
//     Yin y(8, 256);
//     
//     y.difference(in, yinBuffer);
//     
//     BOOST_CHECK_EQUAL(yinBuffer[0], 0);
//     BOOST_CHECK_EQUAL(yinBuffer[1], 16);
//     BOOST_CHECK_EQUAL(yinBuffer[2], 0);
//     BOOST_CHECK_EQUAL(yinBuffer[3], 16);
// }
// 
// BOOST_AUTO_TEST_CASE(difference3)
// {
//     // difference for contaminated fast sinusoid
//     double in[] = { 2, -1, 1, -1, 1, -1, 1, -1 };
//     double yinBuffer[4];
//     
//     Yin y(8, 256);
//     
//     y.difference(in, yinBuffer);
//     
//     BOOST_CHECK_EQUAL(yinBuffer[0], 0);
//     BOOST_CHECK_EQUAL(yinBuffer[1], 21);
//     BOOST_CHECK_EQUAL(yinBuffer[2], 1);
//     BOOST_CHECK_EQUAL(yinBuffer[3], 21);
// }
// 
// BOOST_AUTO_TEST_CASE(difference4)
// {
//     // 
//     double in[2048];
//     double yinBuffer[1024];
//     for (size_t i = 0; i < 2048; i = i + 2) {
//         in[i] = 1;
//         in[i+1] = -1;
//     }
//     Yin y(2048, 44100);
//     
//     y.difference(in, yinBuffer);
//     
//     BOOST_CHECK(yinBuffer[0] >= 0);
//     BOOST_CHECK(yinBuffer[1] >= 0);
// }
// 
// BOOST_AUTO_TEST_CASE(cumulativeDifference1)
// {
//     // test against matlab implementation
//     double yinBuffer[] = {0, 21, 1, 21};
//     
//     Yin y(8, 256);
//     
//     y.cumulativeDifference(yinBuffer);
//     
//     BOOST_CHECK_EQUAL(yinBuffer[0], 1);
//     BOOST_CHECK_EQUAL(yinBuffer[1], 1); // this is just chance!
//     BOOST_CHECK(yinBuffer[2] >= 0.0909 && yinBuffer[2] < 0.091);
//     BOOST_CHECK(yinBuffer[3] >= 1.4651 && yinBuffer[3] < 1.466);
// 
// }



// A cumulative mean normalised difference function of 1 everywhere but
// at dips of the given depth, each a few lags wide
static vector<double> dips(int size, const vector<pair<int, double> > &at)
{
    vector<double> buffer(size, 1.0);
    for (const auto &d : at) {
        buffer[d.first - 1] = (1.0 + d.second) / 2;
        buffer[d.first] = d.second;
        buffer[d.first + 1] = (1.0 + d.second) / 2;
    }
    return buffer;
}

BOOST_AUTO_TEST_CASE(betaThresholdPrior)
{
    vector<double> prior = YinUtil::betaThresholdPrior(2, 4);
    BOOST_REQUIRE_EQUAL(prior.size(), 100u);
    double sum = 0;
    size_t mode = 0;
    for (size_t i = 0; i < prior.size(); ++i) {
        sum += prior[i];
        if (prior[i] > prior[mode]) mode = i;
    }
    BOOST_CHECK_CLOSE(sum, 1.0, 1e-9);
    // Beta(2, 4) has its mode at 1/4: threshold 0.25
    BOOST_CHECK_EQUAL(mode, 24u);
}

BOOST_AUTO_TEST_CASE(orderedRuleFavoursTheFirstDip)
{
    // A voice's period at lag 100, and the deeper dip at twice it that
    // a period-doubled voice (or two voices' common period) makes:
    // every dip gives the deeper the more, the order gives the first
    YinUtil yinUtil(512);
    vector<double> buffer = dips(512, { { 100, 0.15 }, { 200, 0.02 } });
    vector<double> prior = YinUtil::betaThresholdPrior(2, 4);
    vector<double> every = yinUtil.yinProb(&buffer[0], 2);
    vector<double> ordered = yinUtil.yinProbOrdered(&buffer[0], prior, 1.0, 0.08);
    BOOST_CHECK_GT(every[200], every[100]);
    BOOST_CHECK_GT(ordered[100], ordered[200]);
    BOOST_CHECK_GT(ordered[200], 0.0);
    for (int tau = 0; tau < 512; ++tau) {
        if (tau != 100 && tau != 200) BOOST_CHECK_EQUAL(ordered[tau], 0.0);
    }
}

BOOST_AUTO_TEST_CASE(orderedRuleSharesTheThresholds)
{
    // Each threshold's probability goes to the dips below it, and 0.01
    // of it to the deepest dip if none is
    YinUtil yinUtil(512);
    vector<double> buffer = dips(512, { { 100, 0.205 } });
    vector<double> prior = YinUtil::betaThresholdPrior(2, 4);
    vector<double> p = yinUtil.yinProbOrdered(&buffer[0], prior, 1.0, 0.08);
    double expected = 0;
    for (int i = 0; i < 100; ++i) {
        expected += (0.01 + i * 0.01 > 0.205 ? 1.0 : 0.01) * prior[i];
    }
    BOOST_CHECK_CLOSE(p[100], expected, 1e-9);
}

BOOST_AUTO_TEST_CASE(orderedRuleMergesRipples)
{
    // Two dips 4% of a period apart are one valley, at its deepest
    YinUtil yinUtil(512);
    vector<double> buffer = dips(512, { { 100, 0.10 }, { 104, 0.08 } });
    vector<double> prior = YinUtil::betaThresholdPrior(2, 4);
    vector<double> merged = yinUtil.yinProbOrdered(&buffer[0], prior, 1.0, 0.08);
    BOOST_CHECK_EQUAL(merged[100], 0.0);
    BOOST_CHECK_GT(merged[104], 0.0);
    vector<double> apart = yinUtil.yinProbOrdered(&buffer[0], prior, 1.0, 0.0);
    BOOST_CHECK_GT(apart[100], apart[104]);
}

BOOST_AUTO_TEST_SUITE_END()

