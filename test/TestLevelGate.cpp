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

#include "LevelGate.h"

#include <boost/test/unit_test.hpp>

#include <cmath>

// An RMS level of this many dB. The gate's reference is the middle of
// a quarter-dB bin, so the levels here are: -20.125 dB is one.
static float level(double db)
{
    return float(std::pow(10.0, db / 20.0));
}

BOOST_AUTO_TEST_SUITE(TestLevelGate)

BOOST_AUTO_TEST_CASE(steadyLevelIsOpen)
{
    LevelGate gate(13, 100);
    for (int i = 0; i < 300; ++i) {
        BOOST_CHECK_EQUAL(gate.process(level(-20.125), true), 1.0);
    }
}

BOOST_AUTO_TEST_CASE(fadesBelowTheFloor)
{
    // 13 dB below the reference is the floor; the gate fades from 1
    // there to 0 6 dB further down
    LevelGate gate(13, 1000);
    for (int i = 0; i < 500; ++i) gate.process(level(-20.125), true);
    BOOST_CHECK_EQUAL(gate.process(level(-20.125 - 12), true), 1.0);
    BOOST_CHECK_CLOSE(gate.process(level(-20.125 - 16), true), 0.5, 1e-3);
    BOOST_CHECK_EQUAL(gate.process(level(-20.125 - 30), true), 0.0);
    BOOST_CHECK_EQUAL(gate.process(0.f, true), 0.0);
}

BOOST_AUTO_TEST_CASE(referenceIsTheSingingNotItsPeaks)
{
    // One frame in twenty 20 dB louder (accents, plosives) does not
    // raise the reference: it is the 90th percentile
    LevelGate gate(13, 1000);
    for (int i = 0; i < 500; ++i) {
        gate.process(level(i % 20 == 0 ? -0.125 : -20.125), true);
    }
    BOOST_CHECK_EQUAL(gate.process(level(-20.125 - 12), true), 1.0);
}

BOOST_AUTO_TEST_CASE(windowForgets)
{
    // After a loud passage longer ago than the window, a quieter one
    // is its own reference
    LevelGate gate(13, 100);
    for (int i = 0; i < 100; ++i) gate.process(level(-10.125), true);
    BOOST_CHECK_EQUAL(gate.process(level(-40.125), true), 0.0);
    for (int i = 0; i < 100; ++i) gate.process(level(-40.125), true);
    BOOST_CHECK_EQUAL(gate.process(level(-40.125), true), 1.0);
}

BOOST_AUTO_TEST_CASE(restDoesNotLowerTheReference)
{
    // Singing, then a rest of all but five frames of the window: the
    // silence is not sung, so the reference is still the singing's, and
    // the echo of the next phrase is below it
    LevelGate gate(13, 100);
    for (int i = 0; i < 100; ++i) gate.process(level(-20.125), true);
    for (int i = 0; i < 95; ++i) gate.process(level(-90.125), false);
    BOOST_CHECK_EQUAL(gate.process(level(-40.125), true), 0.0);
}

BOOST_AUTO_TEST_CASE(noSingingNoReference)
{
    // A window with nothing sung in it has nothing to be below: the
    // first frame sung after a long rest, however quiet, is open
    LevelGate gate(13, 100);
    for (int i = 0; i < 100; ++i) gate.process(level(-20.125), true);
    for (int i = 0; i < 100; ++i) gate.process(level(-90.125), false);
    BOOST_CHECK_EQUAL(gate.process(level(-40.125), true), 1.0);
}

BOOST_AUTO_TEST_CASE(unsungFramesAreGatedNotCounted)
{
    // Noise between phrases, 25 dB under the singing and nine frames in
    // ten of the window: it is gated against the singing's reference,
    // and however much of it there is, it leaves the reference where it
    // is, so that the echo after it is gated as well
    LevelGate gate(13, 100);
    for (int i = 0; i < 9; ++i) gate.process(level(-20.125), true);
    for (int i = 0; i < 90; ++i) {
        BOOST_CHECK_EQUAL(gate.process(level(-45.125), false), 0.0);
    }
    BOOST_CHECK_EQUAL(gate.process(level(-40.125), true), 0.0);
}

BOOST_AUTO_TEST_CASE(resetForgets)
{
    LevelGate gate(13, 1000);
    for (int i = 0; i < 100; ++i) gate.process(level(-10.125), true);
    gate.reset();
    BOOST_CHECK_EQUAL(gate.process(level(-40.125), true), 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
