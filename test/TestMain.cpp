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

// The one main() of the test executable. Boost.Test's header-only
// variant is compiled here, so no Boost library has to be built or
// linked; the other test files include the ordinary header.

#define BOOST_TEST_MODULE pyin

#include <boost/test/included/unit_test.hpp>
