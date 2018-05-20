#!/bin/bash
# Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
# See COPYING and LICENSE files at project root for more details.
make
make test
./bin/test
./bin/blitzdg

