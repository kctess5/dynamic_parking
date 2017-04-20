#!/bin/bash

convert -delay 3 -loop 0 ./build/policy_sequence/*.png policy_animated.gif
convert -delay 3 -loop 0 ./build/sequence/*.png animated.gif
convert -delay 3 -loop 0 ./build/slices/*.png animated_slices.gif
convert -delay 3 -loop 0 ./build/policy_slices/*.png policy_slices.gif
convert -delay 3 -loop 0 ./build/useful_sequence/*.png useful_sequence.gif