#!/bin/bash

convert -delay 5 -loop 0 ./build/policy_sequence/*.png policy_animated.gif
convert -delay 5 -loop 0 ./build/sequence/*.png animated.gif