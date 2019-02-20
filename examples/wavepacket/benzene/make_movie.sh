#!/bin/bash

mencoder "mf://*.png" -mf fps=10 -o movie.avi -ovc lavc -ovc lavc -lavcopts vcodec=mpeg4:mbd=1:vbitrate=10000

