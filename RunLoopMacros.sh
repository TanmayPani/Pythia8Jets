#!/bin/bash

root -l <<EOF
.L JetMaker.C
JetMaker t;
t.Loop();
EOF
