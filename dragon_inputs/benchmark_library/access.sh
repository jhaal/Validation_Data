#!/bin/sh
pos=$1/../../libraries/l_endian
if [ -f "$pos"/draglibJef2p2_v5p1.gz ] && [ ! -f "$pos"/draglibJef2p2_v5p1 ]; then
  gunzip -k "$pos"/draglibJef2p2_v5p1.gz
fi
ln -s "$pos"/draglibJef2p2_v5p1 DLIB_J2
if [ -f "$pos"/draglibendfb8r1_v5p1.gz ] && [ ! -f "$pos"/draglibendfb8r1_v5p1 ]; then
  gunzip -k "$pos"/draglibendfb8r1_v5p1.gz
fi
ln -s "$pos"/draglibendfb8r1_v5p1 DLIB_E81
ln -sf "$1"/data/VRF_MATCH_J2_R20_BU00p0_proc/*.c2m .
echo "VRF_MATCH_J2_R20_BU00p0 access script terminated"
