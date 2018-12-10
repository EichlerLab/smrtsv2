#!/usr/bin/env bash

echo "Refused to run interpreter \"python\". SMRT-SV must call \"python2\" or \"python3\"" >&2
echo "explicitly. If this message appears from SMRT-SV, please report this as a bug." >&2

exit 1
