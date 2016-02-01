#!/bin/sh
autoreconf -i
if test -f isl/autogen.sh; then
	(cd isl; ./autogen.sh)
fi
if test -f polylib/autogen.sh; then
    (cd polylib; ./autogen.sh)
fi
