#!/usr/bin/env bash

grep -P "^A" $1 | cut -f 2-7 | awk '{printf "%s\t%i\t%i\t%s\t%s\t%i\t%i\n", $1,$2,$2+$6,$3,$4,$5,$6}'
