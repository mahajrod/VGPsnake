#!/usr/bin/env bash

grep -P "^S" $1 | cuf -f 2,5 | sed 's/rd:i://'
