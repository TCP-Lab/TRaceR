#!/bin/bash

# --- Global settings ----------------------------------------------------------

# Strict mode options
#set -e           # "exit-on-error" shell option
set -u           # "no-unset" shell option
#set -o pipefail  # exit on within-pipe error
#set -o errtrace  # ERR trap inherited by shell functions

# For a friendlier use of colors in Bash
red=$'\e[1;31m' # Red
grn=$'\e[1;32m' # Green
yel=$'\e[1;33m' # Yellow
blu=$'\e[1;34m' # Blue
mag=$'\e[1;35m' # Magenta
cya=$'\e[1;36m' # Cyan
end=$'\e[0m'
