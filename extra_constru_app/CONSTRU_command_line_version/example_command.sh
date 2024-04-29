#!/usr/bin/env bash

./CONSTRU_command_line_version.r -m Test_part1.csv -d Test_part2.csv -r 'CYT.Tertiles' -f 'Surv( OS.TIME..yrs. , OS.Event..0.censored..1.death. ) ~ metagene' -o Test_CONSTRU.csv
