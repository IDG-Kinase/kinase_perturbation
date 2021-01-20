#!/bin/bash

nice -n 19 ionice -c 3 ./run_salmon_SUM159.sh
nice -n 19 ionice -c 3 ./run_salmon_SUM159_24vs48.sh
nice -n 19 ionice -c 3 ./run_salmon_Zawistowski.sh
nice -n 19 ionice -c 3 ./run_salmon_HCC1806.sh
nice -n 19 ionice -c 3 ./run_salmon_HCC1806_part2.sh
nice -n 19 ionice -c 3 ./run_salmon_HCC1806_part3.sh

nice -n 19 ionice -c 3 ./run_salmon_SUM229.sh
nice -n 19 ionice -c 3 ./run_salmon_MDAMB231.sh
nice -n 19 ionice -c 3 ./run_salmon_MDAMB468.sh
