#!/bin/bash

set -xe

if [[ "$1" == "global" ]]; then

  scale=$(echo "global")

else

  scale=$(echo "local")

fi

for range in $(cat scripts/yearrange)

  do

    start_date=$(echo $range | sed 's/|.*//')
    end_date=$(echo $range | sed 's/.*|//')

    bash -i scripts/env_download.sh $start_date $end_date $scale

done
