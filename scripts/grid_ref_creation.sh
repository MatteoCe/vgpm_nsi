#!/bin/bash

set -ue

for num in $(seq 1 3)

  do
    
    if [[ $num -eq 1 ]]
      
      then

        coord=$(echo -170 90 -180 80)

      elif [[ $num -eq 2 ]]

        then

          coord=$(echo -170 -90 -180 -80)

        else

          coord=$(echo 170 -90 180 -80)

      fi

    gdal_create -of GTiff \
                -ot UInt16 \
                -a_nodata -9999 \
                -burn 0 \
                -outsize 250 250 \
                -a_srs "EPSG:4326" \
                -a_ullr $coord \
                -co COMPRESS=LZW \
                scripts/grid_ref/grid_ref_$num.tif 

  done

