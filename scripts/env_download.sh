#!/bin/bash

set -xue

# salinity data will be obtained from the Global Ocean Physics Reanalysis
# dataset (https://doi.org/10.48670/moi-00021) providing monthly or daily data.
# Daily data will be donwloaded and 8day mean computed to resamble the temp,
# par, chl and da variables downloaded from the Ocean Productivity website,
# which will be used to compute the NPP values. Long-term averages and ranges
# will be computed from 8day extremes and averaged in each year. Thus, groups
# of 8 day daily data hourly means from the first of january 2012 to the 31 of
# december of 2020 will be downloaded from the copernicus website.

start_date=$(date -d "$1" +%Y-%m-%d)
end_date=$(date -d "$2 +8 days" +%Y-%m-%d)

# set copernicus login and password
cop_log=$(echo "COPERNICUS_USERNAME")
cop_pass=$(echo "COPERNICUS_PASSWORD")

mkdir tmp
touch tmp/dates

dynamic_date=$(echo $start_date)

while [[ "$dynamic_date" < "$end_date" ]]; do

  previous_year=$(date -d "${dynamic_date}" +%Y)
  next_year=$(date -d "${dynamic_date} +8 days" +%Y)

  if [[ "$next_year" > "$previous_year" ]]; then

    echo $dynamic_date "T00:00:00" | tr -d ' ' >>tmp/dates

    dynamic_date=$(date -d "${next_year}-01-01" +%Y-%m-%d)

  else

    echo $dynamic_date "T00:00:00" | tr -d ' ' >>tmp/dates

    dynamic_date=$(date -d "${dynamic_date} +8 days" +%Y-%m-%d)

  fi

done

# copernicus data on salinity will be downloaded using copernicus-marine CLI
# framework

mkdir tmp/raw
mkdir tmp/processed
mkdir tmp/processed/week
mkdir tmp/processed/year
mkdir variables/

nlines=$(seq 1 $(cat tmp/dates | wc -l) | head -n -1 | tail -n 1)

for line in $(seq 1 $nlines); do

  nextline=$(echo $line + 1 | bc)

  begin_date=$(cat tmp/dates | sed "${line}q;d")
  final_date=$(date -d "$(cat tmp/dates | sed "${nextline}q;d") 1 day ago" +%Y-%m-%d)

  name=$(paste -d "_" <(echo $begin_date | sed 's/T.*//') <(echo $final_date | sed 's/T.*//'))

  copernicusmarine subset \
    --username ${cop_log} \
    --password ${cop_pass} \
    --force-download \
    --dataset-id cmems_mod_glo_phy_my_0.083deg_P1D-m \
    --dataset-version 202311 \
    --variable so \
    --start-datetime ${begin_date} \
    --end-datetime ${final_date} \
    --minimum-longitude -180 \
    --maximum-longitude 180 \
    --minimum-latitude -80 \
    --maximum-latitude 90 \
    --minimum-depth 0 \
    --maximum-depth 0.49402499198913574 \
    --output-filename tmp/raw/${name}

  numberOfBands=$(gdalinfo tmp/raw/${name}.nc | grep 'Band' | wc -l)

  echo "Translating" $name "week.."

  seq 1 $numberOfBands | tr ' ' '\n' |
    parallel -j $numberOfBands 'gdal_translate -a_srs epsg:4326 -b {} -sds tmp/raw/'"$name"'.nc -of XYZ -co "COLUMN_SEPARATOR=\t" tmp/processed/week/'"$name"'_{}'

  for i in $(seq 1 $numberOfBands); do

    next_name=$(paste -d "_" <(echo $name) <(echo $i))

    if [ $i -eq 1 ]; then

      mv tmp/processed/week/${next_name} tmp/processed/week/${name}

    else

      paste -d "\t" <(cat tmp/processed/week/${name}) <(cat tmp/processed/week/${next_name} | cut -f 3) >tmp/processed/week/tmpout && mv -f tmp/processed/week/tmpout tmp/processed/week/${name}

      rm -f tmp/processed/week/${next_name}

      if [ $i -eq $numberOfBands ]; then
        
        # Turn copernicus daily data to 8.day period and append to file with 
        # other 8.day periods
        Rscript scripts/env_process.R 1 $name

        rm -f tmp/processed/week/${name}

      else

        continue

      fi

    fi

  done

  rm -f tmp/raw/${name}.nc

  previous_year=$(date -d "$(echo $begin_date | sed 's/T.*//')" +%Y)
  next_year=$(date -d "$(echo $begin_date | sed 's/T.*//') +8 days" +%Y)

  if [[ "$next_year" > "$previous_year" ]]; then

    # Calculate range and mean values for all 8.day periods in a full year.
    Rscript scripts/env_process.R 2 $previous_year

    rm -f tmp/processed/year/${previous_year}

  else

    continue

  fi

  if [ $line -eq $nlines ]; then

    # this is disabled as it will be done at the end of all downloads
    # remove # in following line and comment "continue"

    continue

    #first_year=$(date -d "${start_date}" +%Y)
    #final_year=$(date -d "$(echo $begin_date | sed 's/T.*//')" +%Y)

    #year_range=$(paste -d "_" <(echo $first_year) <(echo $final_year))

    #Rscript scripts/env_process.R 3 $year_range

    #rm -f variables/sss_final

  else

    continue

  fi

done

# now we can proceed with the Ocean Productivity data, to donaload ancillary
# data on Sea Surface Temperature, Chlorophyll, Photosynthetically Available
# Radiation and Diffure attenuation coefficient, in order to get their long-term
# averages and ranges for the period 2012-2020, calculate and get the same
# statistics for the Net Primary Productivity and calculate the Normalized
# Seasonality index. In order to optimize the procedure, data on 8day means of
# all these parameters must be downloaded and processed altogether.

# In order to do that a "netrc" file must be created in the home directory with
# the nasa's profile username and password.

http_address_1=$(echo "http://orca.science.oregonstate.edu/data/2x4/8day/")
http_address_2=$(echo ".modis.r2022/hdf/")

# it will start by downloading the tar files in the http_manifest file
# corresponding to a single year and the environmental variables sst, chl and
# par:

touch tmp/years

dynamic_date=$(echo $start_date)

dynamic_year=$(date -d "${dynamic_date}" +%Y)
final_year=$(date -d "${end_date} +8 days" +%Y)

while [[ "$dynamic_year" < "$final_year" ]]; do

  echo $dynamic_year >>tmp/years

  dynamic_year=$(date -d "${dynamic_year}-01-01 +1 years" +%Y)

done

# now the list of years can be used in a for loop

for year in $(cat tmp/years); do

  touch tmp/http_manifest

  for var in $(echo sst chl par k490); do

    paste -d "" <(echo $http_address_1) <(echo $var) <(echo $http_address_2) <(echo $var) <(echo ".m.") <(echo $year) <(echo ".tar") >>tmp/http_manifest

  done

  correct_download_number=$(echo 184)
  download_number=$(echo 0)

  while [[ "$download_number" < "$correct_download_number" ]]; do

    wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --no-check-certificate --content-disposition -P tmp/raw/ -i tmp/http_manifest

    rm -f tmp/http_manifest

    for utar in $(ls tmp/raw/*tar); do

      tar -xvf ${utar} -C tmp/raw/

      rm -f ${utar}

    done

    download_number=$(ls tmp/raw/* | tr '\t' '\n' | wc -l)

  done

  for week in $(ls tmp/raw/* | sed 's/\.hdf.*//' | sed 's/.*\.//' | sort | uniq); do

    begin_date=$(date -d "${year}-01-01 +$(echo $week | cut -c 5-8 | paste - <(echo - 1) | bc) days" +%Y-%m-%d)

    year_b=$(date -d "${begin_date}" +%Y)
    year_f=$(date -d "${begin_date} +7 days" +%Y)

    if [[ "$year_f" > "$year_b" ]]; then

      final_date=$(date -d "${year_f}-01-01 -1 days" +%Y-%m-%d)

    else

      final_date=$(date -d "${begin_date} +7 days" +%Y-%m-%d)

    fi

    name=$(paste -d "_" <(echo $begin_date) <(echo $final_date))

    gzip -d $(ls -D tmp/raw/* | grep $week)

    echo sst chl par k490 | tr ' ' '\n' |
      parallel -j 3 'gdal_translate -a_srs epsg:4087 -sds tmp/raw/{}.'"$week"'.hdf -of XYZ -co "COLUMN_SEPARATOR=\t" tmp/processed/week/'"$name"'_{}'

    rm -f tmp/raw/*.hdf

    # Calculate Net Primary Productivity and append 8.day data to year file
    Rscript scripts/env_process.R 4 $name

    rm -f tmp/processed/week/${name}*

  done
  
  # Get Normalized Seasonality Index, total NPP for a year and mean and range 
  # values of all parameters.

  Rscript scripts/env_process.R 5 $year

  rm -f tmp/processed/year/${year}*

done

first_year=$(date -d "${start_date}" +%Y)
final_year=$(date -d "$(echo $begin_date | sed 's/T.*//')" +%Y)

year_range=$(paste -d "_" <(echo $first_year) <(echo $final_year))

if [ ! -f scripts/timeframes ]; then

  echo $1 "|" $2 | sed 's/ //g' >scripts/timeframes

fi

for times in $(cat scripts/timeframes); do

  start_year=$(date -d "$(echo $times | sed 's/|.*//')" +%Y)
  end_year=$(date -d "$(echo $times | sed 's/.*|//') 1 day ago" +%Y)

  year_range=$(echo $start_year "_" $end_year | sed 's/ //g')

  # Calculate long-term averages of all statistics for the entire year
  Rscript scripts/env_process.R 3 $year_range

  # if included as a third argument to the script, the environmental variables
  # will be interpolated in the entire globe, using 225 regions as different areas
  # for the interpolations, thus increasing speed.

  if [[ "$3" == "global" ]]; then

    rm -f variables/*_final

    # choose number of regions
    nregs=$(echo 225)

    # create extent of regions for interpolation of data. The value given
    # to regions_generator.r will only be the approximate number of regions
    # as it will be rounded to numbers whose square root give an integer
    xmin=$(Rscript scripts/regions_generator.r $nregs xmin)
    xmax=$(Rscript scripts/regions_generator.r $nregs xmax)
    ymin=$(Rscript scripts/regions_generator.r $nregs ymin)
    ymax=$(Rscript scripts/regions_generator.r $nregs ymax)

    mkdir variables/shp
    mkdir variables/shp/regions
    mkdir variables/tif

    for reg in $(seq 1 $(echo $xmin | tr ' ' '\n' | wc -l)); do

      mkdir variables/shp/regions/region_${reg}

      cat scripts/support_files/sample.geojson |
        sed "s/xmin/$(echo $xmin | cut -d ' ' -f $reg)/" |
        sed "s/xmax/$(echo $xmax | cut -d ' ' -f $reg)/" |
        sed "s/ymin/$(echo $ymin | cut -d ' ' -f $reg)/" |
        sed "s/ymax/$(echo $ymax | cut -d ' ' -f $reg)/" >variables/shp/regions/region_${reg}/points.geojson

      ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:4326 variables/shp/regions/region_${reg}/extent_${reg}.shp \
        variables/shp/regions/region_${reg}/points.geojson

    done

    # reset correct number of regions
    nregs=$(ls variables/shp/regions/ | sed 's/region_//' | wc -l)

    python3 scripts/gdal_scripts/ogrmerge.py \
      -single \
      -f 'ESRI Shapefile' \
      -o variables/shp/regions/extent_general \
      -nln extent_general \
      variables/shp/regions/*/extent_*.shp

    for file in $(ls -d variables/* | grep -v shp); do

      var_name=$(echo $file | sed 's/.*\///' | sed 's/\_.*//')

      year_range=$(echo $file | sed 's/.*\///' | sed 's/'${var_name}'\_//')

      num_stats=$(echo $(cat $file | head -n 1 | tr '\t' '\n' | wc -l) - 2 | bc)

      for cols in $(seq 1 $num_stats); do

        num_col=$(echo $cols + 2 | bc)

        stats_name=$(cat $file | head -n 1 | tr '\t' '\n' | sed "${num_col}q;d")

        # create csv for each statistics of each variable and remove NA values
        cat $file | cut -f 1,2,$num_col | grep -v "NA" | tr '\t' ',' >${var_name}_${stats_name}_${year_range}.csv

        # start the creation of a shapefile
        mkdir variables/shp/${var_name}_${stats_name}
        ogr2ogr -f "ESRI Shapefile" variables/shp/${var_name}_${stats_name}/${var_name}_${stats_name}_${year_range}.dbf ${var_name}_${stats_name}_${year_range}.csv

        # edit vrt sample file
        cat scripts/support_files/sample.vrt |
          sed 's/OGRVRTLayer name=""/OGRVRTLayer name="'${var_name}'_'${stats_name}'_'${year_range}'"/' |
          sed 's/filepathhere/'${var_name}'_'${stats_name}'_'${year_range}'\.csv/' |
          sed 's/Field name=""/Field name="'${stats_name}'"/' >variables/shp/${var_name}_${stats_name}/${var_name}_${stats_name}_${year_range}.vrt

        # create final shapefile
        ogr2ogr -f "ESRI Shapefile" variables/shp/${var_name}_${stats_name}/${var_name}_${stats_name}_${year_range}.shp \
          variables/shp/${var_name}_${stats_name}/${var_name}_${stats_name}_${year_range}.vrt

        rm -f ${var_name}_${stats_name}_${year_range}.csv

        for reg in $(ls -d variables/shp/regions/* | grep -v general | sed 's/.*_//'); do

          mkdir variables/shp/${var_name}_${stats_name}/clip_${reg}

        done

        seq 1 $nregs | parallel -j 16 'python3 scripts/gdal_scripts/ogr_layer_algebra.py CLIP \
              -input_ds variables/shp/'"$var_name"'_'"$stats_name"'/'"$var_name"'_'"$stats_name"'_'"$year_range"'.shp \
              -method_ds variables/shp/regions/region_{}/extent_{}.shp \
              -output_ds variables/shp/'"$var_name"'_'"$stats_name"'/clip_{} \
              -output_lyr '"$var_name"'_'"$stats_name"'_'"$year_range"'_clipped_{}'

        set +e

        for reg in $(ls -d variables/shp/regions/* | grep -v general | sed 's/.*_//'); do

          if [[ $(ogrinfo -al -so variables/shp/${var_name}_${stats_name}/clip_$reg/${var_name}_${stats_name}_${year_range}_clipped_$reg.shp | grep "Feature Count:" | sed 's/.* //') -lt 4 ]]; then

            continue

          else

            saga_cmd statistics_kriging 0 \
              -POINTS variables/shp/${var_name}_${stats_name}/clip_$reg/${var_name}_${stats_name}_${year_range}_clipped_$reg.shp \
              -VAR_MODEL "2.20114 + 2.75005 * x" \
              -TARGET_USER_SIZE 0.04 \
              -SEARCH_RANGE 0 \
              -SEARCH_RADIUS 0.65 \
              -SEARCH_POINTS_ALL 1 \
              -SEARCH_POINTS_MIN 16 \
              -SEARCH_POINTS_MAX 20 \
              -PREDICTION variables/tif/${var_name}_${stats_name}_${year_range}_reg_$reg

          fi

        done

        seq 1 $nregs | parallel -j $nregs 'gdal_translate variables/tif/'"$var_name"'_'"$stats_name"'_'"$year_range"'_reg_{}.sdat \
                variables/tif/'"$var_name"'_'"$stats_name"'_'"$year_range"'_reg_{}.tif'

        set -e

        rm -f variables/tif/*.mgrd variables/tif/*.prj variables/tif/*.sdat.aux variables/tif/*.sdat* variables/tif/*.sgrd

        python3 scripts/gdal_scripts/gdal_merge.py \
          -o variables/tif/${var_name}_${stats_name}.tif \
          variables/tif/${var_name}_${stats_name}_${year_range}_reg_*.tif

        rm -f variables/tif/*_reg_*

      done

    done

  # else, if global is not specified as a third argument, the interpolation will
  # be performed only on specific regions, corresponding to the paper's dataset
  # locations

  else

    # temporarily move the final files in tmp/
    mv -f variables/*_final tmp/

    # Now proceed with the interpolation, which will be performed on separate
    # regions corresponding to the location of the samples included in the paper's
    # dataset.

    # create extent of regions for interpolation of data
    xmin=$(echo 2 11 23 32 37 39.5 26 18 8 -6 -5 11.5 13 -160 137.5 161.5)
    xmax=$(echo 7 16 28 38 42 44.5 31 23 13 -1 0 16.5 18 -155 142.5 166.5)
    ymin=$(echo 40 42 33 25 18 14 41 53.5 53.5 46 41.5 65.5 76 19 -69 -77)
    ymax=$(echo 45 47 38 31 23 19 46 58.5 60.5 51 46.5 70.5 81 24 -64 -72)

    mkdir variables/shp
    mkdir variables/shp/regions
    mkdir variables/tif

    for reg in $(seq 1 $(echo $xmin | tr ' ' '\n' | wc -l)); do

      mkdir variables/shp/regions/region_${reg}

      cat scripts/support_files/sample.geojson |
        sed "s/xmin/$(echo $xmin | cut -d ' ' -f $reg)/" |
        sed "s/xmax/$(echo $xmax | cut -d ' ' -f $reg)/" |
        sed "s/ymin/$(echo $ymin | cut -d ' ' -f $reg)/" |
        sed "s/ymax/$(echo $ymax | cut -d ' ' -f $reg)/" >variables/shp/regions/region_${reg}/points.geojson

      ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:4326 variables/shp/regions/region_${reg}/extent_${reg}.shp \
        variables/shp/regions/region_${reg}/points.geojson
      
      # extract the GSHHS
      tar -xf scripts/shapefiles_global/GSHHS.tar.xz -C scripts/shapefiles_global/

      python3 scripts/gdal_scripts/ogr_layer_algebra.py CLIP \
        -input_ds scripts/shapefiles_global/GSHHS.shp \
        -method_ds variables/shp/regions/region_${reg}/extent_${reg}.shp \
        -output_ds variables/shp/regions/region_${reg}/GSHHS_${reg} \
        -output_lyr GSHHS_${reg}

      python3 scripts/gdal_scripts/ogr_layer_algebra.py ERASE \
        -input_ds variables/shp/regions/region_${reg}/extent_${reg}.shp \
        -method_ds variables/shp/regions/region_${reg}/GSHHS_${reg}/GSHHS_${reg}.shp \
        -output_ds variables/shp/regions/region_${reg}/clip_${reg} \
        -output_lyr clip_${reg}

    done

    python3 scripts/gdal_scripts/ogrmerge.py \
      -single \
      -f 'ESRI Shapefile' \
      -o variables/shp/regions/extent_general \
      -nln extent_general \
      variables/shp/regions/*/extent_*.shp

    for file in $(ls -d variables/* | grep -v shp); do

      var_name=$(echo $file | sed 's/.*\///' | sed 's/\_.*//')

      year_range=$(echo $file | sed 's/.*\///' | sed 's/'${var_name}'\_//')

      num_stats=$(echo $(cat $file | head -n 1 | tr '\t' '\n' | wc -l) - 2 | bc)

      for cols in $(seq 1 $num_stats); do

        num_col=$(echo $cols + 2 | bc)

        stats_name=$(cat $file | head -n 1 | tr '\t' '\n' | sed "${num_col}q;d")

        # create csv for each statistics of each variable and remove NA values
        cat $file | cut -f 1,2,$num_col | grep -v "NA" | tr '\t' ',' >${var_name}_${stats_name}_${year_range}.csv

        # start the creation of a shapefile
        mkdir variables/shp/${var_name}_${stats_name}
        ogr2ogr -f "ESRI Shapefile" variables/shp/${var_name}_${stats_name}/${var_name}_${stats_name}_${year_range}.dbf ${var_name}_${stats_name}_${year_range}.csv

        # edit vrt sample file
        cat scripts/support_files/sample.vrt |
          sed 's/OGRVRTLayer name=""/OGRVRTLayer name="'${var_name}'_'${stats_name}'_'${year_range}'"/' |
          sed 's/filepathhere/'${var_name}'_'${stats_name}'_'${year_range}'\.csv/' |
          sed 's/Field name=""/Field name="'${stats_name}'"/' >variables/shp/${var_name}_${stats_name}/${var_name}_${stats_name}_${year_range}.vrt

        # create final shapefile
        ogr2ogr -f "ESRI Shapefile" variables/shp/${var_name}_${stats_name}/${var_name}_${stats_name}_${year_range}.shp \
          variables/shp/${var_name}_${stats_name}/${var_name}_${stats_name}_${year_range}.vrt

        rm -f ${var_name}_${stats_name}_${year_range}.csv

        for reg in $(ls -d variables/shp/regions/* | grep -v general | sed 's/.*_//'); do

          mkdir variables/shp/${var_name}_${stats_name}/clip_${reg}

        done

        seq 1 16 | parallel -j 16 'python3 scripts/gdal_scripts/ogr_layer_algebra.py CLIP \
                  -input_ds variables/shp/'"$var_name"'_'"$stats_name"'/'"$var_name"'_'"$stats_name"'_'"$year_range"'.shp \
                  -method_ds variables/shp/regions/region_{}/clip_{}/clip_{}.shp \
                  -output_ds variables/shp/'"$var_name"'_'"$stats_name"'/clip_{} \
                  -output_lyr '"$var_name"'_'"$stats_name"'_'"$year_range"'_clipped_{}'

        set +e

        seq 1 16 | parallel -j 16 'saga_cmd  statistics_kriging 0 \
                    -POINTS variables/shp/'"$var_name"'_'"$stats_name"'/clip_{}/'"$var_name"'_'"$stats_name"'_'"$year_range"'_clipped_{}.shp \
                    -VAR_MODEL "2.20114 + 2.75005 * x" \
                    -TARGET_USER_SIZE 0.04 \
                    -PREDICTION variables/tif/'"$var_name"'_'"$stats_name"'_'"$year_range"'_reg_{}'

        set -e

        seq 1 16 | parallel -j 16 'gdal_translate variables/tif/'"$var_name"'_'"$stats_name"'_'"$year_range"'_reg_{}.sdat \
                    variables/tif/'"$var_name"'_'"$stats_name"'_'"$year_range"'_reg_{}.tif'

        rm -f variables/tif/*.mgrd variables/tif/*.prj variables/tif/*.sdat.aux variables/tif/*.sdat* variables/tif/*.sgrd

        bash -i scripts/grid_ref_creation.sh

        python3 scripts/gdal_scripts/gdal_merge.py \
          -o variables/tif/${var_name}_${stats_name}.tif \
          variables/tif/${var_name}_${stats_name}_${year_range}_reg_*.tif scripts/grid_ref/*.tif

        rm -f variables/tif/*_reg_*

      done

    done

    mkdir $times
    mv variables/tif $times/
    rm -fr variables/*
    mv -f tmp/*_final variables/

  fi

done
