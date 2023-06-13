#!/bin/bash

set -xue

# This bash script can be run from the terminal by providing an argument that must be either 8day or month. As the manifest
# file only reports the links to the 8day files, if not changed this will need to be run exclusively with the "8day"
# argument. Future versions of the script will allow to choose between the two and automatically download the chosen type
# of files.
# to use the script, run:
# bash -i main.sh 8day

# Following the indications of the Ocean Productivity website (https://sites.science.oregonstate.edu/ocean.productivity/),
# the manifest file with the http links to the tar files of 8day at ~9km resolution at the equator of the ancillary data for
# the calculation of the VGPM of net primary productivity was created, from 2003 to 2014.
# This file can be used by wget to download all the wanted files. 

# Although not strictly necessary, in unix a "netrc" file can be created in the home directory with the nasa's profile
# username and password.

# echo "machine urs.earthdata.nasa.gov login matteocecchetto@gmail.com password Scar2017" > ~/.netrc ; > ~/.urs_cookies
# chmod  0600 ~/.netrc

# Download all MODIS ancillary data for chlorophyll (chl), photosynthetically available radiation (par) and sea surface
# temperature (sst) for each year (from 2003 to 2014)
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --no-check-certificate --content-disposition -i data/manifest

# move the downloaded files into a new directory "tarfiles"
mkdir tarfiles
mv *tar tarfiles/

# extract the content of each tar file, being the hdf files for each 8day period
for utar in $(ls tarfiles/*tar); do tar -xvf ${utar}; done

# store the codes for the ancillary data files into a variable. These codes are "chl", "sst" and "par"
parameters=$(ls *.gz | cut -c 1-3 | uniq)

# create a directory for temporary files. This is done to easily eliminate the files not necessary after each step, thus
# avoiding to occupy too much disk storage, as the ancillary data represent the entire world

mkdir processed
mkdir processed/tmp

# create an empty file that will store the translation from numeric day (e.g. 2003009) to a date format (e.g. 2003-01-09)
touch processed/tmp/dating

# The argument provided by running the script will be stored as a variable called type, allowing easy recognition throughout
# the script text
type=$(echo $1)

# store the files' names (which are the numeric data format, e.g. 2003009)
ls *.gz | sed -e 's/.hdf.gz//' | cut -c 5- > processed/filelist

# translate and process the hdf files. Initially, each file name, corresponding to a date number format, will be translated
# to a date format and stored into a list, in the dating file. The translation will be performed differently if the files 

for i in $(cat processed/filelist);

do
	if [[ "$1" == "month" ]]
	then
	
	y=$(echo ${i} | cut -c 1-4)
	d=$(echo ${i} | cut -c 5-7)
	
	dating=$(date -d "${y}-1-1 +${d} days" +%Y-%m)
	
	echo $dating > processed/tmp/dating
	
	elif [[ "$1" == "8day" ]]
	then
	
		d=$(echo ${i} | cut -c 5-7)

		if [[ $d -eq 361 ]]
		then
			
			d1=$(echo ${i}-1 | cut -c 5-9 | bc)
			d2=$(echo 001)
			
			y1=$(echo ${i} | cut -c 1-4)
	
			date1=$(date -d "${y1}-1-1 +${d1} days" +%Y-%m-%d)
			date2=$(date -d "${date1} +8 days" +%Y-%m-%d)
			
			dating=$(echo ${date1}_${date2})
			
			echo $dating > processed/tmp/dating
			
		else
		
			d1=$(echo ${i}-1 | cut -c 5-9 | bc)
			d2=$(echo ${d1}+8 | bc)
			
			y=$(echo ${i} | cut -c 1-4)
			
			date1=$(date -d "${y}-1-1 +${d1} days" +%Y-%m-%d)
			date2=$(date -d "${y}-1-1 +${d2} days" +%Y-%m-%d)
			
			dating=$(echo ${date1}_${date2})
			
			echo $dating > processed/tmp/dating
		
		fi

	fi

# each file in the main directory will be moved to the "processed" directory.
mv *${i}* processed/

# extract the file.
gzip -d processed/*.gz

# for each date (8day period) 3 jobs, one for each parameter (chl, par and sst), will be run in parallel and the hdf file
# converted to a tsv file, which will have three colums, one for latitude, one for longitude and another for the values of
# the parameter echoed in the piped parallellized gdal_translate command

echo ${parameters} | tr ' ' '\n' |\
	parallel -j 3 'gdal_translate -a_srs epsg:4087 -sds processed/{}.'"$i"'.hdf -of XYZ -co "COLUMN_SEPARATOR=\t" processed/'"$i"'_{}'

# The hdf file can be deleted now, in order to save space
rm -f processed/*.hdf

# Now, for each 8day period file, there are three files, one for chl, one for par and another for sst. In order to process
# them in R, it is necessary to aggregate columnwise the three different values at each lat/lon point. This will create
# a single file, corresponding to an 8day period, with all the information needed for the calculation of the VGPM net
# primary productivity model.
param=$(ls processed/* | tr ' ' '\n' | grep $i | cut -d '_' -f 2)

	for k in $(echo ${param} | tr ' ' '\n');
		
		do
	
		linenum=$(echo ${param} | tr ' ' '\n' | grep -n ${k} | cut -d ':' -f1)
		maxlinenum=$(echo ${param} | tr ' ' '\n' | wc -l)
		
		if [[ $linenum -eq $maxlinenum ]]
		then
			
			rm -f processed/${i}_${k}
			break
		else
		
			if [[ $linenum -eq 1 ]]
			then
				
				nextlinenum=$(echo "$linenum+1" | bc)
				nextparam=$(echo $param | tr ' ' '\n' | sed ''"$nextlinenum"'!d' -)
				
				paste <(cat processed/${i}_${k}) <(cat processed/${i}_${nextparam} | cut -f 3) > processed/${i}_${linenum}
			else
				
				prevlinenum=$(echo "$linenum-1" | bc)
				nextlinenum=$(echo "$linenum+1" | bc)
				nextparam=$(echo $param | tr ' ' '\n' | sed ''"$nextlinenum"'!d' -)
				
				paste <(cat processed/${i}_${prevlinenum}) <(cat processed/${i}_${nextparam} | cut -f 3) > processed/${i}_${linenum}
				
				rm -f processed/${i}_${prevlinenum}
				
			fi
		fi
		
		# remove the original, single values files
		rm -f processed/${i}_${k}

	done

# rename the files leaving only the annotation for the 8day period
rename 's/_[0-9]*/_tot/' processed/${i}_*
rm -f processed/${i}_[0-9]*

# store all dates in a variable
dating=$(cat processed/tmp/dating)

# append the dates to each line of the newly created files, thus obtaining a new variable that refers to it
sed "s|$|\t${dating}|" processed/${i}_tot > processed/${i}_tot_date
rm -f processed/${i}_tot

# remove lines with no data of any variable (-9999 values)
sed '/-9999/d' processed/${i}_tot_date > processed/${i}_tot_date_clean
rm -f processed/${i}_tot_date

# create header giving the names of the columns on each file
echo longitude latitude ${param} date | tr ' ' '\t' | cat - processed/${i}_tot_date_clean > processed/${i}_tot_date_clean_head
rm -f processed/${i}_tot_date_clean

# Now, calculate the VGPM model using an R script. This script works in two parts, the first for the calculation of the
# VGPM model of Net Primary Productivity, the second by calculating the seasonality index of the NPP. This can be performed
# by providing two arguments, 1 (first part) or 2 (second part).

Rscript vgpm_nsi.R 1 $i $type

# remove previous files, to save space
rm -f processed/${i}_tot_date_clean_head

# set a check to inspect at which point the analyses are at, by looking at the line number of the 8day file at the filelist
# file, which is used by the loop to iterate on the files. If the processed 8day file is the last of the list, the loop can
# stop, and the second part of the R script "vgpm_nsi.R" can proceed
firstlooplinenum=$(cat processed/filelist | sort | uniq | grep -n ${i} | cut -d ':' -f1)
firstloopmaxlinenum=$(cat processed/filelist | sort | uniq | wc -l)
		
	if [[ $firstlooplinenum -eq $firstloopmaxlinenum ]]
	then
		break

	else
		continue
	fi

done

# create a directory for the final file
mkdir processed/final
# remove all files in the temporary directory
rm -f -r processed/tmp/

# This will run the scond part of the R script, calculating the average NSI index for the year range chosen
Rscript vgpm_nsi.R 2 $type


