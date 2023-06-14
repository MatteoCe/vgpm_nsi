# NPP Normalized Seasonality Index calculator

This repository store the scripts used in the paper "XXX" to calculate the Normalized Seasonality Index (NSI) of the Net Primary Productivity (NPP), following<sup>1</sup>. This index is computed on the Vertically Generalized Production Model (VGPM), following codes and instructions at the Ocean Productivity website (http://sites.science.oregonstate.edu/ocean.productivity/).

The scripts download, process and calculate the VGPM model from the ancillary data indicated by the Ocean Productivity website, providing a final average NSI index, for the specified time range and decimal coordinates bounding box. Both 8day and monthly data can be processed.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Contact](#Contact)
- [References](#References)

## Installation

The two scripts, "main" and "vgpm nsi", are written in bash and R respectively, and rely on a limited number of programs and R packages. In order to use the scripts, _GDAL_<sup>2</sup> (Geospatial Data Abstraction Library), _GNU parallel_ and the R packages _tidyverse_<sup>3</sup> and _parallel_<sup>4</sup> must be installed. The scripts must be run from the bash terminal, after assuring that R can be run from it.

## Usage

The scripts can be run at once by navigating to a directory that must include both scripts and the data sub-directory, which will store the two files "manifest", with the list of http links to the ancillary data for chlorophyll, photosynthetically active radiation and sea surface temperatures, and "region_coord", a tsv file with the bounding box for the geographical region of interest. After editing both files as interested (if monthly data are required modify the manifest list including "month" instead of 8day in each link), the scripts can be run with this command:

	bash -i main.sh 8day		# change the argument to "month" if monthly data are required

The image below was obtained by running the scripts on 8day data from 2009 to 2011 and choosing the following geographic bounding box:

| region | minlon | maxlon | minlat | maxlat |
| -------- | -------- | -------- | -------- | -------- |
| antarctica | -180 | 180 | -80 | -40 |


<p align="center">
<img src="figures/map.png" style="width: 700px; height:700px">
<p>

This image was produced using QGIS<sup>5</sup> and the package Quantarctica<sup>6</sup>. The blue and red lines are two of the major oceanographic fronts, estimated from observed temperature and salinity data<sup>7</sup>.

## License

You are free to modify and use this script in your projects. If you use this script, we kindly request you to cite our paper as follows:

[Author(s)]: "Title of the Paper," Conference/Journal Name, Year. [DOI or URL](link to the paper)

## Contact

Inquiries can be forwarded to matteocecchetto@gmail.com

## References

1 - Brown, C. W., Schollaert Uz, S. & Corliss, B. H. Seasonality of oceanic primary production 		and its interannual variability from 1998 to 2007. Deep Sea Res. Part Oceanogr. Res. Pap. 90, 166–175 (2014).

2 - GDAL/OGR Contributors. GDAL/OGR Geospatial Data Abstraction Software Library. (2022).

3 - Wickham, H. et al. Welcome to the {tidyverse}. J. Open Source Softw. 4, 1686 (2019).

4 - R Core Team. R: A Language and Environment for Statistical Computing. (2022)

5 - QGIS Development Team. “QGIS Geographic Information System.” Open source geospatial foundation project, 2023. http://qgis.osgeo.org.

6 - Matsuoka, Kenichi, Anders Skoglund, George Roth, Jean de Pomereu, Huw Griffiths, Robert Headland, Brad Herried, Katsuro Katsumata, Anne Le Brocq, and Kathy Licht. “Quantarctica, an Integrated Mapping Environment for Antarctica, the Southern Ocean, and Sub-Antarctic Islands.” Environmental Modelling & Software 140 (2021): 105015.

7 - Alejandro H. Orsi, Thomas Whitworth III, and Worth D. Nowlin Jr (1995), On the meridional extent and fronts of the Antarctic Circumpolar Current., Deep-Sea Research, 42, 5, 641-673 