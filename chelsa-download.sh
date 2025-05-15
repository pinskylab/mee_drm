#!/bin/bash

set -e
set -u
set -o pipefail

BASE_URL_TEMPLATE="https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/monthly/tas/CHELSA_tas_MM_YY_V.2.1.tif"
SHAPEFILE="data/birds/shape/window.shp"
OUTPUT_DIR="data/chelsa"

mkdir -p "$OUTPUT_DIR"
if [ ! -f "$SHAPEFILE" ]; then
    echo "Error: Shapefile not found at $SHAPEFILE"
    exit 1
fi

MONTHS=(01 02 03 04 05 06 07 08 09 10 11 12)
# YEARS=({1980..2019})
YEARS=({1982..2018})

for y in "${YEARS[@]}"; do
    for m in "${MONTHS[@]}"; do
        CURRENT_URL="${BASE_URL_TEMPLATE/MM/${m}}"
        CURRENT_URL="${CURRENT_URL/YY/${y}}"

        TEMP_DOWNLOAD_FILENAME="full_tas_${m}_${y}.tif"
        PROCESSED_FILENAME="tas_${m}_${y}.tif"

        TEMP_DOWNLOAD_PATH="${OUTPUT_DIR}/${TEMP_DOWNLOAD_FILENAME}"
        PROCESSED_FILE_PATH="${OUTPUT_DIR}/${PROCESSED_FILENAME}"

        echo "  Downloading: $CURRENT_URL to $TEMP_DOWNLOAD_PATH"

        if [ ! -f "$PROCESSED_FILE_PATH" ]; then
           if wget --quiet "$CURRENT_URL" -O "$TEMP_DOWNLOAD_PATH"; then
               echo "  Download successful."

               # Process the file with gdalwarp
               echo "  Processing: $TEMP_DOWNLOAD_PATH to $PROCESSED_FILE_PATH"
               if gdalwarp -cutline "$SHAPEFILE" -crop_to_cutline \
                           "$TEMP_DOWNLOAD_PATH" "$PROCESSED_FILE_PATH"; then
                   echo "  Processing successful."
                   echo "Removing temporary file: $TEMP_DOWNLOAD_PATH"
                   rm "$TEMP_DOWNLOAD_PATH"
               else
                   echo "  Error: gdalwarp failed for $TEMP_DOWNLOAD_PATH. Temporary file kept for inspection."
               fi
           else
               echo "  Error: wget failed to download $CURRENT_URL. Skipping this file."
           fi
           echo "----------------------------------------"
        fi
    done
    echo "  Calculating yearly average for ${y}."
    gdal_calc -A ${OUTPUT_DIR}/*_${y}.tif --outfile="${OUTPUT_DIR}/tas_${y}_avg.tif" --calc="average(A,axis=0)"
    # echo "  Calculating yearly minimum for ${y}."
    # gdal_calc -A ${OUTPUT_DIR}/*_${y}.tif --outfile="${OUTPUT_DIR}/tas_${y}_min.tif" --calc="min(A,axis=0)"
    # echo "  Calculating yearly maximum for ${y}."
    # gdal_calc -A ${OUTPUT_DIR}/*_${y}.tif --outfile="${OUTPUT_DIR}/tas_${y}_max.tif" --calc="max(A,axis=0)"
    rm ${OUTPUT_DIR}/*_${y}.tif
done
