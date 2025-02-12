#!/bin/bash
set -eE -o functrace

# introduction: to download files from Google Drive using gdrive3
# "gdrive files download 1Tf8iHaUQG2yQ1ei_29t5PTIFDHD3sQzg  --recursive" should download a folder
# and all files inside it.
# but sometime, it ran into network error:
# Error: Failed read from stream: error reading a body from connection: connection reset

# simply return the script multiple times.

# use "gdrive files list --full-name" to list the folder ID and names
#1bmuY_BUEOIrumtM0DkwsGnfPON_eZH9y              panArctic200km_S2_SR_HARMONIZED_20220701_2022830_images         folder              2025-02-10 13:57:18
#1KcdzicIOl2aKHyzUCBrmgumP_HWIH7r9              panArctic200km_S2_SR_HARMONIZED_20210701_2021830_images         folder              2025-02-09 15:28:55
#1Tf8iHaUQG2yQ1ei_29t5PTIFDHD3sQzg              panArctic200km_S2_SR_HARMONIZED_20200701_2020830_images         folder              2025-02-08 14:11:02
#1Ot3xO1AZBl9WxAatxmiY0T7ZMc0K8SqS              panArctic200km_S2_SR_HARMONIZED_20200701_2020830_images         folder              2025-02-08 14:11:02

folderID="1bmuY_BUEOIrumtM0DkwsGnfPON_eZH9y"
folderName="panArctic200km_S2_SR_HARMONIZED_20220701_2022830_images"

logfile=${folderName}_log.txt

# Retry loop for downloading the folder
while true; do
  echo "Downloading folder '$folderName' (ID: $folderID)..."

  # Attempt to download the folder
  gdrive files download ${folderID}  --recursive --overwrite 2>&1 | tee ${logfile}

  # Extract the expected file count from the log
  EXPECTED_FILE_COUNT=$(grep -oP 'Found \K[0-9]+' ${logfile} | head -1)

  # Check if the expected file count was extracted
  if [ -z "$EXPECTED_FILE_COUNT" ]; then
    echo "Error: Could not extract the expected file count. Retrying..."
    sleep 5
    continue
  fi

  # Count the number of files in the folder
  FILE_COUNT=$(find "$folderName" -type f | wc -l)
  echo "Current file count: $FILE_COUNT (Expected: $EXPECTED_FILE_COUNT)"

  # Break the loop if the file count matches the expected count
  if [ "$FILE_COUNT" -ge "$EXPECTED_FILE_COUNT" ]; then
    echo "All $FILE_COUNT files downloaded successfully."
    break
  else
    echo "File count mismatch: Retrying..."
    sleep 300
  fi
done





