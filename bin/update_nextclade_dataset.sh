#!/bin/bash

version="0.20230412"
USAGE="
This is the script for updating the included nextclade database.
This is not part of the Cecret workflow and is kept here for 
UPHL's purposes. Also, 'nextclade' MUST be in path.

Usage:
update_nextclade_dataset.sh
"

echo "$USAGE"

if [ -z $(which nextclade) ]
then
    echo "$(date): FATAL : nextclade could not be found."
    exit 1
fi

echo "$(date): Downloading dataset from nextclade." && \
nextclade dataset get --name sars-cov-2 --output-zip sars.zip && \
echo "$(date): Moving sars.zip to data" && \
mv sars.zip data/sars.zip && \
echo "$(date): Finished."

