#!/bin/sh

# Use to delete unnecessary files for back from MiSeq run.
# Working directory must be MiSeq run directory. Example: 140604_M00963_0131_000000000-A8VEU

rm *_Read*.txt
rm -rf ./Config
rm -rf ./Images &
rm -rf ./Logs &
rm -rf ./Recipe
rm -rf ./Thumbnail_Images &

rm -rf ./Data/RTALogs &
rm -rf ./Data/TileStatus &
rm ./Data/FractionUsedLog.txt
rm ./Data/ImageSize.dat

rm -rf ./Data/Intensities/L001 &
rm -rf ./Data/Intensities/Offsets &
rm ./Data/Intensities/config.xml
rm ./Data/Intensities/RTAConfiguration.xml

#rm -rf ./Data/Intensities/BaseCalls/Alignment &
rm -rf ./Data/Intensities/BaseCalls/L001 &
rm -rf ./Data/Intensities/BaseCalls/Matrix &
rm -rf ./Data/Intensities/BaseCalls/Phasing &

rm ./Data/Intensities/BaseCalls/config.xml

wait

echo "Deletions have finshed!!!"

#
#  Created by Stuber, Tod P - APHIS on 2014-06-06.
#


##########
