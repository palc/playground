#!/usr/bin/env bash

# First argument is a tree file to be updated

# Get the number of tree samples
count=`nw_labels $1 | wc -l`
filename=`echo $1 | sed 's/\.tre//'`
outfile="${filename}.nex"
printf "#NEXUS\n" > $outfile
printf "begin taxa;\n" >> $outfile
printf "dimensions ntax=${count};\n" >> $outfile
printf "taxlabels\n" >> $outfile
nw_labels $1 >> $outfile
printf ";\n" >> $outfile
printf "end;\n\n" >> $outfile
printf "begin trees;\n" >> $outfile
printf "\ttree tree_1 = [&R]\n" >> $outfile
cat $1 >> $outfile
printf "end;\n\n" >> $outfile
printf 'begin figtree;\n' >> $outfile
printf '\tset appearance.backgroundColorAttribute="Default";\n' >> $outfile
printf '\tset appearance.backgroundColour=#ffffff;\n' >> $outfile
printf '\tset appearance.branchColorAttribute="User selection";\n' >> $outfile
printf '\tset appearance.branchColorGradient=false;\n' >> $outfile
printf '\tset appearance.branchLineWidth=1.0;\n' >> $outfile
printf '\tset appearance.branchMinLineWidth=0.0;\n' >> $outfile
printf '\tset appearance.branchWidthAttribute="Fixed";\n' >> $outfile
printf '\tset appearance.foregroundColour=#000000;\n' >> $outfile
printf '\tset appearance.hilightingGradient=false;\n' >> $outfile
printf '\tset appearance.selectionColour=#2d3680;\n' >> $outfile
printf '\tset branchLabels.colorAttribute="User selection";\n' >> $outfile
printf '\tset branchLabels.displayAttribute="Branch times";\n' >> $outfile
printf '\tset branchLabels.fontName="sansserif";\n' >> $outfile
printf '\tset branchLabels.fontSize=8;\n' >> $outfile
printf '\tset branchLabels.fontStyle=0;\n' >> $outfile
printf '\tset branchLabels.isShown=false;\n' >> $outfile
printf '\tset branchLabels.significantDigits=4;\n' >> $outfile
printf '\tset layout.expansion=0;\n' >> $outfile
printf '\tset layout.layoutType="RADIAL";\n' >> $outfile
printf '\tset layout.zoom=0;\n' >> $outfile
printf '\tset legend.attribute=null;\n' >> $outfile
printf '\tset legend.fontSize=10.0;\n' >> $outfile
printf '\tset legend.isShown=false;\n' >> $outfile
printf '\tset legend.significantDigits=4;\n' >> $outfile
printf '\tset nodeBars.barWidth=4.0;\n' >> $outfile
printf '\tset nodeBars.displayAttribute=null;\n' >> $outfile
printf '\tset nodeBars.isShown=false;\n' >> $outfile
printf '\tset nodeLabels.colorAttribute="User selection";\n' >> $outfile
printf '\tset nodeLabels.displayAttribute="Node ages";\n' >> $outfile
printf '\tset nodeLabels.fontName="sansserif";\n' >> $outfile
printf '\tset nodeLabels.fontSize=8;\n' >> $outfile
printf '\tset nodeLabels.fontStyle=0;\n' >> $outfile
printf '\tset nodeLabels.isShown=false;\n' >> $outfile
printf '\tset nodeLabels.significantDigits=4;\n' >> $outfile
printf '\tset nodeShape.colourAttribute="User selection";\n' >> $outfile
printf '\tset nodeShape.isShown=false;\n' >> $outfile
printf '\tset nodeShape.minSize=10.0;\n' >> $outfile
printf '\tset nodeShape.scaleType=Width;\n' >> $outfile
printf '\tset nodeShape.shapeType=Circle;\n' >> $outfile
printf '\tset nodeShape.size=4.0;\n' >> $outfile
printf '\tset nodeShape.sizeAttribute="Fixed";\n' >> $outfile
printf '\tset polarLayout.alignTipLabels=false;\n' >> $outfile
printf '\tset polarLayout.angularRange=0;\n' >> $outfile
printf '\tset polarLayout.rootAngle=0;\n' >> $outfile
printf '\tset polarLayout.rootLength=100;\n' >> $outfile
printf '\tset polarLayout.showRoot=true;\n' >> $outfile
printf '\tset radialLayout.spread=0.0;\n' >> $outfile
printf '\tset rectilinearLayout.alignTipLabels=false;\n' >> $outfile
printf '\tset rectilinearLayout.curvature=0;\n' >> $outfile
printf '\tset rectilinearLayout.rootLength=100;\n' >> $outfile
printf '\tset scale.offsetAge=0.0;\n' >> $outfile
printf '\tset scale.rootAge=1.0;\n' >> $outfile
printf '\tset scale.scaleFactor=1.0;\n' >> $outfile
printf '\tset scale.scaleRoot=false;\n' >> $outfile
printf '\tset scaleAxis.automaticScale=true;\n' >> $outfile
printf '\tset scaleAxis.fontSize=8.0;\n' >> $outfile
printf '\tset scaleAxis.isShown=false;\n' >> $outfile
printf '\tset scaleAxis.lineWidth=1.0;\n' >> $outfile
printf '\tset scaleAxis.majorTicks=1.0;\n' >> $outfile
printf '\tset scaleAxis.origin=0.0;\n' >> $outfile
printf '\tset scaleAxis.reverseAxis=false;\n' >> $outfile
printf '\tset scaleAxis.showGrid=true;\n' >> $outfile
printf '\tset scaleBar.automaticScale=true;\n' >> $outfile
printf '\tset scaleBar.fontSize=10.0;\n' >> $outfile
printf '\tset scaleBar.isShown=true;\n' >> $outfile
printf '\tset scaleBar.lineWidth=1.0;\n' >> $outfile
printf '\tset scaleBar.scaleRange=0.0;\n' >> $outfile
printf '\tset tipLabels.colorAttribute="User selection";\n' >> $outfile
printf '\tset tipLabels.displayAttribute="Names";\n' >> $outfile
printf '\tset tipLabels.fontName="sansserif";\n' >> $outfile
printf '\tset tipLabels.fontSize=8;\n' >> $outfile
printf '\tset tipLabels.fontStyle=0;\n' >> $outfile
printf '\tset tipLabels.isShown=true;\n' >> $outfile
printf '\tset tipLabels.significantDigits=4;\n' >> $outfile
printf '\tset trees.order=false;\n' >> $outfile
printf '\tset trees.orderType="increasing";\n' >> $outfile
printf '\tset trees.rooting=false;\n' >> $outfile
printf '\tset trees.rootingType="User Selection";\n' >> $outfile
printf '\tset trees.transform=false;\n' >> $outfile
printf '\tset trees.transformType="cladogram";\n' >> $outfile
printf 'end;' >> $outfile

# created 2016-12-09 stuber
