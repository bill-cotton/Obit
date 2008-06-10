# Do full processing of the Pband data on the Galactic Center from May 03
# Read data to OTF multi channel file
#DONE source ReadPband.csh
# Edit multi channel file
#DONE source doPbandEdit.csh
# Average in Frequency
#DONE test/SplitOTF -input input/SplitOTFP.in
#
# have to edit Array Geometry file here by hand in fv
# not any moreecho "have to edit Array Geometry file here by hand in fv"
#
# Basic (Atmospheric) calibration
test/AtmCorOTF -input input/AtmCorOTFP.in
# Convert Soln 1 to Cal table 1
test/Soln2Cal -input input/Soln2CalP1.in
# First calibration Image
test/ImageOTF -input input/ImageOTFP1.in
# Flag data not part of the GC imaging
test/EditOTF -input input/EditOTFP2.in
# Residual calibration
test/ResidCalOTF -input input/ResidCalOTFP.in
# Convert Soln 1 + Cal 1 to Cal table 2
test/Soln2Cal -input input/Soln2CalP2.in
# Second calibration Image
test/ImageOTF -input input/ImageOTFP2.in
