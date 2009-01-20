# Do full processing of the Xband data on the Galactic Center from June 03
#DONE echo source ReadXBandDay1.csh 
#DONE source ReadXBandDay1.csh 
#DONE echo source ReadXBandDay2.csh 
#DONE source ReadXBandDay2.csh 
#DONE echo source ReadXBandDay3.csh 
#DONE source ReadXBandDay3.csh 
# Basic (Atmospheric) calibration
#DONE echo test/AtmCorOTF -input input/AtmCorOTFXDay1.in
#DONE test/AtmCorOTF -input input/AtmCorOTFXDay1.in
#DONE echo test/AtmCorOTF -input input/AtmCorOTFXDay2.in
#DONE test/AtmCorOTF -input input/AtmCorOTFXDay2.in
#DONE echo test/AtmCorOTF -input input/AtmCorOTFXDay3.in
#DONE test/AtmCorOTF -input input/AtmCorOTFXDay3.in
# Split/calibrate
#DONE echo test/SplitOTF -input input/SplitOTFXDay1.in
#DONE touch FITSdata/GCXbandOTF.fits
#DONE rm FITSdata/GCXbandOTF.fits
#DONE test/SplitOTF -input input/SplitOTFXDay1.in
#DONE echo test/SplitOTF -input input/SplitOTFXDay2.in
#DONE touch FITSdata/GCXcalDay2OTF.fits
#DONE rm FITSdata/GCXcalDay2OTF.fits
#DONE test/SplitOTF -input input/SplitOTFXDay2.in
#DONE echo test/SplitOTF -input input/SplitOTFXDay3.in
#DONE touch FITSdata/GCXcalDay3OTF.fits
#DONE rm FITSdata/GCXcalDay3OTF.fits
#DONE test/SplitOTF -input input/SplitOTFXDay3.in
# Combine Day1 and Day2
#DONE echo test/ConcatOTF -input input/ConcatOTFX1.in
#DONE test/ConcatOTF -input input/ConcatOTFX1.in
#DONE echo test/ConcatOTF -input input/ConcatOTFX2.in
#DONE test/ConcatOTF -input input/ConcatOTFX2.in
#DONE echo Now edit index and target files in 
# Flag data not part of the GC imaging
#DONE source doXbandEdit.csh
# Copy pristine version of data
/bin/cp FITSdata/GCXbandOTF.save FITSdata/GCXbandOTF.fits
# First calibration Image
echo test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal1.fits -GAINUSE -1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal1.fits -GAINUSE -1 -FLAGVER -1
# Residual calibration with result of last time
echo test/ResidCalOTF -input input/ResidCalOTFX.in -image GCXcal1.fits -calType Filter -SOLINT 90. -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX.in -image GCXcal1.fits -calType Filter -SOLINT 90.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Convert Soln 1 to Cal table 
echo test/Soln2Cal -input input/Soln2CalX.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
test/Soln2Cal -input input/Soln2CalX.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
# Second calibration Image
echo test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal2.fits -GAINUSE 1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal2.fits -GAINUSE 1 -FLAGVER -1
# Residual calibration pass 2
echo test/ResidCalOTF -input input/ResidCalOTFX.in -image GCXcal2.fits -calType Filter -SOLINT 90.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX.in -image GCXcal2.fits -calType Filter -SOLINT 90.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
# Convert Soln 2 + Cal 1  to Cal table 2
echo test/Soln2Cal -input input/Soln2CalX.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
test/Soln2Cal -input input/Soln2CalX.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
# Third calibration Image
echo test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal3.fits -GAINUSE 2 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal3.fits -GAINUSE 2 -FLAGVER -1
# Residual calibration pass 3
echo test/ResidCalOTF -input input/ResidCalOTFX.in -image GCXcal3.fits --calType Filter SOLINT 45.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX.in -image GCXcal3.fits -calType Filter -SOLINT 45.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
# Convert Soln 3 + Cal 2 to Cal table 3
echo test/Soln2Cal -input input/Soln2CalX.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
test/Soln2Cal -input input/Soln2CalX.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal4.fits -GAINUSE 3 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal4.fits -GAINUSE 3 -FLAGVER -1

# Residual calibration pass 4 - Filter calibration
echo test/ResidCalOTF -input input/ResidCalOTFX.in -image GCXcal4.fits -calType Filter -SOLINT 45.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX.in -image GCXcal4.fits -calType Filter -SOLINT 45.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
# Convert Soln 4 + Cal 3 to Cal table 4
echo test/Soln2Cal -input input/Soln2CalX.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
test/Soln2Cal -input input/Soln2CalX.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
# Fifth calibration Image
echo test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal5.fits -GAINUSE 4 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX.in -outfile \!GCXcal5.fits -GAINUSE 4 -FLAGVER -1
#
# Split/calibrate with shorter time - Filter calibration
echo test/SplitOTF -input input/SplitOTFXAll.in
touch FITSdata/GCXbandCalOTF.fits
rm FITSdata/GCXbandCalOTF.fits
test/SplitOTF -input input/SplitOTFXAll.in
echo test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal5.fits -calType Filter -SOLINT 10.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal5.fits -calType Filter -SOLINT 10.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Convert Soln 1 to Cal table 1
echo test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal6.fits -GAINUSE 1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal6.fits -GAINUSE 1 -FLAGVER -1
#
# another calibration pass
echo test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal6.fits -calType Filter -SOLINT 10.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal6.fits -calType Filter -SOLINT 10.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
# Convert Soln 2 + Cal 1 to Cal table 2
echo test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal7.fits -GAINUSE 2 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal7.fits -GAINUSE 2 -FLAGVER -1
#
# another calibration pass - use boxcar
echo test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal7.fits -calType Filter -SOLINT 10.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal7.fits -calType Filter -SOLINT 10.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
# Convert Soln 3 + Cal 2 to Cal table 3
echo test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
# Fifth calibration Image
echo test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal8.fits -GAINUSE 3 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal8.fits -GAINUSE 3 -FLAGVER -1

# another calibration pass - use boxcar
echo test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal8.fits -calType Filter -SOLINT 5.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal8.fits -calType Filter -SOLINT 5.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
# Convert Soln 4 + Cal 3 to Cal table 4
echo test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
# Sixth calibration Image
echo test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal9.fits -GAINUSE 4 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal9.fits -GAINUSE 4 -FLAGVER -1

# another calibration pass - use Filter
echo test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal9.fits -calType Filter -SOLINT 5.0 -MINEL 5.0 -GAINUSE 4 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFX2.in -image GCXcal9.fits -calType Filter -SOLINT 5.0 -MINEL 5.0 -GAINUSE 4 -FLAGVER 1
# Convert Soln 5 + Cal 4 to Cal table 5
echo test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 5 -CALIN 4 -CALOUT 5
test/Soln2Cal -input input/Soln2CalX2.in -SOLNUSE 5 -CALIN 4 -CALOUT 5
# Sixth calibration Image
echo test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal10.fits -GAINUSE 5 -FLAGVER -1
test/ImageOTF -input input/ImageOTFX2.in -outfile \!GCXcal10.fits -GAINUSE 5 -FLAGVER -1


