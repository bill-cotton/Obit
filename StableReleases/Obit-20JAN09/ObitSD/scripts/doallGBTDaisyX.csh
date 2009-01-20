# Do full processing of the GBT X band daisy test Dec 22, 2003
#DONE echo test/DCR2OTF -input input/DCR2OTFDaisyX.in -scan 2003_12_22_22:25:41
#DONE test/DCR2OTF -input input/DCR2OTFDaisyX.in -scan 2003_12_22_22:25:41
# Basic (Atmospheric) calibration
#DONE echo test/AtmCorOTF -input input/AtmCorOTFDaisyX.in
#DONE test/AtmCorOTF -input input/AtmCorOTFDaisyX.in
# Split/calibrate
#DONE echo test/SplitOTF -input input/SplitOTFDaisyX.in
#DONE touch FITSdata/GBTDaisyX2OTF.fits
#DONE rm FITSdata/GBTDaisyX2OTF.fits
#DONE test/SplitOTF -input input/SplitOTFDaisyX.in
# Copy pristine version of data
/bin/cp FITSdata/GBTDaisyX2OTF.fits.save FITSdata/GBTDaisyX2OTF.fits

# First calibration Image 
echo test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal1.fits -GAINUSE -1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal1.fits -GAINUSE -1 -FLAGVER -1
# Residual calibration with result of last time
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX.in -image GBTDaisyXcal1.fits -calType Filter -MINFLUX 0.4 -SOLINT 90. -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX.in -image GBTDaisyXcal1.fits -calType Filter -MINFLUX 0.4 -SOLINT 90.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Convert Soln 1 to Cal table 
echo test/Soln2Cal -input input/Soln2CalDaisyX.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
test/Soln2Cal -input input/Soln2CalDaisyX.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
# Second calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal2.fits -GAINUSE 1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal2.fits -GAINUSE 1 -FLAGVER -1

# Residual calibration pass 2
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX.in -image GBTDaisyXcal2.fits -calType Filter -MINFLUX 0.4 -SOLINT 90.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX.in -image GBTDaisyXcal2.fits -calType Filter -MINFLUX 0.4 -SOLINT 90.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
# Convert Soln 2 + Cal 1  to Cal table 2
echo test/Soln2Cal -input input/Soln2CalDaisyX.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
test/Soln2Cal -input input/Soln2CalDaisyX.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
# Third calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal3.fits -GAINUSE 2 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal3.fits -GAINUSE 2 -FLAGVER -1

# Residual calibration pass 3
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX.in -image GBTDaisyXcal3.fits --calType Filter -MINFLUX 0.3 -SOLINT 45.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX.in -image GBTDaisyXcal3.fits -calType Filter -MINFLUX 0.3  -SOLINT 45.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
# Convert Soln 3 + Cal 2 to Cal table 3
echo test/Soln2Cal -input input/Soln2CalDaisyX.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
test/Soln2Cal -input input/Soln2CalDaisyX.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal4.fits -GAINUSE 3 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal4.fits -GAINUSE 3 -FLAGVER -1

# Residual calibration pass 4 - Filter calibration
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX.in -image GBTDaisyXcal4.fits -calType Filter -MINFLUX 0.3 -SOLINT 45.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX.in -image GBTDaisyXcal4.fits -calType Filter -MINFLUX 0.3 -SOLINT 45.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
# Convert Soln 4 + Cal 3 to Cal table 4
echo test/Soln2Cal -input input/Soln2CalDaisyX.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
test/Soln2Cal -input input/Soln2CalDaisyX.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
# Fifth calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal5.fits -GAINUSE 4 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX.in -outfile \!GBTDaisyXcal5.fits -GAINUSE 4 -FLAGVER -1
#
# Split/calibrate with shorter time - Filter calibration
echo test/SplitOTF -input input/SplitOTFXAll.in
touch FITSdata/GBTDaisyXCalOTF.fits
rm FITSdata/GBTDaisyXCalOTF.fits
test/SplitOTF -input input/SplitOTFDaisyXAll.in
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal5.fits -calType Filter -MINFLUX 0.2 -SOLINT 10.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal5.fits -calType Filter -MINFLUX 0.2 -SOLINT 10.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Convert Soln 1 to Cal table 1
echo test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal6.fits -GAINUSE 1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal6.fits -GAINUSE 1 -FLAGVER -1
#
# another calibration pass
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal6.fits -calType Filter -MINFLUX 0.1 -SOLINT 10.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal6.fits -calType Filter -MINFLUX 0.1 -SOLINT 10.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
# Convert Soln 2 + Cal 1 to Cal table 2
echo test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal7.fits -GAINUSE 2 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal7.fits -GAINUSE 2 -FLAGVER -1
#
# another calibration pass - use boxcar
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal7.fits -calType Filter -MINFLUX 0.1 -SOLINT 10.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal7.fits -calType Filter -MINFLUX 0.1 -SOLINT 10.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
# Convert Soln 3 + Cal 2 to Cal table 3
echo test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
# Fifth calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal8.fits -GAINUSE 3 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal8.fits -GAINUSE 3 -FLAGVER -1

# another calibration pass - use boxcar
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal8.fits -calType Filter -MINFLUX 0.1 -SOLINT 5.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal8.fits -calType Filter -MINFLUX 0.1 -SOLINT 5.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
# Convert Soln 4 + Cal 3 to Cal table 4
echo test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
# Sixth calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal9.fits -GAINUSE 4 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal9.fits -GAINUSE 4 -FLAGVER -1

# another calibration pass - use Filter
echo test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal9.fits -calType Filter -MINFLUX 0.05 -SOLINT 5.0 -MINEL 5.0 -GAINUSE 4 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFDaisyX2.in -image GBTDaisyXcal9.fits -calType Filter -MINFLUX 0.05 -SOLINT 5.0 -MINEL 5.0 -GAINUSE 4 -FLAGVER 1
# Convert Soln 5 + Cal 4 to Cal table 5
echo test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 5 -CALIN 4 -CALOUT 5
test/Soln2Cal -input input/Soln2CalDaisyX2.in -SOLNUSE 5 -CALIN 4 -CALOUT 5
# Sixth calibration Image
echo test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal10.fits -GAINUSE 5 -FLAGVER -1
test/ImageOTF -input input/ImageOTFDaisyX2.in -outfile \!GBTDaisyXcal10.fits -GAINUSE 5 -FLAGVER -1


