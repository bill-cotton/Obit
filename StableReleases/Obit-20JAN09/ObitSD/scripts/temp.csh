# Partial processing of GBT GC C band data
echo "******** Edit Index table to remove all **********************************************"
# Flag data not part of the GC imaging
#DONE source doCbandEdit.csh
# First calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal1.fits -GAINUSE -1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal1.fits -GAINUSE -1 -FLAGVER -1
# Residual calibration was GCCDay2.fits
echo test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCstart.fits -calType Filter-SOLINT 60.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCstart.fits -calType Filter -SOLINT 60.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Convert Soln 1 to Cal table 
echo test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
# Second calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal2.fits -GAINUSE 1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal2.fits -GAINUSE 1 -FLAGVER -1
# Residual calibration pass 2
echo test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal2.fits -calType Filter -SOLINT 60.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal2.fits -calType Filter -SOLINT 60.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
# Convert Soln 2 + Cal 1  to Cal table 2
echo test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
# Third calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal3.fits -GAINUSE 2 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal3.fits -GAINUSE 2 -FLAGVER -1
# Residual calibration pass 3
echo test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal3.fits -calType Filter -SOLINT 60.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal3.fits -calType Filter -SOLINT 60.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
# Convert Soln 3 + Cal 2 to Cal table 3
echo test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal4.fits -GAINUSE 3 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal4.fits -GAINUSE 3 -FLAGVER -1
# Residual calibration pass 4
echo test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal4.fits -calType Filter -SOLINT 60.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal4.fits -calType Filter -SOLINT 60.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
# Convert Soln 4 + Cal 3 to Cal table 4
echo test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
# Fifth calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal5.fits -GAINUSE 4 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal5.fits -GAINUSE 4 -FLAGVER -1
#
# Split/calibrate with shorter time 
echo test/SplitOTF -input input/SplitOTFCAll.in
touch FITSdata/GCCbandCalOTF.fits
rm FITSdata/GCCbandCalOTF.fits
test/SplitOTF -input input/SplitOTFCAll.in
echo "******** Edit Index table to remove all **********************************************"

echo test/ResidCalOTF -input input/ResidCalOTFC2.in -image GCCcal5.fits -calType Filter -SOLINT 5.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC2.in -image GCCcal5.fits -calType Filter -SOLINT 5.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Convert Soln 1 to Cal table 1
echo test/Soln2Cal -input input/Soln2CalC2.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
test/Soln2Cal -input input/Soln2CalC2.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFC2.in -outfile \!GCCcal6.fits -GAINUSE 1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC2.in -outfile \!GCCcal6.fits -GAINUSE 1 -FLAGVER -1
#
# another calibration pass
echo test/ResidCalOTF -input input/ResidCalOTFC2.in -image GCCcal6.fits -calType Filter -SOLINT 5.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC2.in -image GCCcal6.fits -calType Filter -SOLINT 16.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
# Convert Soln 2 + Cal 1 to Cal table 2
echo test/Soln2Cal -input input/Soln2CalC2.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
test/Soln2Cal -input input/Soln2CalC2.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
# calibration Image
echo test/ImageOTF -input input/ImageOTFC2.in -outfile \!GCCcal7.fits -GAINUSE 2 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC2.in -outfile \!GCCcal7.fits -GAINUSE 2 -FLAGVER -1
#
# another calibration pass
echo test/ResidCalOTF -input input/ResidCalOTFC2.in -image GCCcal7.fits -calType Offset -SOLINT 10.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC2.in -image GCCcal7.fits -calType Offset -SOLINT 10.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
# Convert Soln 3 + Cal 2 to Cal table 3
echo test/Soln2Cal -input input/Soln2CalC2.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
test/Soln2Cal -input input/Soln2CalC2.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
# calibration Image
echo test/ImageOTF -input input/ImageOTFC2.in -outfile \!GCCcal8.fits -GAINUSE 3 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC2.in -outfile \!GCCcal8.fits -GAINUSE 3 -FLAGVER -1
#
# another calibration pass
echo test/ResidCalOTF -input input/ResidCalOTFC2.in -image GCCcal8.fits -calType Offset -SOLINT 5.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC2.in -image GCCcal8.fits -calType Offset -SOLINT 5.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
# Convert Soln 3 + Cal 2 to Cal table 3
echo test/Soln2Cal -input input/Soln2CalC2.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
test/Soln2Cal -input input/Soln2CalC2.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
# calibration Image
echo test/ImageOTF -input input/ImageOTFC2.in -outfile \!GCCcal9.fits -GAINUSE 4 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC2.in -outfile \!GCCcal9.fits -GAINUSE 4 -FLAGVER -1
