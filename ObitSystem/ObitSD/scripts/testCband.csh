# Do full processing of the Cband data on the Galactic Center from June 03
#DONE echo source ReadCBandDay1.csh 
#DONE source ReadCBandDay1.csh 
#DONE echo source ReadCBandDay2.csh 
#DONE source ReadCBandDay2.csh 
# Basic (Atmospheric) calibration
#DONE echo test/AtmCorOTF -input input/AtmCorOTFCDay1.in
#DONE test/AtmCorOTF -input input/AtmCorOTFCDay1.in
#DONE echo test/AtmCorOTF -input input/AtmCorOTFCDay2.in
#DONE test/AtmCorOTF -input input/AtmCorOTFCDay2.in
# Split/calibrate
#DONE echo test/SplitOTF -input input/SplitOTFCDay1.in
#DONE test/SplitOTF -input input/SplitOTFCDay1.in
#DONE echo test/SplitOTF -input input/SplitOTFCDay2.in
#DONE test/SplitOTF -input input/SplitOTFCDay2.in
# Combine Day1 and Day2
#DONE echo test/ConcatOTF -input input/ConcatOTFC.in
#DONE test/ConcatOTF -input input/ConcatOTFC.in
#DONE echo Now edit index and target filed in 
# Flag data not part of the GC imaging
#DONE source doCbandEdit.csh
# First calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal1.fits -GAINUSE -1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal1.fits -GAINUSE -1 -FLAGVER -1
# Residual calibration
echo test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCDay2.fits -SOLINT 15.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCDay2.fits -SOLINT 15.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Convert Soln 1 to Cal table 
echo test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 1 -CALIN -1 -CALOUT 1
# Second calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal2.fits -GAINUSE 1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal2.fits -GAINUSE 1 -FLAGVER -1
# Residual calibration pass 2
echo test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal2.fits -SOLINT 30.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal2.fits -SOLINT 30.0 -MINEL 5.0 -GAINUSE 1 -FLAGVER 1
# Convert Soln 2 + Cal 1  to Cal table 2
echo test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 2 -CALIN 1 -CALOUT 2
# Third calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal3.fits -GAINUSE 2 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal3.fits -GAINUSE 2 -FLAGVER -1
# Residual calibration pass 3
echo test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal3.fits -SOLINT 45.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal3.fits -SOLINT 45.0 -MINEL 5.0 -GAINUSE 2 -FLAGVER 1
# Convert Soln 3 + Cal 2 to Cal table 3
echo test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 3 -CALIN 2 -CALOUT 3
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal4.fits -GAINUSE 3 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal4.fits -GAINUSE 3 -FLAGVER -1
# Residual calibration pass 4
echo test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal4.fits -SOLINT 60.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFC.in -image GCCcal4.fits -SOLINT 60.0 -MINEL 5.0 -GAINUSE 3 -FLAGVER 1
# Convert Soln 4 + Cal 3 to Cal table 4
echo test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
test/Soln2Cal -input input/Soln2CalC.in -SOLNUSE 4 -CALIN 3 -CALOUT 4
# Fifth calibration Image
echo test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal5.fits -GAINUSE 4 -FLAGVER -1
test/ImageOTF -input input/ImageOTFC.in -outfile \!GCCcal5.fits -GAINUSE 4 -FLAGVER -1
