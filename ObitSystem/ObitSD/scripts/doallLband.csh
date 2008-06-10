# Do full processing of the Lband data on the Galactic Center from May 03
# Read both Harvey's and our data to OTF file
#DONE source ReadLBandBoth.csh 
# Basic (Atmospheric) calibration
test/AtmCorOTF -input input/AtmCorOTFL.in
# Split/calibrate
echo test/SplitOTF -input input/SplitOTFL.in
test/SplitOTF -input input/SplitOTFL.in


# ******************************* Non cumulative calibration
# First calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal1.fits -GAINUSE -1 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal1.fits -GAINUSE -1 -FLAGVER -1
# Residual calibration
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal1.fits -SOLINT 120.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal1.fits -SOLINT 120.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Second calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal2.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal2.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 2
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal2.fits -SOLINT 90.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal2.fits -SOLINT 60.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Third calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal3.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal3.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 3
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal3.fits -SOLINT 60.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal3.fits -SOLINT 60.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Fourth calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal4.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal4.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 4
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal4.fits -SOLINT 45.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal4.fits -SOLINT 45.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# Fifth calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal5.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal5.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 5
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal5.fits -SOLINT 30.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal5.fits -SOLINT 30.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal6.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal6.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 6
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal6.fits -SOLINT 15.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal6.fits -SOLINT 15.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal7.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal7.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 7
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal7.fits -SOLINT 15.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal7.fits -SOLINT 15.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal8.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal8.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 8
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal8.fits -SOLINT 8.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal8.fits -SOLINT 8.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal9.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal9.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 9
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal9.fits -SOLINT 5.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal9.fits -SOLINT 5.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal10.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal10.fits -GAINUSE 0 -FLAGVER -1
#
# Residual calibration pass 10
echo test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal10.fits -SOLINT 5.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
test/ResidCalOTF -input input/ResidCalOTFL.in -image GCLcal10.fits -SOLINT 5.0 -MINEL 5.0 -GAINUSE -1 -FLAGVER 1
# calibration Image
echo test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal11.fits -GAINUSE 0 -FLAGVER -1
test/ImageOTF -input input/ImageOTFL.in -outfile \!GCLcal11.fits -GAINUSE 0 -FLAGVER -1
