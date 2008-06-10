C     $Id$
LOCAL INCLUDE 'INPUT.INC'
C                                       Declarations for inputs
      INCLUDE 'INCS:PAOOF.INC'
      INCLUDE 'INCS:PUVD.INC'
      INTEGER NPARMS
      PARAMETER (NPARMS=34)
      INTEGER AVTYPE(NPARMS), AVDIM(2,NPARMS)
      CHARACTER AVNAME(NPARMS)*8
LOCAL END
LOCAL INCLUDE 'INPUTDATA.INC'
C                                       DATA statments defining input
C                                       parameters.
C                      1        2         3          4 
      DATA AVNAME /'INNAME', 'INCLASS', 'INSEQ', 'INDISK', 
C            5          6         7            8
     *   'SOURCES', 'CALCODE', 'TIMERANG', 'DOSTOKES', 
C            9          10       11       12        13
     *   'DOCALIB', 'GAINUSE', 'DOBAND', 'BPVER', 'SMOOTH', 
C          14          15   
     *   'FLAGVER', 'DOPOL',
C           16          17        18          19        20
     *   'OUTNAME', 'OUTCLASS', 'OUT2CLAS', 'OUTSEQ', 'OUTDISK', 
C            21         22
     *   'CELLSIZE', 'IMSIZE',
C           23         24        25         26
     *   'NFIELD', 'FLDSIZE', 'RASHIFT', 'DECSHIFT', 
C           27         28         29    
     *   'UVTAPER', 'UVRANGE', 'ROBUST',
C          30        31        32        33
     *   'UVWTFN', 'UVSIZE', 'UVBOX', 'UVBXFN',
C           34
     *   'BADDISK'/
C                    1       2       3       4 
      DATA AVTYPE /OOACAR, OOACAR, OOAINT, OOAINT, 
C            5       6      7      8
     *    OOACAR, OOACAR, OOARE, OOACAR, 
C            9    10      11      12     13
     *   OOARE, OOAINT, OOAINT, OOAINT, OOARE,
C           14     15
     *   OOAINT, OOAINT, 
C          16     17      18       19      20
     *   OOACAR, OOACAR, OOACAR, OOAINT, OOAINT,
C          21     22
     *   OOARE, OOAINT,
C          23      24       25     26
     *   OOAINT,  OOAINT, OOARE, OOARE, 
C          27      28     29
     *   OOARE,  OOARE, OOARE, 
C          30      31      32     33 
     *   OOACAR,  OOARE,  OOAINT, OOAINT,
C          34
     *   OOAINT/
C                   1    2    3    4 
      DATA AVDIM /12,1, 6,1, 1,1, 1,1, 
C           5    6    7    8
     *   16,30, 4,1, 8,1, 4,1, 
C         9   10    11  12    13
     *   1,1, 1,1, 1,1, 1,1, 3,1,
C        14    15
     *   1,1, 1,1, 
C        16     17   18   19   20
     *   12,1, 6,1, 6,1, 1,1, 1,1, 
C        21   22
     *   2,1, 2,1,
C        23      24        25       26
     *   1,1, 2,MAXAFL, MAXAFL,1, MAXAFL,1,
C        27   28   29 
     *   2,1, 2,1, 1,1, 
C        30   31   32    33
     *   2,1, 2,1, 1,1, 1,1,
C          34
     *   10,1/
LOCAL END
      PROGRAM OBTST
C-----------------------------------------------------------------------
C! Testbed for Obit/AIPS connection
C# Task AP Imaging calibration
C-----------------------------------------------------------------------
C;  Copyright (C) 2003
C;  Associated Universities, Inc. Washington DC, USA.
C;
C;  This program is free software; you can redistribute it and/or
C;  modify it under the terms of the GNU General Public License as
C;  published by the Free Software Foundation; either version 2 of
C;  the License, or (at your option) any later version.
C;
C;  This program is distributed in the hope that it will be useful,
C;  but WITHOUT ANY WARRANTY; without even the implied warranty of
C;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C;  GNU General Public License for more details.
C;
C;  You should have received a copy of the GNU General Public
C;  License along with this program; if not, write to the Free
C;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
C;  MA 02139, USA.
C;
C;  Correspondence concerning this software should be addressed as follows:
C;         Internet email: bcotton@nrao.edu
C;         Postal address: W. D. Cotton
C;                         National Radio Astronomy Observatory
C;                         520 Edgemont Road
C;                         Charlottesville, VA 22903-2475 USA
C-----------------------------------------------------------------------
C   Obit/AIPS testbed
C-----------------------------------------------------------------------
      INCLUDE 'INCS:PUVD.INC'
      INCLUDE 'INCS:PSTD.INC'
      CHARACTER PRGM*6, INPUT*32, UVMULT*32,  UVSCR*32,
     *   IMAGE(MAXAFL)*32, BEAM(MAXAFL)*32
      INTEGER  IRET, BUFF1(256), NFIELD, CHAN, NCHAN, IMCHAN, I, OFF
      HOLLERITH NAME1(8), NAME2(8), ARR1(8*MAXAFL), ARR2(8*MAXAFL)
      INCLUDE 'INCS:DFIL.INC'
      INCLUDE 'INCS:DMSG.INC'
      INCLUDE 'INCS:DDCH.INC'
      INCLUDE 'INCS:DHDR.INC'
      DATA PRGM /'OBTST '/
C-----------------------------------------------------------------------
C                                       Startup
      CALL OBTSIN (PRGM, NFIELD, INPUT, UVMULT, UVSCR, IMAGE, BEAM,
     *   IRET)
      IF (IRET.GT.0) GO TO 995
C                                       Start Obit
      CALL OBITSU (IRET)
      IF (IRET.NE.0) GO TO 990
C                                       Calibrate to scratch file
C                                       Convert Strings to HOLLERITH for
C                                       C
      CALL CHR2H (32, UVMULT, 1, NAME1)
      CALL CHR2H (32, UVSCR, 1, NAME2)
      CALL OBUVCP(NAME1, NAME2, IRET)
      IF (IRET.NE.0) GO TO 990
C                                       Uniform weight
      CALL CHR2H (32, UVSCR, 1, NAME1)
      CALL CHR2H (32, IMAGE(1), 1, NAME2)
      CALL OBUFWT (NAME1, NAME2, IRET)
      IF (IRET.NE.0) GO TO 990
C                                       Images
      CALL CHR2H (32, UVSCR, 1, NAME1)
      OFF = 1
      DO 20 I = 1,NFIELD
         CALL CHR2H (32, IMAGE(I), OFF, ARR1)
         CALL CHR2H (32, BEAM(I), OFF, ARR2)
         OFF = OFF + 32
 20      CONTINUE
      CALL OBIMUV (NAME1, 0, NFIELD, IMAGE, BEAM, .TRUE., .TRUE., CHAN, 
     *   NCHAN, IMCHAN, IRET)
      IF (IRET.NE.0) GO TO 990
 990  CONTINUE
C                                       Shutdown Obit
      CALL OBSHTX
C                                       Close down files, etc.
 995  CALL DIE (IRET, BUFF1)
C
 999  STOP
      END
      SUBROUTINE OBTSIN (PRGN, NFIELD, INPUT, UVMULT, UVSCR, IMAGE,
     *   BEAM, IERR) 
C-----------------------------------------------------------------------
C   OBTSIN gets input parameters for OBTST and creates the input uvdata
C   object. 
C   Set some defaults.
C   Inputs:
C      PRGN    C*6  Program name
C   Output:
C      NFIELD  I    Number of fields
C      INPUT   C*?  Task inputs object
C      UVMULT  C*?  Input multisource uv data object.
C      UVSCR   C*?  Scratch uv data object.
C      IMAGE   C(*)*? Array of Image objects
C      BEAM    C(*)*? Array of Beam objects
C      IERR    I    Error code: 0 => ok
C-----------------------------------------------------------------------
      CHARACTER PRGN*6, INPUT*32, UVMULT*32, UVSCR*32, 
     *   IMAGE(*)*32, BEAM(*)*32
      INTEGER   NFIELD, IERR
C
      INCLUDE 'INCS:PUVD.INC'
      INCLUDE 'INCS:PSTD.INC'
      INTEGER   NKEY,NKEY2
C                                       NKEY=no. adverbs to copy to
C                                       UVMULT
      PARAMETER (NKEY=15)
      INTEGER   DIM(7), TYPE, JERR, SUOFF, FGVER, IMSI(2),
     *   IMSIZE(2,MAXFLD)
C                                       NKEY2=no. adverbs to copy to
C                                       IMAGE,BEAM
      PARAMETER (NKEY2=10)
      INTEGER   FLDSIZ(2,MAXAFL), MAXNX, MAXNY, MINNX, MINNY,
     *   IDUM, NSUBA, SUBA, I, DOCAL
      LOGICAL   ISMULT, LTEMP
      REAL  RASH(MAXAFL), DECSH(MAXAFL)
      CHARACTER INK(NKEY)*8, OUTK(NKEY)*32, 
     *   INK2(NKEY2)*8, OUTK2(NKEY2)*32,  
     *   STOKES*4, FGTAB*32, CDUM*1,
     *   ANTAB*32,INAME*12, ONAME*12, CLASS*6,
     *   CLTYPE*3 

      INCLUDE 'INCS:DFIL.INC'
      INCLUDE 'INCS:DMSG.INC'
      INCLUDE 'INPUT.INC'
      INCLUDE 'INPUTDATA.INC'
C                                       Adverbs for UVMULT
C                    1        2         3          4 
      DATA INK /'INNAME', 'INCLASS', 'INSEQ', 'INDISK', 
C             5           6          7          8        9
     *   'SOURCES', 'TIMERANG', 'DOSTOKES', 'DOCALIB', 'GAINUSE',
C            10        11          12  
     *   'FLAGVER', 'UVRANGE', 'DOPOL',
C           13        14       15
     *   'DOBAND', 'BPVER', 'SMOOTH'/
C                    1       2        3       4 
      DATA OUTK /'NAME', 'CLASS', 'IMSEQ', 'DISK', 
C             5           6          7          8        9
     *   'SOURCES', 'TIMERANG', 'DOSTOKES', 'DOCAIBL', 'GAINUSE',
C            10        11          12  
     *   'FLAGVER', 'UVRANGE', 'DOPOL',
C           13        14       15
     *   'DOBAND', 'BPVER', 'SMOOTH'/
C                                       Adverbs for IMAGE, BEAM
C                     1         2          3         4  
      DATA INK2 /'OUTSEQ', 'OUTDISK', 'UVTAPER', 'ROBUST',
C           5         6         7        8          9          10
     *   'UVWTFN', 'UVSIZE', 'UVBOX', 'UVBXFN', 'CELLSIZE', 'IMSIZE'/
C                    1       2       3          4  
      DATA OUTK2 /'IMSEQ', 'DISK', 'UVTAPER', 'ROBUST',
C           5         6         7        8          9          10
     *   'UVWTFN', 'UVSIZE', 'UVBOX', 'UVBXFN', 'CELLSIZE', 'IMSIZE'/
C-----------------------------------------------------------------------
C                                       Startup
      INPUT = 'Input'
      CALL AV2INP (PRGN, NPARMS, AVNAME, AVTYPE, AVDIM, INPUT, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       BADDISK
      CALL OGET (INPUT, 'BADDISK', TYPE, DIM, IBAD, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Create UVMULT
      UVMULT = 'Input multisource UVdata'
      CALL CREATE (UVMULT, 'UVDATA', IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Copy adverbs to object
      CALL IN2OBJ (INPUT, NKEY, INK, OUTK, UVMULT, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Open UVMULT to be sure it's OK.
      CALL OOPEN (UVMULT, 'READ', IERR)
      IF (IERR.NE.0) GO TO 999
      CALL OCLOSE (UVMULT, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       DOCALIB 
      CALL OGET (INPUT, 'DOCALIB', TYPE, DIM, DOCAL, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
      LTEMP = DOCAL .GT. 0
      CALL OPUT (UVMULT, 'DOWTCL', OOALOG, DIM, LTEMP, CDUM,
     *   IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Stokes' parameters
      CALL OGET (INPUT, 'DOSTOKES', TYPE, DIM, IDUM, STOKES, IERR)
      IF (IERR.NE.0) GO TO 999
      IF (STOKES.EQ.'    ') STOKES = 'I'
      CALL OPUT (INPUT, 'DOSTOKES', TYPE, DIM, IDUM, STOKES, IERR)
      IF (IERR.NE.0) GO TO 999
      DIM(1) = LEN (STOKES)
      DIM(2) = 1
      CALL OPUT (UVMULT, 'STOKES', OOACAR, DIM, IDUM, STOKES,
     *   IERR) 
      IF (IERR.NE.0) GO TO 999
C                                       Number of subarrays
      ANTAB = 'Temp AN for OBTSIN'
      CALL UV2TAB (UVMULT, ANTAB, 'AN', 1, IERR)
      IF (IERR.NE.0) GO TO 999
      CALL TBLHIV (ANTAB, NSUBA, IERR)
      IF (IERR.NE.0) GO TO 999
      NSUBA = MAX (1, NSUBA)
      CALL TABDES (ANTAB, IERR)
      IF (IERR.NE.0) GO TO 999
      SUBA = 1
      DIM(1) = 1
      DIM(2) = 1
      CALL OPUT (UVMULT, 'SUBARR', OOAINT, DIM, SUBA, CDUM,
     *   IERR)
      IF (IERR.NE.0) GO TO 999
      CALL OPUT (UVMULT, 'NSUBARR', OOAINT, DIM, NSUBA, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Channel increment 1
      DIM(1) = 1
      DIM(2) = 1
      CALL OUVPUT (UVMULT, 'CHINC', OOAINT, DIM, 1, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Is it a multisource file.
      CALL UVDFND (UVMULT , 1, 'SOURCE', SUOFF, JERR)
      ISMULT = ((JERR.EQ.0) .AND. (SUOFF.GT.0))
      DIM(1) = 1
      DIM(2) = 1
      CALL OPUT (UVMULT, 'ISMULTI', OOALOG, DIM, ISMULT, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Check defaults
C                                       Get number of fields (min 0)
      CALL OGET (INPUT, 'NFIELD', TYPE, DIM, NFIELD, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
      IF (NFIELD.LE.0) NFIELD = 1
      CALL OPUT (INPUT, 'NFIELD', TYPE, DIM, NFIELD, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Default IMSIZE (256)
      CALL OGET (INPUT, 'IMSIZE', TYPE, DIM, IMSI, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
      IF (IMSI(1).LE.0.0) IMSI(1) = 256
      IF (IMSI(2).LE.0.0) IMSI(2) = 256
      CALL OPUT (INPUT, 'IMSIZE', TYPE, DIM, IMSI, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       FLDSIZE - largely this is set
C                                       automatically. 
      CALL OGET (INPUT, 'FLDSIZE', TYPE, DIM, FLDSIZ, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
      CALL OPUT (INPUT, 'FLDSIZE', TYPE, DIM, FLDSIZ, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       Fill in image sizes
      DO 10 I = 1,NFIELD
         IMSIZE(1,I) = MAX (IMSI(1), (FLDSIZ(1,I)))
         IMSIZE(2,I) = MAX (IMSI(2), (FLDSIZ(2,I)))
 10   CONTINUE
C                                       Default FLAGVER
      CALL OGET (INPUT, 'FLAGVER', TYPE, DIM, FGVER, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
      IF (FGVER.LE.0.0) THEN
C                                       Find highest FG table.
         FGTAB = 'Temporary FG table for OBTSIN'
         FGVER = 1
         CALL UV2TAB (UVMULT, FGTAB, 'FG', FGVER, IERR)
         IF (IERR.NE.0) GO TO 999
         CALL TBLHIV (FGTAB, FGVER, IERR)
         IF (IERR.NE.0) GO TO 999
         CALL TABDES (FGTAB, IERR)
         IF (IERR.NE.0) GO TO 999
         END IF
      CALL OPUT (INPUT, 'FLAGVER', TYPE, DIM, FGVER, CDUM, IERR)
      IF (IERR.NE.0) GO TO 999
C                                       UVSCR
      UVSCR = 'uv data scratch file'
      CALL CREATE (UVSCR, 'UVDATA', IERR)
      IF (IERR.NE.0) GO TO 999
C                                       clone from UVDATA
      CALL OUVSCR (UVSCR, UVMULT, 0, .FALSE., IERR)
      IF (IERR.NE.0) GO TO 999
C                                 check max field number
      IF (NFIELD .GT. MAXFLD) THEN
         MSGTXT = 'TOO MANY FIELDS REQUESTED'
         IERR = 8
         CALL MSGWRT (8)
         GO TO 999
      END IF
C                                       Create Images and Beams
C                                       Name 
      CALL OGET (INPUT, 'INNAME', TYPE, DIM, IDUM, INAME, IERR)
      IF (IERR.NE.0) GO TO 999
      CALL OGET (INPUT, 'OUTNAME', TYPE, DIM, IDUM, ONAME, IERR)
      IF (IERR.NE.0) GO TO 999
      IF (ONAME .EQ. '    ') ONAME = INAME
C                                       Loop over fields
      DO 100 I = 1,NFIELD
C                                       Create Beam
         WRITE (BEAM(I),1000) ONAME(1:8),I
         CALL CREATE (BEAM(I), 'IMAGE', IERR)
         IF (IERR.NE.0) GO TO 990
C                                       Set name
         DIM(1) = LEN(ONAME)
         DIM(2) = 1
         CALL OPUT (BEAM(I), 'NAME', OOACAR, DIM, IDUM, ONAME, IERR)
         IF (IERR.NE.0) GO TO 990
         CLTYPE = 'BM'
         WRITE (CLASS,1040) CLTYPE, I
         DIM(1) = LEN(CLASS)
         CALL OPUT (BEAM(I), 'CLASS', OOACAR, DIM, IDUM, CLASS, IERR)
         IF (IERR.NE.0) GO TO 990
         DIM(1) = 2
         CALL OPUT (BEAM(I), 'IMSIZE', OOAINT, DIM, IMSIZE(1,I), CDUM,
     *      IERR) 
         IF (IERR.NE.0) GO TO 990
C                                       Create CLEAN image
         WRITE (IMAGE(I),1002) ONAME(1:8),I
         CALL CREATE (IMAGE(I), 'IMAGE', IERR)
         IF (IERR.NE.0) GO TO 990
C                                       Set name
         DIM(1) = LEN(ONAME)
         DIM(2) = 1
         CALL OPUT (IMAGE(I), 'NAME', OOACAR, DIM, IDUM, ONAME, IERR)
         IF (IERR.NE.0) GO TO 990
         CLTYPE = 'IM'
         WRITE (CLASS,1040) CLTYPE, I
         DIM(1) = LEN(CLASS)
         CALL OPUT (IMAGE(I), 'CLASS', OOACAR, DIM, IDUM, CLASS, IERR)
         IF (IERR.NE.0) GO TO 990
C                                         First image gets all sizes
         IF (I.EQ.1) THEN
            DIM(1) = 2
            DIM(2) = NFIELD
         ELSE
            DIM(1) = 2
            DIM(2) = 1
         END IF
         CALL OPUT (IMAGE(I), 'IMSIZE', OOAINT, DIM, IMSIZE(1,I), CDUM,
     *      IERR) 
         IF (IERR.NE.0) GO TO 990
         DIM(2) = 1
C                                       Copy adverbs to object
         CALL IN2OBJ (INPUT, NKEY2, INK2, OUTK2, IMAGE(I), IERR)
         IF (IERR.NE.0) GO TO 990
         CALL IN2OBJ (INPUT, NKEY2, INK2, OUTK2, BEAM(I), IERR)
         IF (IERR.NE.0) GO TO 990
 100  CONTINUE
C                                       Save number of fields on
C                                       first field
      DIM(1) = 1
      DIM(2) = 1
      CALL OPUT (IMAGE(1), 'NFIELD', OOAINT, DIM, NFIELD, CDUM, IERR)
      IF (IERR.NE.0) GO TO 990
C                                       Get Shifts
      CALL OGET (INPUT, 'RASHIFT', TYPE, DIM, RASH, CDUM, IERR)
      IF (IERR.NE.0) GO TO 990
      CALL OGET (INPUT, 'DECSHIFT', TYPE, DIM, DECSH, CDUM, IERR)
      IF (IERR.NE.0) GO TO 990
C                                       Copy shifts to first field
      DIM(1) = MAXFLD
      DIM(2) = 1
      CALL OPUT (IMAGE(1), 'RASHIFT', OOARE, DIM, RASH, CDUM, IERR)
      IF (IERR.NE.0) GO TO 990
      CALL OPUT (IMAGE(1), 'DECSHIFT', OOARE, DIM, DECSH, CDUM, IERR) 
      IF (IERR.NE.0) GO TO 990
      GO TO 999
C                                       Error
 990  MSGTXT = 'ERROR Creating objects '
      CALL MSGWRT (8)
C
 999     RETURN
C-----------------------------------------------------------------------
 1000 FORMAT (A,' BEAM field number ',I3.3)
 1002 FORMAT (A,' CLEAN field number ',I3.3)
 1040 FORMAT (A3,I3.3)
      END
      SUBROUTINE OBITSU (IERR)
C-----------------------------------------------------------------------
C   Initialize Obit system using values in COMMON
C   Output:
C      IERR    I    Error code: 0 => ok
C-----------------------------------------------------------------------
      INTEGER IERR
C
      HOLLERITH PGM(2), AIPSD(20*32), FITSD(1)
      INTEGER   I, OFF, LAST
      CHARACTER DIRNAM*128, DDXX*4, EHX*4
      INCLUDE 'INCS:DMSG.INC'
      INCLUDE 'INCS:DDCH.INC'
C-----------------------------------------------------------------------
      IERR = 0
      CALL RFILL (20*32, HBLANK, AIPSD)
      CALL RFILL (1, HBLANK, FITSD)
C                                       Convert strings to HOLLERITH for
C                                       C
      CALL CHR2H (6, TSKNAM, 1, PGM)
C                                       AIPS Disk directory numbers
      OFF = 1
      DO 10 I=1,NVOL
C                                       Disk number as EHEX
         CALL ZEHEX (I, 4, EHX)
C                                       Set name of environment
C                                       variable
         DDXX = 'DA0' // EHX(4:4)
C                                       Get String
         CALL ZTRLOG (4, DDXX, 127, DIRNAM, LAST, IERR)
         IF (IERR.NE.0) GO TO 990
C                                       Add terminal '/'
         DIRNAM(LAST+1:LAST+1) = '/'
C                                       Pack into Hollerith for C
         CALL CHR2H (128, DIRNAM, OFF, AIPSD)
         OFF = OFF + 128
 10      CONTINUE
C                                       Initialize
      CALL OBINTX (PGM, 6, NPOPS, NLUSER,
     *   NVOL, AIPSD, 128,
     *   0, FITSD, 0,
     *   .TRUE., .FALSE., IERR)
      IF (IERR.NE.0) GO TO 990
      GO TO 999
C                                       Error
 990  MSGTXT = 'ERROR Initializing Obit '
      CALL MSGWRT (8)
C
 999     RETURN
         END
