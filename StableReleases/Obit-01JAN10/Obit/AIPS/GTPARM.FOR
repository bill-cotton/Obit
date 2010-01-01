      SUBROUTINE GTPARM (NAME, NPARMS, RQUICK, RPARM, SCRTCH, IERR)
C-----------------------------------------------------------------------
C! starts tasks, getting parameters and task ID number, does accounting
C# Utility
C-----------------------------------------------------------------------
C;  Copyright (C) 1995-1996, 1998
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
C;  Correspondence concerning AIPS should be addressed as follows:
C;         Internet email: aipsmail@nrao.edu.
C;         Postal address: AIPS Project Office
C;                         National Radio Astronomy Observatory
C;                         520 Edgemont Road
C;                         Charlottesville, VA 22903-2475 USA
C-----------------------------------------------------------------------
C   GTPARM obtains the activator task number, obtains the transmitted
C   parameters, initializes the message common, and outputs the message
C   'task NAME begins'.  It also handles startup accounting.
C   Inputs:
C      NAME    C*6        Task name
C      NPARMS  I          Number of real variables wanted
C   Outputs:
C      RQUICK  L          T => release POPs as soon as possible
C                         F => wait until you have finished
C      RPARM   R(NPARMS)  Parameters received
C      SCRTCH  I(256)     Scratch buffer
C      IERR    I          Error code: 0 -> ok
C                                     1 -> initiator not found
C                                     2 -> disk troubles
C                                     3 -> initiator zeroed
C-----------------------------------------------------------------------
      CHARACTER NAME*6
      INTEGER   NPARMS, SCRTCH(256), IERR
      REAL      RPARM(*)
      LOGICAL   RQUICK
C
      INTEGER   FINDTD, ISCR, I, IER, J, LUNTD, NBL, JBL
      REAL      RSCR
      LOGICAL   T, F
      CHARACTER PHNAME*48, PTYPE*2, TNAME*6, OPT*4
      INCLUDE 'INCS:DMSG.INC'
      INCLUDE 'INCS:DDCH.INC'
      EQUIVALENCE (RSCR, ISCR)
      DATA LUNTD, PTYPE /13, 'TD'/
      DATA T, F /.TRUE.,.FALSE./
C-----------------------------------------------------------------------
      NBL = 4
C                                       set initial values
      IERR = 0
      MSGCNT = -1
      TSKNAM = NAME
      NPOPS = 1
      NLUSER = 1
      MSGSUP = 0
C                                       open task data file
      CALL ZPHFIL (PTYPE, 1, 0, 4, PHNAME, IER)
      CALL ZOPEN (LUNTD, FINDTD, 1, PHNAME, F, T, T, IER)
      IF (IER.EQ.0) GO TO 10
         IERR = 2
         GO TO 999
C                                       read directory record
 10   CALL ZFIO ('READ', LUNTD, FINDTD, 1, SCRTCH, IER)
      IF (IER.NE.0) GO TO 980
      NBL = MAX (4, SCRTCH(256))
C                                       locate initiator
C Hack this screws up big time when linked with modern software 
C restrict to POPS 1.
CHack      CALL WHOAMI (NAME, NPOPS, IER)
      NPOPS = 1
C
      J = 5 * (NPOPS-1)
      IF (IER.NE.0) GO TO 20
         CALL H2CHR (6, 1, SCRTCH(J+1), TNAME)
         IF (NAME.EQ.TNAME) GO TO 30
 20   IERR = 1
      IF ((SCRTCH(J+1).EQ.0) .AND. (SCRTCH(J+2).EQ.0)) IERR = 3
      GO TO 990
C                                       Clear task name from file
 30   SCRTCH(J+1) = 0
      SCRTCH(J+2) = 0
      CALL ZFIO ('WRIT', LUNTD, FINDTD, 1, SCRTCH, IER)
C                                       set values, read pars
      JBL = NBL*(NPOPS-1) + 2
      CALL ZFIO ('READ', LUNTD, FINDTD, JBL, SCRTCH, IER)
      IF (IER.NE.0) GO TO 980
C                                       User number
      NLUSER = SCRTCH(1)
C                                       Get TV and TK assignments
      NTVDEV = SCRTCH(2)
      NTKDEV = SCRTCH(3)
C                                       Get Message Kill parameter
      MSGKIL = SCRTCH(4)
C                                       Is this batch job
      ISBTCH = SCRTCH(5)
C                                       Is debug needed?
      DBGAIP = SCRTCH(6)
C                                       quick return
      J = 9
      ISCR = SCRTCH(J)
      RQUICK = RSCR.LE.0.0
      IF (NPOPS.GT.NINTRN+1) ISBTCH = 32000
      RQUICK = RQUICK .AND. (NPOPS.LE.NINTRN) .AND. (ISBTCH.NE.32000)
C                                       Get version ID
      J = J + 1
      CALL H2CHR (4, 1, SCRTCH(J), VERNAM)
C                                       any machine-dependent parms
      CALL ZSETUP
C                                       Get task parms
      IF (NPARMS.LE.0) GO TO 60
         DO 50 I = 1,NPARMS
            J = J + 1
            IF (J.LE.256) GO TO 45
               J = 1
               JBL = JBL + 1
               CALL ZFIO ('READ', LUNTD, FINDTD, JBL, SCRTCH, IER)
               IF (IER.NE.0) GO TO 980
 45         ISCR = SCRTCH(J)
            RPARM(I) = RSCR
 50         CONTINUE
C                                       close file
 60   CALL ZCLOSE (LUNTD, FINDTD, IER)
C                                       Start accounting
      CALL ACOUNT (1)
C                                       message: start
      WRITE (MSGTXT,1060) TSKNAM, RLSNAM
      CALL MSGWRT (5)
C                                       clean out TC file
      CALL GTTELL (-1, OPT, RPARM, SCRTCH, IER)
      GO TO 999
C                                       Errors
 980  IERR = 2
 990  CALL ZCLOSE (LUNTD, FINDTD, IER)
C
 999  RETURN
C-----------------------------------------------------------------------
 1060 FORMAT ('Task ',A6,' (release of ',A7,') begins')
      END
