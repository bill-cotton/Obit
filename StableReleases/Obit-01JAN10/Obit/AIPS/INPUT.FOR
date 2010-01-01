C   Inputs Class Module
C-----------------------------------------------------------------------
C! Object Oriented AIPS Fortran "INPUT" class library
C# Map-util Utility Object-Oriented
C-----------------------------------------------------------------------
C;  Copyright (C) 1995-1996, 1999
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
C   Class data:
C      Passed Adverbs
C
C   Public functions:
C      AV2INP (prgn, nparm, parm, type, dim, out, ierr)
C         Does AIPS startup and copies AIPS adverbs to an Inputs object
C         restarts AIPS unless DOCRT present and true
C      AV2INT (prgn, nparm, parm, type, dim, out, ierr)
C         For interactive tasks: does AIPS startup and copies AIPS
C         adverbs to an Inputs object, does not restart AIPS
C      IN2OBJ (in, nkey, inkey, outkey, out, ierr)
C         Copies a list of keywords to object out with possible renaming
C      INTELL (nparm, parm, type, dim, out, optell, ierr)
C         Obtain "TELL" instructions as an INPUT object.
C      INGET (name, keywrd, type, dim, value, valuec, ierr)
C         Return adverb value.
C      INPUTT (name, keywrd, type, dim, value, valuec, ierr)
C         Store adverb value.
C      RESAIP (iret, ierr)
C         Resume AIPS now with code return code iret
C   Private functions:
C      INSTRT (prgn, nparm, parm, type, dim, out, ierr)
C         Does AIPS startup and copies AIPS adverbs to an Inputs object
C-----------------------------------------------------------------------
LOCAL INCLUDE 'INPUT.INC'
      INTEGER   MAXLEN
C                                       MAXLEN = max length of strings
      PARAMETER (MAXLEN = 1000)
      CHARACTER STRING*(MAXLEN)
LOCAL END
      SUBROUTINE AV2INP (PRGN, NPARM, PARM, TYPE, DIM, OUT, IERR)
C-----------------------------------------------------------------------
C   Public function
C   Does AIPS startup. Fetches AIPS parameters and saves them in an
C   Input Object.  Also initilizes the object manager and restarts AIPS.
C      Useage note:  If there is a real or integer adverb named DOCRT
C   with a positive value, POPS will not be restarted at the beginning
C   of the task but QUICK will be set false and POPS will be restarted
C   by DIE.
C   Inputs:
C      PRGN    C*6     Task name
C      NPARM   I       Number of AIPS adverbs.
C      PARM    C(?)*8  Names of adverbs.
C      TYPE    I(?)    Types of the adverbs (2=real, 3=character,
C                      4=integer, 5=logical (integer>0 => true.))
C      DIM     I(2,?)  First and second dimensions of adverb arrays.
C      OUT     C*?     Name of the output Input array
C   Outputs:
C      IERR   I  Return error code 0=>OK else failed to create class.
C-----------------------------------------------------------------------
      INTEGER   NPARM
      CHARACTER PRGN*6, PARM(NPARM)*(*), OUT*(*)
      INTEGER   TYPE(NPARM), DIM(2,NPARM), IERR
C
      INTEGER   JERR
      INCLUDE 'INCS:DFIL.INC'
      INCLUDE 'INCS:CLASSIO.INC'
C-----------------------------------------------------------------------
C                                       Do start up stuff
      CALL INSTRT (PRGN, NPARM, PARM, TYPE, DIM, OUT, IERR)
      JERR = 8
      IF (IERR.EQ.0) JERR = 0
C                                       Restart AIPS
      IF (RQUICK) CALL RELPOP (JERR, SBUFF, IERR)
C                                       Previous error condition exists?
      IF (JERR.NE.0) IERR = JERR
C
 999  RETURN
      END
      SUBROUTINE AV2INT (PRGN, NPARM, PARM, TYPE, DIM, OUT, IERR)
C-----------------------------------------------------------------------
C   Public function
C   Does AIPS startup for an interactive task.. Fetches AIPS parameters
C   and saves them in an Input object.  Also initilizes the object
C   manager.  AIPS (POPS) will be restarted by DIE.
C   Inputs:
C      PRGN    C*6     Task name
C      NPARM   I       Number of AIPS adverbs.
C      PARM    C(?)*8  Names of adverbs.
C      TYPE    I(?)    Types of the adverbs (2=real, 3=character,
C                      4=integer, 5=logical (integer>0 => true.))
C      DIM     I(2,?)  First and second dimensions of adverb arrays.
C      OUT     C*?     Name of the output Input array
C   Outputs:
C      IERR   I  Return error code 0=>OK else failed to create class.
C-----------------------------------------------------------------------
      INTEGER   NPARM
      CHARACTER PRGN*6, PARM(NPARM)*(*), OUT*(*)
      INTEGER   TYPE(NPARM), DIM(2,NPARM), IERR
C
      INCLUDE 'INCS:DFIL.INC'
C-----------------------------------------------------------------------
C                                       Do start up stuff
      CALL INSTRT (PRGN, NPARM, PARM, TYPE, DIM, OUT, IERR)
C
 999  RETURN
      END
      SUBROUTINE IN2OBJ (IN, NKEY, INKEY, OUTKEY, OUT, IERR)
C-----------------------------------------------------------------------
C   Copies a list keyword/values from an Inputs object to another object
C   with possible renaming of the keyword.
C   Inputs:
C      IN     C*?     Inputs object (should also work for any object)
C      NKEY   I       Number of keywords
C      INKEY  C(?)*8  Keywords to copy
C      OUTKEY C(?)*?  Keywords to copy, if blank use INKEY
C      OUT    C*?     Output object
C   Outputs:
C      IERR   I       Return code, 0=OK else failed.
C-----------------------------------------------------------------------
      INTEGER   NKEY, IERR
      CHARACTER IN*(*), INKEY(*)*(*), OUTKEY(*)*(*), OUT*(*)
C
      INTEGER   LOOP, TYPE, DIM(7), BUFNO, OBJIN, OBJOUT, CLASNO, JERR
      CHARACTER NEWKEY*32, CNAME*8
      INCLUDE 'INCS:OBJPARM.INC'
      INCLUDE 'INCS:CLASSIO.INC'
      INCLUDE 'INPUT.INC'
      INCLUDE 'INCS:DMSG.INC'
C-----------------------------------------------------------------------
C                                       Object numbers
      CALL OBNAME (IN, OBJIN, IERR)
      IF (IERR.NE.0) GO TO 990
      CALL OBNAME (OUT, OBJOUT, IERR)
      IF (IERR.NE.0) GO TO 990
C                                       Output Class
      CALL OBCLAS (OBJOUT, CLASNO, CNAME, IERR)
      IF (IERR.NE.0) GO TO 990
C                                       Open output object, use buffer
C                                       for temporary storage.
      CALL OBOPEN (OUT, IERR)
      IF (IERR.NE.0) GO TO 990
      CALL OBINFO (OUT, BUFNO, IERR)
      IF (IERR.NE.0) GO TO 980
C                                       Copy
      DO 100 LOOP = 1,NKEY
         NEWKEY = OUTKEY(LOOP)
         IF (NEWKEY.EQ.' ') NEWKEY = INKEY(LOOP)
         CALL OBGET (OBJIN, INKEY(LOOP), TYPE, DIM, OBUFFR(1,BUFNO),
     *      STRING, IERR)
         IF (IERR.NE.0) GO TO 980
         CALL OPUT (OUT, NEWKEY, TYPE, DIM, OBUFFR(1,BUFNO),
     *      STRING, IERR)
         IF (IERR.NE.0) GO TO 980
 100     CONTINUE
C                                       Close I/O
 980  CALL OBCLOS (OUT, JERR)
      IF (IERR.EQ.0) IERR = JERR
C                                       Error
 990  IF (IERR.NE.0) THEN
         MSGTXT = 'IN2OBJ: ERROR COPYING FROM:' // IN
         CALL MSGWRT (7)
         MSGTXT = 'IN2OBJ:                 TO:' // OUT
         CALL MSGWRT (7)
         END IF
C
 999  RETURN
      END
      SUBROUTINE INTELL (NPARM, PARM, TYPE, DIM, OUT, OPTELL,
     *   IERR)
C-----------------------------------------------------------------------
C   Public function
C   Fetches AIPS parameters from a TELL operations and saves them in an
C   Input Object.  Also return, the TELL operation type.
C   Inputs:
C      NPARM   I       Number of AIPS adverbs.
C      PARM    C(?)*8  Names of adverbs.
C      TYPE    I(?)    Types of the adverbs (2=real, 3=character,
C                      4=integer, 5=logical (integer>0 => true.))
C      DIM     I(2,?)  First and second dimensions of adverb arrays.
C      OUT     C*?     Name of the output INPUT array
C      OPTELL  C*4     Operation type:
C   Outputs:
C                      '    ' => no new instruction.
C                      'ABOR' => declare a failure and quit
C                      'QUIT' => done enough, shut done gracefully.
C                      'CHAN' => new parameters
C      IERR    I       Return error code 0=>OK else failed.
C-----------------------------------------------------------------------
      INTEGER   NPARM
      CHARACTER PARM(NPARM)*(*), OUT*(*), OPTELL*4
      INTEGER   TYPE(NPARM), DIM(2,NPARM), IERR
C
      INTEGER   LOOP, NWORDS, NUM, BUFNO, IROUND, LOOP2, NCHAR, NSTRIN,
     *   POINT, PNTW, PNTC, OBJNUM, DUMMY, JERR
      CHARACTER CDUMMY*1
      INCLUDE 'INCS:DMSG.INC'
      INCLUDE 'INCS:CLASSIO.INC'
      INCLUDE 'INCS:PAOOF.INC'
      INCLUDE 'INPUT.INC'
      INTEGER   IBUFFR(BUFSIZ,MAXIO)
      LOGICAL   LBUFFR(BUFSIZ,MAXIO), EXIST
      EQUIVALENCE (IBUFFR, LBUFFR, OBUFFR)
C-----------------------------------------------------------------------
      IERR = 0
      OPTELL = '    '
C                                       Create output if necessary
      CALL OBFEXS (OUT, EXIST, IERR)
      IF (IERR.NE.0) GO TO 900
      IF (.NOT.EXIST) THEN
         CALL OBCREA (OUT, 'INPUTS', IERR)
         IF (IERR.NE.0) GO TO 900
         END IF
      CALL OBNAME (OUT, OBJNUM, IERR)
      IF (IERR.NE.0) GO TO 900
C                                       Find number of AIPS parameters
      NWORDS = 0
      DO 100 LOOP = 1,NPARM
         IF (TYPE(LOOP).EQ.OOARE) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
         ELSE IF (TYPE(LOOP).EQ.OOACAR) THEN
            NUM = MAX (1, DIM(1,LOOP))
            NUM = ((NUM - 1) / 4 + 1) * MAX (1, DIM(2,LOOP))
         ELSE IF (TYPE(LOOP).EQ.OOAINT) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
         ELSE IF (TYPE(LOOP).EQ.OOALOG) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
         ELSE
C                                       Illegal type
            IERR = 1
            WRITE (MSGTXT,1000) TYPE(LOOP), PARM(LOOP)
            CALL MSGWRT (7)
            GO TO 900
            END IF
         NWORDS = NWORDS + NUM
 100     CONTINUE
C                                       Open object to get buffer for
C                                       GTTELL.
      CALL OBOPEN (OUT, IERR)
      IF (IERR.NE.0) GO TO 900
      CALL OBINFO (OUT, BUFNO, IERR)
      IF (IERR.NE.0) GO TO 900
C                                       Read adverbs
      CALL GTTELL (NWORDS, OPTELL, OBUFFR(1,BUFNO), SBUFF, IERR)
      IF (IERR.EQ.0) OPTELL = '    '
      IF (IERR.LT.3) IERR = 0
      IF (OPTELL.NE.'CHAN') GO TO 980
C                                       Copy Adverbs to inputs
      POINT = 1
      DO 200 LOOP = 1,NPARM
         NUM = 0
C                                       Reals
         IF (TYPE(LOOP).EQ.OOARE) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
            CALL OBPUT (OBJNUM, PARM(LOOP), OOARE, DIM(1,LOOP),
     *         OBUFFR(POINT,BUFNO), CDUMMY, IERR)
            IF (IERR.NE.0) GO TO 900
C                                       Characters
         ELSE IF (TYPE(LOOP).EQ.OOACAR) THEN
            NCHAR = MAX (1, DIM(1,LOOP))
            NSTRIN = MAX (1, DIM(2,LOOP))
            NUM = ((NCHAR - 1) / 4 + 1) * MAX (1, DIM(2,LOOP))
C                                       Check buffer size for strings.
            IF ((NCHAR*NSTRIN).GT.MAXLEN) THEN
               IERR = 1
               WRITE (MSGTXT,1100) PARM(LOOP)
               CALL MSGWRT (7)
               GO TO 900
               END IF
            PNTW = POINT
            PNTC = 1
            DO 140 LOOP2 = 1,NSTRIN
               CALL H2CHR (NCHAR, 1, OBUFFR(PNTW,BUFNO),
     *            STRING(PNTC:PNTC+NCHAR-1))
               PNTW = PNTW + ((NCHAR - 1) / 4 + 1)
               PNTC = PNTC + NCHAR
 140           CONTINUE
            CALL OBPUT (OBJNUM, PARM(LOOP), OOACAR, DIM(1,LOOP), DUMMY,
     *         STRING, IERR)
            IF (IERR.NE.0) GO TO 900
C                                       Integers
         ELSE IF (TYPE(LOOP).EQ.OOAINT) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
            DO 150 LOOP2 = 1,NUM
               IBUFFR(POINT+LOOP2-1,BUFNO) =
     *            IROUND (OBUFFR(POINT+LOOP2-1,BUFNO))
 150           CONTINUE
            CALL OBPUT (OBJNUM, PARM(LOOP), OOAINT, DIM(1,LOOP),
     *         IBUFFR(POINT,BUFNO), CDUMMY, IERR)
            IF (IERR.NE.0) GO TO 900
C                                       Logical
         ELSE IF (TYPE(LOOP).EQ.OOALOG) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
            DO 160 LOOP2 = 1,NUM
               LBUFFR(POINT+LOOP2-1,BUFNO) =
     *            OBUFFR(POINT+LOOP2-1,BUFNO) .GT. 0.0
 160           CONTINUE
            CALL OBPUT (OBJNUM, PARM(LOOP), OOALOG, DIM(1,LOOP),
     *         LBUFFR(POINT,BUFNO), CDUMMY, IERR)
            IF (IERR.NE.0) GO TO 900
            END IF
         POINT = POINT + NUM
 200     CONTINUE
C                                       Error?
 900  IF (IERR.NE.0) THEN
         MSGTXT = 'INTELL: ERROR DURING TASK INITIALIZATION'
         CALL MSGWRT (7)
         END IF
C                                       Close I/O
 980  CALL OBCLOS (OUT, JERR)
      IF (IERR.EQ.0) IERR = JERR
C
 999  RETURN
C-----------------------------------------------------------------------
 1000 FORMAT ('INTELL: ILLEGAL TYPE CODE ',I3,' FOR ',A)
 1100 FORMAT ('INTELL: BLEW INTERNAL BUFFER FOR STRINGS ON ', A)
      END
      SUBROUTINE INGET (NAME, KEYWRD, TYPE, DIM, VALUE, VALUEC, IERR)
C-----------------------------------------------------------------------
C   Public function
C   Returns the dimensionality and value(s) of the value associated with
C   a given keyword.
C   Inputs:
C      NAME     C*?   The name of the object.
C      KEYWRD   C*?   The name of the keyword.
C   Outputs:
C      TYPE     I     Data type: 1=D, 2=R, 3=C, 4=I, 5=L
C      DIM      I(*)  Dimensionality of value, an axis dimension of zero
C                     means that that dimension and higher are
C                     undefined.
C      VALUE    R(*)  The value associated with keyword.
C      VALUEC   C*?   Associated value (character)
C      IERR     I     Error return code, 0=OK
C-----------------------------------------------------------------------
      CHARACTER NAME*(*), KEYWRD*(*), VALUEC*(*)
      INTEGER   TYPE, DIM(*), VALUE(*), IERR
C
      INTEGER   OBJNUM
      INCLUDE 'INCS:DMSG.INC'
C-----------------------------------------------------------------------
      IERR = 0
C                                       Lookup NAME
      CALL OBNAME (NAME, OBJNUM, IERR)
      IF (IERR.NE.0) GO TO 990
      CALL OBGET (OBJNUM, KEYWRD, TYPE, DIM, VALUE, VALUEC, IERR)
      IF (IERR.NE.0) GO TO 990
      GO TO 999
C                                       Error
 990  MSGTXT = 'INGET: ERROR ACCESSING:' // NAME
      CALL MSGWRT (7)
C
 999  RETURN
      END
      SUBROUTINE INPUTT (NAME, KEYWRD, TYPE, DIM, VALUE, VALUEC, IERR)
C-----------------------------------------------------------------------
C   Public function
C   Stores the integer (array) associated with a given keyword.
C   A new keyword will be created if necessary.
C   Inputs:
C      NAME     C*?   The name of the object.
C      KEYWRD   C*?   The name of the keyword.
C      TYPE     I     Data type: 1=D, 2=R, 3=C, 4=I, 5=L
C      DIM      I(*)  Dimensionality of value, an axis dimension of zero
C                     means that that dimension and higher are
C                     undefined.
C      VALUE    I(*)  The value associated with keyword.
C      VALUEC   C*?   Associated value (character)
C   Outputs:
C      IERR     I     Error return code, 0=OK
C-----------------------------------------------------------------------
      CHARACTER NAME*(*), KEYWRD*(*), VALUEC*(*)
      INTEGER   TYPE, DIM(*), VALUE(*), IERR
C
      INTEGER   OBJNUM
      INCLUDE 'INCS:DMSG.INC'
C-----------------------------------------------------------------------
      IERR = 0
C                                       Lookup NAME
      CALL OBNAME (NAME, OBJNUM, IERR)
      IF (IERR.NE.0) GO TO 990
      CALL OBPUT (OBJNUM, KEYWRD, TYPE, DIM, VALUE, VALUEC, IERR)
      IF (IERR.NE.0) GO TO 990
      GO TO 999
C                                       Error
 990  MSGTXT = 'INPUTT: ERROR ACCESSING:' // NAME
      CALL MSGWRT (7)
C
 999  RETURN
      END
      SUBROUTINE RESAIP (IRET, IERR)
C-----------------------------------------------------------------------
C   Public function
C   Resumes AIPS (POPS) sometime after the initial start of the task.
C   Usage note:  For tasks started by AV2INT, this routine allows you to
C   resume AIPS after all interactive functions are over,but before the
C   end of the task.
C   Inputs:
C      IRET   I   Return code to pass to AIPS
C   Outputs:
C      IERR   I  Return error code 0=>OK else failed to resume
C-----------------------------------------------------------------------
      INTEGER   IRET, IERR
C
      INCLUDE 'INCS:DMSG.INC'
      INCLUDE 'INCS:DDCH.INC'
      INCLUDE 'INCS:DFIL.INC'
      INCLUDE 'INCS:CLASSIO.INC'
C-----------------------------------------------------------------------
C                                       Only if not already resumed
C                                       and not batch: Restart AIPS
      IF ((.NOT.RQUICK) .AND. (NPOPS.LE.NINTRN) .AND. (ISBTCH.NE.32000))
     *   THEN
         CALL RELPOP (IRET, SBUFF, IERR)
         RQUICK = .TRUE.
         END IF
C
 999  RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Private functions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INSTRT (PRGN, NPARM, PARM, TYPE, DIM, OUT, IERR)
C-----------------------------------------------------------------------
C   Private function
C   Does AIPS startup. Fetches AIPS parameters and saves them in an
C   Input Object.  Also initilizes the object manager and restarts AIPS.
C      Useage note:  If there is a real or integer adverb named DOCRT
C   with a positive value, POPS will not be restarted at the beginning
C   of the task but QUICK will be set false and POPS will be restarted
C   by DIE.
C   Inputs:
C      PRGN    C*6     Task name
C      NPARM   I       Number of AIPS adverbs.
C      PARM    C(?)*8  Names of adverbs.
C      TYPE    I(?)    Types of the adverbs (2=real, 3=character,
C                      4=integer, 5=logical (integer>0 => true.))
C      DIM     I(2,?)  First and second dimensions of adverb arrays.
C      OUT     C*?     Name of the output Input array
C   Outputs:
C      IERR    I       Return error code 0=>OK else failed to create
C                      class.
C-----------------------------------------------------------------------
      INTEGER   NPARM
      CHARACTER PRGN*6, PARM(NPARM)*(*), OUT*(*)
      INTEGER   TYPE(NPARM), DIM(2,NPARM), IERR
C
      INTEGER   LOOP, NWORDS, NUM, BUFNO, IROUND, LOOP2, NCHAR, NSTRIN,
     *   POINT, PNTW, PNTC, JERR, OBJNUM, TY, DM(3), IDOCRT, MSGSAV,
     *   DUMMY
      REAL      RDOCRT
      LOGICAL   INTER
      CHARACTER CDUMMY*1
      INCLUDE 'INCS:DMSG.INC'
      INCLUDE 'INCS:DFIL.INC'
      INCLUDE 'INCS:DDCH.INC'
      INCLUDE 'INCS:DHDR.INC'
      INCLUDE 'INCS:OBJPARM.INC'
      INCLUDE 'INCS:CLASSIO.INC'
      INCLUDE 'INCS:PAOOF.INC'
      INCLUDE 'INPUT.INC'
      INTEGER   IBUFFR(BUFSIZ,MAXIO)
      LOGICAL   LBUFFR(BUFSIZ,MAXIO)
      EQUIVALENCE (IBUFFR, LBUFFR, OBUFFR), (IDOCRT, RDOCRT)
C-----------------------------------------------------------------------
      IERR = 0
C                                       Init for AIPS, disks, ...
      CALL ZDCHIN (.TRUE.)
      CALL VHDRIN
C                                       Initialize /CFILES/
      NSCR = 0
      NCFILE = 0
C                                       In case of error before RELPOP
C                                       call.
      JERR = 8
C                                       Initialize Object manager
      CALL OBINIT (IERR)
      IF (IERR.NE.0) GO TO 900
C                                       Create Input object
      CALL OBCREA (OUT, 'INPUTS', IERR)
      IF (IERR.NE.0) GO TO 900
      CALL OBNAME (OUT, OBJNUM, IERR)
      IF (IERR.NE.0) GO TO 900
C                                       Find number of AIPS parameters
      NWORDS = 0
      DO 100 LOOP = 1,NPARM
         IF (TYPE(LOOP).EQ.OOARE) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
         ELSE IF (TYPE(LOOP).EQ.OOACAR) THEN
            NUM = MAX (1, DIM(1,LOOP))
            NUM = ((NUM - 1) / 4 + 1) * MAX (1, DIM(2,LOOP))
         ELSE IF (TYPE(LOOP).EQ.OOAINT) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
         ELSE IF (TYPE(LOOP).EQ.OOALOG) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
         ELSE
C                                       Illegal type
            IERR = 1
            WRITE (MSGTXT,1000) TYPE(LOOP), PARM(LOOP)
            CALL MSGWRT (7)
            GO TO 900
            END IF
         NWORDS = NWORDS + NUM
 100     CONTINUE
C                                       Open object to get buffer for
C                                       GTPARM.
      CALL OBOPEN (OUT, IERR)
      IF (IERR.NE.0) GO TO 900
      CALL OBINFO (OUT, BUFNO, IERR)
      IF (IERR.NE.0) GO TO 900
C                                       Read adverbs
      CALL GTPARM (PRGN, NWORDS, RQUICK, OBUFFR(1,BUFNO), SBUFF, IERR)
      IF (IERR.NE.0) THEN
         IF (IERR.EQ.1) GO TO 900
         RQUICK = .FALSE.
         WRITE (MSGTXT,1001) IERR
         CALL MSGWRT (8)
         END IF
C                                       Copy Adverbs to inputs
      POINT = 1
      DO 200 LOOP = 1,NPARM
         NUM = 0
C                                       Reals
         IF (TYPE(LOOP).EQ.OOARE) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
            CALL OBPUT (OBJNUM, PARM(LOOP), OOARE, DIM(1,LOOP),
     *         OBUFFR(POINT,BUFNO), CDUMMY, IERR)
            IF (IERR.NE.0) GO TO 900
C                                       Characters
         ELSE IF (TYPE(LOOP).EQ.OOACAR) THEN
            NCHAR = MAX (1, DIM(1,LOOP))
            NSTRIN = MAX (1, DIM(2,LOOP))
            NUM = ((NCHAR - 1) / 4 + 1) * MAX (1, DIM(2,LOOP))
C                                       Check buffer size for strings.
            IF ((NCHAR*NSTRIN).GT.MAXLEN) THEN
               IERR = 1
               WRITE (MSGTXT,1100) PARM(LOOP)
               CALL MSGWRT (7)
               GO TO 900
               END IF
            PNTW = POINT
            PNTC = 1
            DO 140 LOOP2 = 1,NSTRIN
               CALL H2CHR (NCHAR, 1, OBUFFR(PNTW,BUFNO),
     *            STRING(PNTC:PNTC+NCHAR-1))
               PNTW = PNTW + ((NCHAR - 1) / 4 + 1)
               PNTC = PNTC + NCHAR
 140           CONTINUE
            CALL OBPUT (OBJNUM, PARM(LOOP), OOACAR, DIM(1,LOOP), DUMMY,
     *         STRING, IERR)
            IF (IERR.NE.0) GO TO 900
C                                       Integers
         ELSE IF (TYPE(LOOP).EQ.OOAINT) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
            DO 150 LOOP2 = 1,NUM
               IBUFFR(POINT+LOOP2-1,BUFNO) =
     *            IROUND (OBUFFR(POINT+LOOP2-1,BUFNO))
 150           CONTINUE
            CALL OBPUT (OBJNUM, PARM(LOOP), OOAINT, DIM(1,LOOP),
     *         IBUFFR(POINT,BUFNO), CDUMMY, IERR)
            IF (IERR.NE.0) GO TO 900
C                                       Logical
         ELSE IF (TYPE(LOOP).EQ.OOALOG) THEN
            NUM = MAX (1, DIM(1,LOOP)) * MAX (1,DIM(2,LOOP))
            DO 160 LOOP2 = 1,NUM
               LBUFFR(POINT+LOOP2-1,BUFNO) =
     *            OBUFFR(POINT+LOOP2-1,BUFNO) .GT. 0.0
 160           CONTINUE
            CALL OBPUT (OBJNUM, PARM(LOOP), OOALOG, DIM(1,LOOP),
     *         LBUFFR(POINT,BUFNO), CDUMMY, IERR)
            IF (IERR.NE.0) GO TO 900
            END IF
         POINT = POINT + NUM
 200     CONTINUE
C                                       Set RQUICK by value of DOCRT if
C                                       present.

      MSGSAV = MSGSUP
      MSGSUP = 32000
      CALL OBGET (OBJNUM, 'DOCRT', TY, DM, IDOCRT, CDUMMY, IERR)
      MSGSUP = MSGSAV
      IF (IERR.EQ.0) THEN
         INTER = ((TY.EQ.4).AND.(IDOCRT.GE.1)) .OR.
     *      ((TY.EQ.2).AND.(RDOCRT.GE.0.01))
         RQUICK = RQUICK .AND. (.NOT.INTER)
         END IF
      IERR = 0
      JERR = 0
C                                       Error?
 900  IF (IERR.NE.0) THEN
         MSGTXT = 'INSTRT: ERROR DURING TASK INITIALIZATION'
         CALL MSGWRT (7)
         END IF
C                                       Previous error condition exists?
      IF (JERR.NE.0) THEN
         IERR = JERR
         GO TO 999
         END IF
C                                       Close I/O
      CALL OBCLOS (OUT, IERR)
      IF (IERR.NE.0) GO TO 999
C
 999  RETURN
C-----------------------------------------------------------------------
 1000 FORMAT ('INSTRT: ILLEGAL TYPE CODE ',I3,' FOR ',A)
 1001 FORMAT ('INSTRT: ERROR',I3,' OBTAINING INPUT PARAMETERS')
 1100 FORMAT ('INSTRT: BLEW INTERNAL BUFFER FOR STRINGS ON ', A)
      END
