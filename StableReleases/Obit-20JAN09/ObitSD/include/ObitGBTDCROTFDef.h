/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;  This program is free software; you can redistribute it and/or    */
/*;  modify it under the terms of the GNU General Public License as   */
/*;  published by the Free Software Foundation; either version 2 of   */
/*;  the License, or (at your option) any later version.              */
/*;                                                                   */
/*;  This program is distributed in the hope that it will be useful,  */
/*;  but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*;  GNU General Public License for more details.                     */
/*;                                                                   */
/*;  You should have received a copy of the GNU General Public        */
/*;  License along with this program; if not, write to the Free       */
/*;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*;  MA 02139, USA.                                                   */
/*;                                                                   */
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitGBTDCROTF structure.    */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitGBTDCROTFListDef.h
 * ObitGBTDCROTFDef structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Output OTF */
ObitOTF *outOTF;
/** Target name */
gchar Name[48];   
/** Correction in sec to add to GBT time labels  */
odouble timeCorr; 
/**number of antenna time samples */
olong nAntTime;   
/** Array of antenna times */
odouble *AntDMJD; 
/** Array of Antenna RA J2000 values */
odouble *AntRA;   
/** Array of Antenna Dec J2000 values */
odouble *AntDec;  
/** Reference Julian date */
odouble refMJD;   
/** Integration time in days */
odouble integTime;
/** Start time of Scan in days */
odouble startTime;
/** End time of Scan in days */
odouble endTime;
/** target number */
ofloat target;
/** scan number */
ofloat scan;
/** First record number (1-rel) in scan */
olong  startRec;
/** End record number (1-rel) in scan */
olong  endRec;
/** Number of feeds  */
olong  nfeed;
/** Number of frequencies */
olong  nchan;
/** Number of Stokes */
olong  nstok;
/** True if data beam switched */
gboolean isBS;
/**  Polarization states per detector, 'R', 'L', 'X', 'Y' */
gchar  poln[10];
/** Channel increment */
ofloat deltaFreq;
/**Frequency reference pixel */
ofloat refPixel;
/** reference frequency (Hz) */
odouble refFrequency;
/** IF information structure*/
ObitGBTIFInfo* IFdata;
/** State information structure */
ObitGBTDCRStateInfo* StateData;
/** Beam offset information structure */
ObitGBTBeamOffInfo* BeamOffData;
/** Az and el offsets to apply to pointing positions */
ofloat azOff, elOff;
