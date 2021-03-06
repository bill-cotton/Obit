/* $Id$   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2012                                              */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitTablePD  structure          */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitTablePDDef.h
 * ObitTablePD structure members for derived classes.
 */
#include "ObitTableDef.h"  /* Parent class definitions */
/** Revision number of the table definition. */
oint  revision;
/** The number of orthogonal polarizations. */
oint  numPol;
/** The number of IFs */
oint  numIF;
/** Number of frequency channels */
oint  numChan;
/** The number of antennas in table. */
oint  numAnt;
/** Polarization parameterazation type, 'ORI-ELP', 'APPROX', 'VLBI' */
gchar  polType[MAXKEYCHARTABLEPD];
/** Column offset for Antenna number in table record */
olong  antNoOff;
/** Physical column number for Antenna number in table record */
olong  antNoCol;
/** Column offset for Subarray number. in table record */
olong  SubAOff;
/** Physical column number for Subarray number. in table record */
olong  SubACol;
/** Column offset for Frequency ID in table record */
olong  FreqIDOff;
/** Physical column number for Frequency ID in table record */
olong  FreqIDCol;
/** Column offset for Reference antenna Poln # NO_POL in table record */
olong  RefAntOff;
/** Physical column number for Reference antenna Poln # NO_POL in table record */
olong  RefAntCol;
/** Column offset for Right-left phase difference per channel/IF, channel varing fastest. in table record */
olong  RLPhaseOff;
/** Physical column number for Right-left phase difference per channel/IF, channel varing fastest. in table record */
olong  RLPhaseCol;
/** Column offset for First parameter, real or ellipticity per channel/IF for poln  # 1, in table record */
olong  Real1Off;
/** Physical column number for First parameter, real or ellipticity per channel/IF for poln  # 1, in table record */
olong  Real1Col;
/** Column offset for Second parameter, imaginary or orientation per channel/IF for poln  # 1, in table record */
olong  Imag1Off;
/** Physical column number for Second parameter, imaginary or orientation per channel/IF for poln  # 1, in table record */
olong  Imag1Col;
/** Column offset for First parameter, real or ellipticity per channel/IF for poln  # 2, in table record */
olong  Real2Off;
/** Physical column number for First parameter, real or ellipticity per channel/IF for poln  # 2, in table record */
olong  Real2Col;
/** Column offset for Second parameter, imaginary or orientation per channel/IF for poln  # 2, in table record */
olong  Imag2Off;
/** Physical column number for Second parameter, imaginary or orientation per channel/IF for poln  # 2, in table record */
olong  Imag2Col;
