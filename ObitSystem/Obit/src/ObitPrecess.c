/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2010                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include <sys/types.h>
#include <time.h>
#include <math.h>
#include "ObitPrecess.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPrecess.c
 * ObitPrecess module function definitions.
 *
 * Functions to presess celestial coordinates.
 */

/*---------------Private function prototypes----------------*/
/** Private: Precess between apparent and J2000 epoch positions */
static void jpreces(double JD, ofloat equin, double deldat, olong dir,
		    gboolean gr, odouble obspos[3], ofloat polar[2],
		    odouble *RAMean, odouble *DecMean, odouble *RAApp, odouble *DecApp);

/** Private: Compute the rotation matrix for precession */
static void jprenu(olong dir, double JD, ofloat equinox, gboolean doNut, 
		   odouble rotMat[3][3]);

/** Private: Compute aberation and GR light bending */
static void jaber(odouble JD, ofloat equin, gboolean diurn, 
		  odouble rhoGeo, odouble phiGeo, odouble rLST,
		  odouble poso[3], odouble velo[3]);

/** Private:  Correct rectangular position for polar motion */
static void jpolar(olong dir, ofloat polarx, ofloat polary, odouble pos[3]);

/** Private:  Computes nutation from IAU 1980 series */
static void jnut(odouble JD, odouble *delPsi, odouble *delEps);

/** Private:  Earth orbit ephemeris */
static void EarthEphem(odouble JD, ofloat equin, odouble dvb[3], odouble dpb[3], 
		       odouble dvh[3], odouble dph[3]);
/*----------------------Public functions---------------------------*/

/**
 * Precess a source with a UVDesc from the position at a standard epoch
 * to the apparent position.
 * \param desc      UV data descriptor
 * \param source    Obit source with position to precess.
 */
void ObitPrecessUVJPrecessApp (ObitUVDesc *desc, ObitSource *source)
{
  odouble JD, RAMean, DecMean;
  odouble obsPos[] = {0.0, 0.0, 0.0};
  ofloat polar[]   = {0.0, 0.0};
  
  /* error checks */
  g_assert (ObitUVDescIsA(desc));
  g_assert (ObitSourceIsA(source));
  
  /* Get Julian Date of observations */
  ObitUVDescDate2JD (desc->obsdat, &JD);
  
  /* Precess */
  RAMean  = DG2RAD*source->RAMean;
  DecMean = DG2RAD*source->DecMean;
  jpreces (JD, desc->equinox, 0.01, 1, FALSE, obsPos, polar,
	   &RAMean, &DecMean, &source->RAApp, &source->DecApp);

  /* Convert to degrees */
  source->RAApp  *= RAD2DG;
  source->DecApp *= RAD2DG;
  
} /* end ObitPrecessUVJPrecessApp */

/**
 * Precess the RA and declination from a UV data header to the 
 * apparent position of date
 * \param desc      UV data descriptor
 * \param RAApp  [out] RA of date [deg]
 * \param DecApp [out] Declination of date [deg]
 */
void ObitPrecessUVRaDecApp (ObitUVDesc *desc, odouble *RAApp, odouble *DecApp)
{
  odouble JD, RAMean, DecMean;
  odouble obsPos[] = {0.0, 0.0, 0.0};
  ofloat polar[]   = {0.0, 0.0};
  
  /* error checks */
  g_assert (ObitUVDescIsA(desc));
  
  /* Get Julian Date of observations */
  ObitUVDescDate2JD (desc->obsdat, &JD);
  
  /* Precess */
  if (desc->crval[desc->jlocr]<0.0) desc->crval[desc->jlocr] += 360.0; /* Patch AIPS++ corruption */
  RAMean  = desc->crval[desc->jlocr]*DG2RAD;
  DecMean = desc->crval[desc->jlocd]*DG2RAD;
  jpreces (JD, desc->equinox, 0.01, 1, FALSE, obsPos, polar,
	   &RAMean, &DecMean, RAApp, DecApp);

  /* Convert to degrees */
  *RAApp  *= RAD2DG;
  *DecApp *= RAD2DG;
  
} /* end ObitPrecessUVRaDecApp */

/**
 * Predict the Greenwich Sidereal Time at UT=0 and the Earth's
 * rotation rate on a given Julian date.
 * \param JD
 * \param GSTUTC0  [out] Apparent GST (hours) at UTC=0 on JD
 * \param Rate     [out] Earth rotation rate turns per day on JD
 */
void ObitPrecessGST0 (odouble JD, odouble *GSTUTC0, odouble *Rate)
{
  odouble UTC, TC, Eps, EqEq , DelPsi, DelEps, JD0, TU, GMSTM;

  UTC = JD - 2400000.5;
  /* Get equation of equinoxes. */
  TC = (JD - 2433282.423) / 36524.21988;
  Eps = -46.850 * TC - 0.0034 * TC * TC + 0.0018 * TC * TC * TC;
  Eps = (84404.84 + Eps) / 3600.0 / 180.0*3.141592658979;
  /* Nutation terms */
  jnut (JD,  &DelPsi, &DelEps);
  EqEq = DelPsi * cos (Eps);

  /* mean GST at midnight */
  JD0 = 2415020.0;
  TU = (JD - JD0) / 36525.0;
  GMSTM = (8640184.542/3600.0) * TU;
  GMSTM = fmod (GMSTM, 24.0);
  GMSTM = (6.0 + 38.0/60.0 + 45.836/3600.0) + GMSTM + (0.0929/3600.0) * TU * TU;

  *GSTUTC0 = GMSTM + (EqEq / (DG2RAD*15.0));
  *Rate = 1.00273790265 + 0.589e-10 * TU;
}  /* end ObitPrecessGST0 */

/*----------------------Private functions---------------------------*/

/**
 * Routine to precess positions using the Julian IAU 1984 conventions,
 * (i.e. J2000 positions).  Optional corrections can be made for
 * relativistic bending of the light by the sun, diurnal aberation and
 * polar motion.  Proper motion and parallax are assumed negligible.
 * 
 * Adapted from the AIPSish JPRECES.FOR
 * \param JD      Julian date of observation (e.g. 2446754.123445)
 * \param equin   Epoch of mean equinox (e.g. 2000.0)
 * \param deldat  Interpolation interval; compute precession etc.
 *                parameters at this interval and do a linear
 *                interpolation (days).
 * \param dir     1 => convert from mean to apparent;
 *               -1 => convert from apparent to mean.
 * \param gr      If true correct apparent position for the general
 *                relativistic bending of light by the sun.
 * \param obspos  Earth centered location of the observations.
 *                If non zero then diurnal aberation corrections are
 *                made for this location.
 *                1 = Latitude (radians)
 *                2 = East longitude
 *                3 = radius from earth center (meter).
 * \param polar  X and Y position of the pole (arcsec.)
 *                If non zero then the apparent position is corrected
 *                for the position of the pole.
 *                Note: this correction is not desirable if the antenna
 *                positions are corrected for polar motion. 
 * \param RAMean  Right ascension at the mean epoch (radians).
 * \param DecMean Declination at mean epoch.
 * \param RAApp   Apparent Right ascension at JD and OBSPOS
 * \param DecApp  Apparent declination.
 */
static void jpreces(double JD, ofloat equin, double deldat, olong dir,
		    gboolean gr, odouble obspos[3], ofloat polar[2],
		    odouble *RAMean, odouble *DecMean, odouble *RAApp, odouble *DecApp)
{
  olong itemp;
  gboolean diurn;
  odouble prnmat[3][3], velo[3], rlst, pos[3], tu, gmst, 
    twopi, time, e[3], ibige, pdote, konst, v[3], beta, pdotv, konst2, rhogeo;
  odouble poso[3]={0.0,0.0,0.0}, out[3]={0.0,0.0,0.0};

  twopi = 2.0 * G_PI;

  /* convert to rectangular coordinates. */
  if (dir > 0) {
    pos[0] = cos (*RAMean) * cos (*DecMean);
    pos[1] = sin (*RAMean) * cos (*DecMean);
    pos[2] = sin (*DecMean);
  } else {
    pos[0] = cos (*RAApp) * cos (*DecApp);
    pos[1] = sin (*RAApp) * cos (*DecApp);
    pos[2] = sin (*DecApp);
  } 

  /* Get Precession matrix */
  jprenu (dir, JD, equin, TRUE, prnmat);

  /* aberation and light bending */
  diurn = (abs (obspos[2])  >  1.0e-5);
  rhogeo = 0.0;
  rlst = 0.0;

  /* following for diurnal aberation. */
  if (diurn) {
    rhogeo = obspos[2] / 6378140.0;
    itemp = JD;
    time = (JD - itemp - 0.5) * twopi * 1.002737778;
    tu = (JD - 2451545.0) / 36525.0;
    gmst = time + 
	((((((-6.2e-6 * tu) + 0.093104) * tu) + 8640184.812866) * tu + 24110.54841) * twopi / 86400.0);
    gmst  = fmod (gmst, twopi);
    rlst = obspos[1] + gmst;
    } 

  /* Get aberation terms */
  jaber (JD, equin, diurn, rhogeo, obspos[0], rlst, poso, velo);

  /* reduce position. */
  if (dir > 0) {
    /* Mean to Apparent. */
    /* Light deflection?  From Astr. Alm. 1986 b40. */
    if (gr) {
      ibige = 1.0 / sqrt (poso[0]*poso[0] + poso[1]*poso[1] + poso[2]*poso[2]);
      e[0] = poso[0] * ibige;
      e[1] = poso[1] * ibige;
      e[2] = poso[2] * ibige;
      pdote = pos[0]* e[0] + pos[1]*e[1] + pos[2]*e[2];
      konst = (1.974e-8 * ibige) / (1.0 + pdote);
      pos[0] = pos[0] + konst * (e[0] - pdote*pos[0]);
      pos[1] = pos[1] + konst * (e[1] - pdote*pos[1]);
      pos[2] = pos[2] + konst * (e[2] - pdote*pos[2]);
    } 

    /* aberation */
    v[0] = velo[0] * 0.0057755;
    v[1] = velo[1] * 0.0057755;
    v[2] = velo[2] * 0.0057755;
    beta = sqrt (1.0 - (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
    pdotv = pos[0]* v[0] + pos[1]*v[1] + pos[2]*v[2];
    konst = (1.0 + (pdotv / (1.0 + beta))) / (1.0 + pdotv);
    konst2 = beta / (1.0 + pdotv);
    pos[0] = pos[0]*konst2 + konst*v[0];
    pos[1] = pos[1]*konst2 + konst*v[1];
    pos[2] = pos[2]*konst2 + konst*v[2];

    /* precession and nutation */
    out[0] = pos[0]*prnmat[0][0] + pos[1]*prnmat[0][1] + pos[2]*prnmat[0][2];
    out[1] = pos[0]*prnmat[1][0] + pos[1]*prnmat[1][1] + pos[2]*prnmat[1][2];
    out[2] = pos[0]*prnmat[2][0] + pos[1]*prnmat[2][1] + pos[2]*prnmat[2][2];

  } else {
    /* Apparent to Mean */
    /* precession and nutation */
    out[0] = pos[0]*prnmat[0][0] + pos[1]*prnmat[0][1] + pos[2]*prnmat[0][2];
    out[1] = pos[0]*prnmat[1][0] + pos[1]*prnmat[1][1] + pos[2]*prnmat[1][2];
    out[2] = pos[0]*prnmat[2][0] + pos[1]*prnmat[2][1] + pos[2]*prnmat[2][2];

    /* aberation */
    v[0] = velo[0] * 0.0057755;
    v[1] = velo[1] * 0.0057755;
    v[2] = velo[2] * 0.0057755;
    beta = sqrt (1.0 - (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
    pdotv = out[0]* v[0] + out[1]*v[1] + out[2]*v[2];
    konst = (1.0 + (pdotv / (1.0 + beta))) / (1.0 + pdotv);
    konst2 = (1.0 + pdotv) / beta;
    out[0] = (out[0] - konst*v[0]) * konst2;
    out[1] = (out[1] - konst*v[1]) * konst2;
    out[2] = (out[2] - konst*v[2]) * konst2;

    /* light deflection?  From Astr. Alm. 1986 b40. */
    if (gr) {
      ibige = 1.0 / sqrt (poso[0]*poso[0] + poso[1]*poso[1] + poso[2]*poso[2]);
      e[0] = poso[0] * ibige;
      e[1] = poso[1] * ibige;
      e[2] = poso[2] * ibige;
      pdote = out[0]* e[0] + out[1]*e[1] + out[2]*e[2];
      konst = (1.974e-8 * ibige) / (1.0 + pdote);
      pos[0] = out[0] - konst * (e[0] - pdote*out[0]);
      pos[1] = out[1] - konst * (e[1] - pdote*out[1]);
      pos[2] = out[2] - konst * (e[2] - pdote*out[2]);
      /* iterate once */
      pdote = pos[0]* e[0] + pos[1]*e[1] + pos[2]*e[2];
      konst = (1.974e-8 * ibige) / (1.0 + pdote);
      out[0] = pos[0] - konst * (e[0] - pdote*pos[0]);
      out[1] = pos[1] - konst * (e[1] - pdote*pos[1]);
      out[2] = pos[2] - konst * (e[2] - pdote*pos[2]);
    } 
  } 

  /* polar motion */
  if ((fabs (polar[0]) + fabs (polar[1])) > 1.0e-20)
    jpolar (dir, polar[0], polar[1], out);

  /* convert to spherical coord. */
  if (dir > 0) {
    /* mean to apparent */
    *RAApp = atan2 (out[1], out[0]);
    if (*RAApp < 0.0) *RAApp = *RAApp + twopi;
    *DecApp = asin (out[2]);

  } else {
    /* apparent to mean */
    *RAMean = atan2 (out[1], out[0]);
    if (*RAMean < 0.0) *RAMean =* RAMean + twopi;
    *DecMean = asin (out[2]);
  } 
} /* end jpreces */

/**
 * Routine to compute the rotation matrix for rectangular source
 * coordinates for precession and nutation using the conventions for
 * the J2000 system.
 *   The precession series is from the Supplement to the 1984
 * Astronomical Almanac except for the terms independent of T which
 * were taken from the 1986 Astronomical Almanac to give agreement
 * with the tabulated values.
 * 
 * Adapted from the AIPSish JPRENU.FOR
 * \param dir     1 => convert from mean to apparent;
 *               -1 => convert from apparent to mean.
 * \param JD      Julian date of observation (e.g. 2446754.123445)
 * \param equin   Epoch of mean equinox (e.g. 2000.0)
 * \param doNut   If true include nutation terms.
 * \param rotMat  Rotation matrix
 */
static void jprenu(olong dir, double JD, ofloat equin, gboolean doNut, 
		   odouble rotMat[3][3])
{
  odouble delpsi, deleps, pr[3][3], nr[3][3], zeta,
    z, theta, t, tl, con, jdref, twopi, cc, czeta, szeta, cz, sz,
    ctheta, stheta, eps, ceps, seps;

  twopi = 2.0 * G_PI;

  cc = twopi / 1296000.0;
  jdref = 1721045.0 + equin * 365.25; /* reference JD (mean equinox) */

  /* general precession */
  t = (jdref - 2451545.0) / 36525.0;
  tl = dir * (JD - jdref) / 36525;
  con = 2306.21796 + (1.39656 - 0.000139*t) * t;
  zeta = cc * (((0.0180 * tl) + (0.30204 - 0.000344*t)) * tl +  con) * tl;
  z = cc * (((0.01836 * tl) + 1.09476 + 0.000066*t) * tl + con) * tl;
  theta = cc * (((-.04176 * tl) - (0.4266 + 0.000217*t)) * tl +
		2004.3108 - (0.85330 + 0.000217*t) * t) * tl;
  czeta = cos (zeta);
  szeta = sin (zeta);
  cz = cos (z);
  sz = sin (z);
  ctheta = cos (theta);
  stheta = sin (theta);

  /* precession rotation matrix */
  pr[0][0] = czeta*ctheta*cz - szeta*sz;
  pr[1][0] = -szeta*ctheta*cz - czeta*sz;
  pr[2][0] = -stheta*cz;
  pr[0][1] = ctheta*ctheta*sz + szeta*cz;
  pr[1][1] = -szeta*ctheta*sz + czeta*cz;
  pr[2][1] = -stheta*sz;
  pr[0][2] = czeta*stheta;
  pr[1][2] = -szeta*stheta;
  pr[2][2] = ctheta;

  /* Nutation */
  if (doNut) {
    jnut (JD, &delpsi, &deleps);
    delpsi = dir * delpsi;
    deleps = dir * deleps;
    eps =  cc * 3600.0 * ((((5.04e-7 * t) - 1.6e-7) * t - 0.0130042) * t + 23.439291);
    ceps = cos (eps);
    seps = sin (eps);
    nr[0][0] = 1.0;
    nr[1][0] = -delpsi * ceps;
    nr[2][0] = -delpsi * seps;
    nr[0][1] = delpsi * ceps;
    nr[1][1] = 1.0;
    nr[2][1] = -deleps;
    nr[0][2] = delpsi * seps;
    nr[1][2] = deleps;
    nr[2][2] = 1.0;

    /* Nutation * Precession */
    rotMat[0][0] =         pr[0][0]  + nr[1][0]*pr[0][1] + nr[2][0]*pr[0][2];
    rotMat[1][0] = nr[0][1]*pr[0][0] +         pr[0][1]  + nr[2][1]*pr[0][2];
    rotMat[2][0] = nr[0][2]*pr[0][0] + nr[1][2]*pr[0][1] + pr[0][2];
    rotMat[0][1] =         pr[1][0]  + nr[1][0]*pr[1][1] + nr[2][0]*pr[1][2];
    rotMat[1][1] = nr[0][1]*pr[1][0] +         pr[1][1]  + nr[2][1]*pr[1][2];
    rotMat[2][1] = nr[0][2]*pr[1][0] + nr[1][2]*pr[1][1] + pr[1][2];
    rotMat[0][2] =         pr[2][0]  + nr[1][0]*pr[2][1] + nr[2][0]*pr[2][2];
    rotMat[1][2] = nr[0][1]*pr[2][0] +         pr[2][1]  + nr[2][1]*pr[2][2];
    rotMat[2][2] = nr[0][2]*pr[2][0] + nr[1][2]*pr[2][1] + pr[2][2];

  } else {
    /* Precession only */
    rotMat[0][0] = pr[0][0];
    rotMat[0][1] = pr[0][1];
    rotMat[0][2] = pr[0][2];
    rotMat[1][0] = pr[1][0];
    rotMat[1][1] = pr[1][1];
    rotMat[1][2] = pr[1][2];
    rotMat[2][0] = pr[2][0];
    rotMat[2][1] = pr[2][1];
    rotMat[2][2] = pr[2][2];
  } 
} /* end jprenu */

/**
 * Routine to compute vectors needed for corrections to true positions
 * due to aberation  and gravitational deflection of light by the sun.
 * Approximations and formulae mostly from Astronomical Almanac and
 * Starlink routine SLAEVP.
 * 
 * Adapted from the AIPSish JABER.FOR
 * \param JD      Julian date of observation (e.g. 2446754.123445)
 * \param equin   Epoch of mean equinox (e.g. 2000.0)
 * \param diurn   true if diurnal aberation correction desired.
 * \param rhoGeo  geocentric radius to the point of observation,
 *                expressed in units of Earth's equatorial radius.
 *                This is used only when diurnal aberration is
 *                applied.
 * \param phiGeo geocentric latitude of the observatory.  This is
 *                used only when diurnal aberration is applied.
 *                (Radians).
 * \param rLST    local sidereal time.  This is used only when
 *                diurnal aberration is applied. (Radians)
 * \param poso    Heliocentric position of observer (AU)
 * \param velo    Barycentric velocity of observer (AU/day)
 */
static void jaber(odouble JD, ofloat equin, gboolean diurn, 
		  odouble rhoGeo, odouble phiGeo, odouble rLST,
		  odouble poso[3], odouble velo[3])
{
  odouble dvb[3], dpb[3], dvh[3], dph[3], con;
  /* eqrau = earth's equatorial radius in au. */
  odouble eqrau = 6378140.0 / 1.49597870e11;
  /* omega = earth's rotation rate in  rotations per day. */
  odouble omega = 360.9856 / 360.0;

  /* earth positions and motions. */
  EarthEphem (JD, equin, dvb, dpb, dvh, dph);

  /* convert units */
  velo[0] = dvb[0] * 86400.0;
  velo[1] = dvb[1] * 86400.0;
  velo[2] = dvb[2] * 86400.0;
  poso[0] = dph[0];
  poso[1] = dph[1];
  poso[2] = dph[2];

  /* diurnal effects */
  if (diurn) {
    /* position */
    con = eqrau * rhoGeo;
    poso[0] = poso[0] + con * cos (phiGeo) * cos (rLST);
    poso[1] = poso[1] + con * cos (phiGeo) * sin (rLST);
    poso[2] = poso[2] + con * sin (phiGeo);
    /* velocity */
    con = con * omega;
    velo[0] = velo[0] - con * cos (phiGeo) * sin (rLST);
    velo[1] = velo[1] + con * cos (phiGeo) * cos (rLST);
  } 
} /* end jaber */
		  
/**
 * Routine to correct apparent positions for the position of the
 * Earth's pole.  Method from Astronomical Almanac 1986 B59.
 * 
 * Adapted from the AIPSish JPOLAR.FOR
 * \param dir     1 => Celestial to terrestrial
 *               -1 => Terrestrial to celestial 
 * \param polarx Polar "X" position (arcsec)
 * \param polary Polar "Y" position (arcsec)
 * \param pos    Position in rectangular coordinates (angle cosines)
 */
static void jpolar(olong dir, ofloat polarx, ofloat polary, odouble pos[3])
{
  odouble cc, cx, sx, cy, sy, t1, t2;

  cc = 4.848136811e-6;
  cx = cos (cc * polarx);
  sx = sin (cc * polarx);
  cy = cos (cc * polary);
  sy = sin (cc * polary);

  /* direction? */
  if (dir > 0) {
    /* celestial to terrestial */
    /* y correction. */
    t1 = cy * pos[1] - sy * pos[2];
    t2 = sy * pos[1] + cy * pos[2];
    pos[1] = t1;
    pos[2] = t2;

    /* x correction. */
    t1 = cx * pos[0] + sx * pos[2];
    t2 = -sx * pos[0] + cx * pos[2];
    pos[0] = t1;
    pos[2] = t2;

  } else {
    /* terrestial to celestial */
    /* x correction. */
    t1 = cx * pos[0] - sx * pos[2];
    t2 = sx * pos[0] + cx * pos[2];
    pos[0] = t1;
    pos[2] = t2;

    /* y correction. */
    t1 = cy * pos[1] + sy * pos[2];
    t2 = -sy * pos[1] + cy * pos[2];
    pos[1] = t1;
    pos[2] = t2;
  } 
} /* end jpolar */

/**
 *  Calculation of nutation, assuming a non-rigid Earth.
 * Uses IAU 1980 precession series from Supplement to the Astronomical
 * Almanac 1984.  For Epoch J2000.
 * 
 * Adapted from the AIPSish JNUT.FOR
 * \param JD     Julian date.
 * \param delPsi Nutation in longitude (radians)
 * \param delEps Nutation in obliquity
 */
static void jnut(odouble JD, odouble *delPsi, odouble *delEps)
{
  olong   i;
  odouble arg[106], delpsi, deleps, mamoon, masun, f, d, omega, t, twopi, cc;
  
  odouble along[106] = {
    -17.1996, 0.2062, 0.0046, 0.0011, -.0003,
    -.0003, -.0002, 0.0001, -1.3187, 0.1426, -.0517,
    0.0217, 0.0129, 0.0048, -.0022, 0.0017, -.0015,
    -.0016, -.0012, -.0006, -.0005, 0.0004, 0.0004,
    -.0004, 0.0001, 0.0001, -.0001, 0.0001, 0.0001,
    -.0001, -.2274, 0.0712, -.0386, -.0301, -.0158,
    0.0123, 0.0063, 0.0063, -.0058, -.0059, -.0051,
    -.0038, 0.0029, 0.0029, -.0031, 0.0026, 0.0021,
    0.0016, -.0013, -.0010, -.0007, 0.0007, -.0007,
    -.0008, 0.0006, 0.0006, -.0006, -.0007, 0.0006,
    -.0005, 0.0005, -.0005, -.0004, 0.0004, -.0004,
    -.0003, 0.0003, -.0003, -.0003, -.0002, -.0003,
    -.0003, 0.0002, -.0002, 0.0002, -.0002, 0.0002,
    0.0002, 0.0001, -.0001, 0.0001, -.0002, -.0001,
    0.0001, -.0001, -.0001, 0.0001, 0.0001, 0.0001,
    -.0001, -.0001, 0.0001, 0.0001, -.0001, 0.0001,
    0.0001, -.0001, -.0001, -.0001, -.0001, -.0001,
    -.0001, -.0001, 0.0001, -.0001, 0.0001};

  double talong[106] = {
    -174.2e-4, 0.2e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.6e-4, -3.4e-4,
    1.2e-4, -.5e-4, 0.1e-4, 0.0e-4, 0.0e-4, -.1e-4, 0.0e-4, 0.1e-4,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    -0.2e-4, 0.1e-4, -.4e-4, 4*0.0e-4, 0.1e-4, -0.1e-4,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  odouble aobli[106] = {
    9.2025, -.0895, -.0024, 0.0, 0.0001, 0.0,
    0.0001, 0.0, 0.5736, 0.0054, 0.0224, -.0095,
    -.0070, 0.0001, 0.0, 0.0, 0.0009, 0.0007, 0.0006,
    0.0003, 0.0003, -.0002, -.0002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0977,
    -.0007, 0.0200, 0.0129, -.0001, -.0053, -.0002,
    -.0033, 0.0032, 0.0026, 0.0027, 0.0016, -.0001,
    -.0012, 0.0013, -.0001, -.0010, -.0008, 0.0007,
    0.0005, 0.0000, -.0003, 0.0003, 0.0003, 0.0000,
    -.0003, 0.0003, 0.0003, -.0003, 0.0003, 0.0000,
    0.0003,  0.00, 0.0, 0.0, 0.0, 0.0, 
    .0001, 0.0001, 0.0001, 0.0001, 0.0001,
    -.0001, 0.0001, -.0001,
    0.0001, 0.0000, -.0001, -.0001, 0.0000, -.0001,
    0.0001, 0.0000, -.0001, 0.0001, 0.0001, 0.0000,
    0.0000, -.0001, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  odouble taobli[106] = {
    8.9e-4, 0.5e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.1e-4, -.1e-4, -.6e-4, 0.3e-4, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    -.5e-4, 0.0, 0.0e-4, -.1e-4, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0};

  olong nlong  = 106;
  olong nobli  = 89;

  /* time in centuries */
  t = (JD - 2451545.) / 36525.;
  twopi = 2.0 * G_PI;
  cc = twopi / 1296000.;
  
  /* Arguments of the series */
  mamoon = fmod ((((((0.064 * t) + 31.310) * t + 715922.633 
		    + 1325.0 * 1296000.0) * t + 485866.733) * cc), twopi);

  masun = fmod ((((((-.012 * t) - 0.577) * t + 1292581.224
		   + 99.0 * 1296000.0) * t + 1287099.804) * cc), twopi);

  f = fmod ((((((0.011 * t) - 13.257) * t + 295263.137
	       + 1342.0 * 1296000.0) * t + 335778.877) * cc), twopi);

  d = fmod ((((((0.019 * t) - 6.891) * t + 1105601.328
	      + 1236.0 * 1296000.0) * t + 1072261.307) * cc), twopi);
  
  omega = fmod ((((((0.008 * t) + 7.455) * t - 482890.539
		   - 5.0 * 1296000.0) * t + 450160.280) * cc), twopi);
  
/* form final arguments */
  arg[0] = omega;
  arg[1] = 2*omega;
  arg[2] = -2*mamoon + 2*f + omega;
  arg[3] = 2*mamoon - 2*f;
  arg[4] = -2*mamoon + 2*f + 2*omega;
  arg[5] = mamoon - masun - d;
  arg[6] = -2*masun + 2*f - 2*d + omega;
  arg[7] = 2*mamoon - 2*f + omega;
  arg[8] = 2*f - 2*d + 2*omega;
  arg[9] = masun;
  arg[10] = masun + 2*f - 2*d + 2*omega;
  arg[11] = -masun + 2*f - 2*d + 2*omega;
  arg[12] = 2*f - 2*d + omega;
  arg[13] = 2*mamoon - 2*d;
  arg[14] = 2*f - 2*d;
  arg[15] = 2*masun;
  arg[16] = masun + omega;
  arg[17] = 2*masun + 2*f - 2*d + 2*omega;
  arg[18] = -masun + omega;
  arg[19] = -2*mamoon + 2*d + omega;
  arg[20] = -masun + 2*f - 2*d + omega;
  arg[21] = 2*mamoon - 2*d + omega;
  arg[22] = masun + 2*f - 2*d + omega;
  arg[23] = mamoon - d;
  arg[24] = 2*mamoon + masun - 2*d;
  arg[25] = -2*f + 2*d + omega;
  arg[26] = masun - 2*f + 2*d;
  arg[27] = masun + 2*omega;
  arg[28] = -mamoon + d + omega;
  arg[29] = masun + 2*f - 2*d;
  arg[30] = 2*f + 2*omega;
  arg[31] = mamoon;
  arg[32] = 2*f + omega;
  arg[33] = mamoon + 2*f + 2*omega;
  arg[34] = mamoon - 2*d;
  arg[35] = -mamoon + 2*f + 2*omega;
  arg[36] = 2*d;
  arg[37] = mamoon + omega;
  arg[38] = -mamoon + omega;
  arg[39] = -mamoon + 2*f + 2*d + 2*omega;
  arg[40] = mamoon + 2*f + omega;
  arg[41] = 2*f + 2*d + 2*omega;
  arg[42] = 2*mamoon;
  arg[43] = mamoon + 2*f - 2*d + 2*omega;
  arg[44] = 2*mamoon + 2*f + 2*omega;
  arg[45] = 2*f;
  arg[46] = -mamoon + 2*f + omega;
  arg[47] = -mamoon + 2*d + omega;
  arg[48] = mamoon - 2*d + omega;
  arg[49] = -mamoon + 2*f + 2*d + omega;
  arg[50] = mamoon + masun - 2*d;
  arg[51] = masun + 2*f + 2*omega;
  arg[52] = -masun + 2*f + 2*omega;
  arg[53] = mamoon + 2*f + 2*d + 2*omega;
  arg[54] = mamoon + 2*d;
  arg[55] = 2*mamoon + 2*f - 2*d + 2*omega;
  arg[56] = 2*d + omega;
  arg[57] = 2*f + 2*d + omega;
  arg[58] = mamoon + 2*f - 2*d + omega;
  arg[59] = -2*d + omega;
  arg[60] = mamoon - masun;
  arg[61] = 2*mamoon + 2*f + omega;
  arg[62] = masun - 2*d;
  arg[63] = mamoon - 2*f;
  arg[64] = d;
  arg[65] = mamoon + masun;
  arg[66] = mamoon + 2*f;
  arg[67] = mamoon - masun + 2*f + 2*omega;
  arg[68] = -mamoon - masun + 2*f + 2*d + 2*omega;
  arg[69] = -2*mamoon + omega;
  arg[70] = 3*mamoon + 2*f + 2*omega;
  arg[71] = -masun + 2*f + 2*d + 2*omega;
  arg[72] = mamoon + masun + 2*f + 2*omega;
  arg[73] = -mamoon + 2*f - 2*d + omega;
  arg[74] = 2*mamoon + omega;
  arg[75] = mamoon + 2*omega;
  arg[76] = 3*mamoon;
  arg[77] = 2*f + d + 2*omega;
  arg[78] = -mamoon + 2*omega;
  arg[79] = mamoon - 4*d;
  arg[80] = -2*mamoon + 2*f + 2*d + 3*omega;
  arg[81] = -mamoon + 2*f + 4*d + 2*omega;
  arg[82] = 2*mamoon - 4*d;
  arg[83] = mamoon + masun + 2*f - 2*d + 2*omega;
  arg[84] = mamoon + 2*f + 2*d + omega;
  arg[85] = -2*mamoon + 2*f + 4*d + 2*omega;
  arg[86] = -mamoon + 4*f + 2*omega;
  arg[87] = mamoon - masun -2*d;
  arg[88] = 2*mamoon + 2*f - 2*d + omega;
  arg[89] = 2*mamoon + 2*f + 2*d + 2*omega;
  arg[90] = mamoon + 2*d + omega;
  arg[91] = 4*f - 2*d + 2*omega;
  arg[92] = 3*mamoon + 2*f - 2*d + 2*omega;
  arg[93] = mamoon + 2*f - 2*d;
  arg[94] = masun + 2*f + omega;
  arg[95] = -mamoon -masun + 2*d + omega;
  arg[96] = -2*f + omega;
  arg[97] = 2*f - d + 2*omega;
  arg[98] = masun + 2*d;
  arg[99] = mamoon - 2*f - 2*d;
  arg[100] = -masun + 2*f + omega;
  arg[101] = mamoon + masun -2*d + omega;
  arg[102] = mamoon - 2*f + 2*d;
  arg[103] = 2*mamoon + 2*d;
  arg[104] = 2*f + 4*d + 2*omega;
  arg[105] = masun + d;

  /* sum in arc sec. */
  /* longitude: */
  delpsi = 0.0;
  for (i= 0; i<nlong; i++) delpsi = delpsi + (along[i] + t * talong[i]) * sin (arg[i]);

  /* obliquity */
  deleps = 0.0;
  for (i= 0; i<nobli; i++) deleps = deleps + (aobli[i] + t * taobli[i]) * cos (arg[i]);

  /* convert to radians */
  *delPsi = delpsi * cc;
  *delEps = deleps * cc;
} /*  end jnut */

 /**
  * Barycentric and heliocentric velocity and position of the Earth.
  *  Accuracy:
  *    The maximum deviations from the JPL DE96 ephemeris are as
  *     follows:
  *     barycentric velocity                  42  cm/s
  *     barycentric position           0.000 046  au
  *     heliocentric velocity                 42  cm/s
  *     heliocentric position          0.000 011  au
  * This routine is adapted from the BARVEL and BARCOR
  * subroutines of P.Stumpff, which are described in Astron. Astrophys. Suppl. Ser. 
  * 41, 1-8 (1980).  Most of the changes are merely cosmetic and do not affect the 
  * results at all.  However, some adjustments have been made so as to give
  * results that refer to the new (IAU 1976 'FK5') equinox  and precession, 
  * although the differences these changes make relative to the results from Stumpff's 
  * original 'FK4' version are smaller than the inherent accuracy of the algorithm.  
  * One  minor shortcoming in the original routines that has NOT been
  * corrected is that better numerical accuracy could be achieved
  * if the various polynomial evaluations were nested.  Note also
  * that one of Stumpff's precession constants differs by 0.001 arcsec
  * from the value given in the Explanatory Supplement to the A.E.
  * P T Wallace   Starlink   March 1986
  * Adapted for use in AIPS.
  * Adapted from fortran SLAEVP.FOR
  * \param JD      Julian date of observation (e.g. 2446754.123445)
  * \param equin   Epoch of mean equinox (e.g. 2000.0)and
  *                equinox of the vectors returned.  If equin <= 0,
  *                all vectors are referred to the mean equator and
  *                equinox (FK5) of date JD.
  * \param dvb     barycentric velocity (AU/s)
  * \param dpb     barycentric position (AU)
  * \param dvh     heliocentric velocity (AU/s)
  * \param dph     heliocentric position (AU)
  */
static void EarthEphem(odouble JD, ofloat equin, odouble dvb[3], odouble dpb[3], 
		       odouble dvh[3], odouble dph[3])
{
  olong ideq, i, j, k;
  ofloat t, tsq, a, pertl=0.0, pertld, pertr, pertrd, cosa, sina, esq, e, param, twoe,
    twog, g, phi=0.0, f, sinf, cosf, phid, psid, pertp, pertpd, tl,
    sinlm, coslm, sigma, b, plon, pomg, pecc, flatm, flat;
  odouble dt, dtsq, dlocal, dml=0.0, deps, dparam, dpsi, d1pdro, drd, drld, dtl, dsinls,
    dcosls, dxhd, dyhd, dzhd, dxbd, dybd, dzbd, dcosep, dsinep,
    dyahd, dzahd, dyabd, dzabd, dr, dxh, dyh, dzh, dxb, dyb, dzb,
    dyah, dzah, dyab, dzab, depj, deqcor;
  ofloat sn[4], forbel[7], sorbel[17], sinlp[4], coslp[4];
  odouble dprema[3][3], w, vw[3];
  odouble dc2pi=6.2831853071796;
  ofloat  cc2pi=6.283185;
  odouble ds2r=0.7272205216643e-04;
  odouble dcsld=1.990987e-07;
  ofloat ccsgd=1.990969e-07;

  /* some constants used in the calculation of the lunar  contribution. */
  ofloat cckm=3.122140e-05;
  ofloat ccmld=2.661699e-06;
  ofloat ccfdi=2.399485e-07;

  /* constants dcfel(i,k) of fast changing elements */
  odouble dcfel[8][3] = {
    {1.7400353e+00, 6.2833195099091e+02, 5.2796e-06},
    {6.2565836e+00, 6.2830194572674e+02,-2.6180e-06},
    {4.7199666e+00, 8.3997091449254e+03,-1.9780e-05},
    {1.9636505e-01, 8.4334662911720e+03,-5.6044e-05},
    {4.1547339e+00, 5.2993466764997e+01, 5.8845e-06},
    {4.6524223e+00, 2.1354275911213e+01, 5.6797e-06},
    {4.2620486e+00, 7.5025342197656e+00, 5.5317e-06},
    {1.4740694e+00, 3.8377331909193e+00, 5.6093e-06}};

  /* constants dceps and ccsel(i,k)  of slowly changing elements */
  /* i=1           i=2           i=3 */
  odouble dceps[3] = {
    4.093198e-01,-2.271110e-04,-2.860401e-08};

  ofloat ccsel[17][3] = {
    {1.675104e-02,-4.179579e-05,-1.260516e-07},
    {2.220221e-01, 2.809917e-02, 1.852532e-05},
    {1.589963e+00, 3.418075e-02, 1.430200e-05},
    {2.994089e+00, 2.590824e-02, 4.155840e-06},
    {8.155457e-01, 2.486352e-02, 6.836840e-06},
    {1.735614e+00, 1.763719e-02, 6.370440e-06},
    {1.968564e+00, 1.524020e-02,-2.517152e-06},
    {1.282417e+00, 8.703393e-03, 2.289292e-05},
    {2.280820e+00, 1.918010e-02, 4.484520e-06},
    {4.833473e-02, 1.641773e-04,-4.654200e-07},
    {5.589232e-02,-3.455092e-04,-7.388560e-07},
    {4.634443e-02,-2.658234e-05, 7.757000e-08},
    {8.997041e-03, 6.329728e-06,-1.939256e-09},
    {2.284178e-02,-9.941590e-05, 6.787400e-08},
    {4.350267e-02,-6.839749e-05,-2.714956e-07},
    {1.348204e-02, 1.091504e-05, 6.903760e-07},
    {3.106570e-02,-1.665665e-04,-1.590188e-07}};

  /* constants of the arguments of  the short-period perturbations */
  /* by the planets:   dcargs(i,k). */
  /* i=1               i=2 */
  odouble dcargs[15][2] = {
    {5.0974222e+00,-7.8604195454652e+02},
    {3.9584962e+00,-5.7533848094674e+02},
    {1.6338070e+00,-1.1506769618935e+03},
    {2.5487111e+00,-3.9302097727326e+02},
    {4.9255514e+00,-5.8849265665348e+02},
    {1.3363463e+00,-5.5076098609303e+02},
    {1.6072053e+00,-5.2237501616674e+02},
    {1.3629480e+00,-1.1790629318198e+03},
    {5.5657014e+00,-1.0977134971135e+03},
    {5.0708205e+00,-1.5774000881978e+02},
    {3.9318944e+00, 5.2963464780000e+01},
    {4.8989497e+00, 3.9809289073258e+01},
    {1.3097446e+00, 7.7540959633708e+01},
    {3.5147141e+00, 7.9618578146517e+01},
    {3.5413158e+00,-5.4868336758022e+02}};

  /* amplitudes ccamps(n,k) of the short-period perturbations. */
  /* n=1          n=2          n=3          n=4          n=5 */
  ofloat ccamps[15][5] = {
    {-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5,-2.490817e-7},
    {-3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5,-1.823138e-7},
    { 6.593466e-7, 1.322572e-5, 9.258695e-6,-4.674248e-7,-3.646275e-7},
    { 1.140767e-5,-2.049792e-5,-4.747930e-6,-2.638763e-6,-1.245408e-7},
    { 9.516893e-6,-2.748894e-6,-1.319381e-6,-4.549908e-6,-1.864821e-7},
    { 7.310990e-6,-1.924710e-6,-8.772849e-7,-3.334143e-6,-1.745256e-7},
    {-2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6,-1.655307e-7},
    {-3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6,-3.736225e-7},
    { 3.442177e-7, 2.671323e-6, 1.832858e-6,-2.394688e-7,-3.478444e-7},
    { 8.702406e-6,-8.421214e-6,-1.372341e-6,-1.455234e-6,-4.998479e-8},
    {-1.488378e-6,-1.251789e-5, 5.226868e-7,-2.049301e-7, 0.0e0},
    {-8.043059e-6,-2.991300e-6, 1.473654e-7,-3.154542e-7, 0.0e0},
    { 3.699128e-6,-3.316126e-6, 2.901257e-7, 3.407826e-7, 0.0e0},
    { 2.550120e-6,-1.241123e-6, 9.901116e-8, 2.210482e-7, 0.0e0},
    {-6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.0e0}};

  /* constants of the secular perturbations in longitude */
  /* ccsec3 and ccsec(n,k). */
  /* n=1           n=2           n=3 */
  ofloat ccsec3 = -7.757020e-08;
  ofloat ccsec[4][3] = {
    {1.289600e-06, 5.550147e-01, 2.076942e+00},
    {3.102810e-05, 4.035027e+00, 3.525565e-01},
    {9.124190e-06, 9.990265e-01, 2.622706e+00},
    {9.793240e-07, 5.508259e+00, 1.559103e+01}};

  /* sidereal rate dcsld in  longitude, rate ccsgd in mean  anomaly. 
     Constants dcargm(i,k)  of the arguments of the perturbations of the motion of 
     the moon. */
  /*     i=1               i=2 */
  odouble dcargm[3][2] = {
    {5.1679830e+00, 8.3286911095275e+03},
    {5.4913150e+00,-7.2140632838100e+03},
    {5.9598530e+00, 1.5542754389685e+04}};

  /* amplitudes ccampm(n,k) of the  perturbations of the moon. */
  /* n=1          n=2           n=3           n=4 */
  ofloat ccampm[3][4] = {
    { 1.097594e-01, 2.896773e-07, 5.450474e-02, 1.438491e-07},
    {-2.223581e-02, 5.083103e-08, 1.002548e-02,-2.291823e-08},
    { 1.148966e-02, 5.658888e-08, 8.249439e-03, 4.063015e-08}};

  /* ccpamv(k) = a*m*dl/dt (planets),  dc1mme=1-mass(earth+moon) */
  ofloat ccpamv[4] = {8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12};
  odouble dc1mme=0.99999696;

  /* ccpam(k)=a*m(planets), ccim=inclination(moon) */
  ofloat ccpam[4] = {4.960906e-3, 2.727436e-3, 8.392311e-4, 1.556861e-3};
  ofloat ccim=8.978749e-2;

  /* Check parameter ideq, and time  arguments . */
  ideq = 0;
  if (equin > 0) ideq = 1;
  dt = (JD-2415020.0) / 36525;
  t =  (ofloat) (dt);
  dtsq = dt*dt;
  tsq =  (ofloat) (dtsq);

  /* Values of all elements for the desired date. */
  for (k= 0; k<8; k++) { /* loop 10 */
    dlocal = fmod ((dcfel[k][0] + dt*dcfel[k][1] + dtsq*dcfel[k][2]), dc2pi);
    if (k == 0) {
      dml = dlocal;
    } else {
      forbel[k-1] =  (ofloat) (dlocal);
    } 
  } /* end loop  L10:  */

  deps = fmod ((dceps[0] + dt*dceps[1] + dtsq*dceps[2]), dc2pi);
  for (k= 0; k<17; k++) 
    sorbel[k] = fmod ((ccsel[k][0] + t*ccsel[k][1] + tsq*ccsel[k][2]), cc2pi);

  /* Replace equivalences */
  e = sorbel[0];
  g = forbel[0];

  /* Secular perturbations in  longitude. */
  for (k= 0; k<4; k++) { /* loop 30 */
    a =  fmod ((ccsec[k][1] + t*ccsec[k][2]), cc2pi);
    sn[k] = sin (a);
  } /* end loop  L30:  */

  /* Periodic perturbations of the emb (earth-moon barycenter). */
  pertl = ccsec[0][0] * sn[0] + ccsec[1][0] * sn[1] +
    (ccsec[2][0] + t * (ccsec3) ) * sn[2] + ccsec[3][0] * sn[3];
  pertld = 0.0;
  pertr  = 0.0;
  pertrd = 0.0;
  for (k= 0; k<15; k++) { /* loop 40 */
    a = fmod ((dcargs[k][0] + dt*dcargs[k][1]), dc2pi);
    cosa = cos (a);
    sina = sin (a);
    pertl = pertl  +  ccamps[k][0]*cosa + ccamps[k][1]*sina;
    pertr = pertr  +  ccamps[k][2]*cosa + ccamps[k][3]*sina;
    if (k < 10) {
      pertld = pertld + (ccamps[k][1]*cosa - ccamps[k][0]*sina)*ccamps[k][4];
      pertrd = pertrd + (ccamps[k][3]*cosa - ccamps[k][2]*sina)*ccamps[k][4];
    } 
  } /* end loop  L40:  */

  /* Elliptic part of the motion of  the emb. */
  esq = e*e;
  dparam = 1 - (odouble) (esq);
  param  =  (ofloat) (dparam);
  twoe   = e + e;
  twog   = g + g;
  phi    = twoe*((1.0 - esq*0.125)*sin (g) + e*0.625*sin (twog) + esq*0.5416667*sin (g + twog) );
  f      = g + phi;
  sinf   = sin (f);
  cosf   = cos (f);
  dpsi   = dparam  /  (1 + (odouble) (e*cosf));
  phid   = twoe*ccsgd*((1.0 + esq*1.5)*cosf + e*(1.25 - sinf*sinf*0.5));
  psid   = ccsgd*e*sinf / sqrt(param);

  /* Perturbed heliocentric motion of  the EMB. */
  d1pdro = 1 + (odouble) (pertr);
  drd    = d1pdro*((odouble) (psid) + dpsi*(odouble) (pertrd));
  drld   = d1pdro*dpsi*(dcsld + (odouble) (phid) + (odouble) (pertld));
  dtl    = fmod ((dml + (odouble) (phi) + (odouble) (pertl)), dc2pi);
  dsinls = sin (dtl);
  dcosls = cos (dtl);
  dxhd   = drd*dcosls - drld*dsinls;
  dyhd   = drd*dsinls + drld*dcosls;

  /* Influence of eccentricity, evection and variation on the geocentric motion of the moon. */
  pertl = 0.0;
  pertld = 0.0;
  pertp = 0.0;
  pertpd = 0.0;
  for (k= 0; k<3; k++) { /* loop 50 */
    a = fmod ((dcargm[k][0] + dt*dcargm[k][1]), dc2pi);
    sina = sin (a);
    cosa = cos (a);
    pertl  = pertl  + ccampm[k][0]*sina;
    pertld = pertld + ccampm[k][1]*cosa;
    pertp  = pertp  + ccampm[k][2]*cosa;
    pertpd = pertpd - ccampm[k][3]*sina;
  } /* end loop  L50:  */

  /* Heliocentric motion of the  earth. */
  tl = forbel[1] + pertl;
  sinlm = sin (tl);
  coslm = cos (tl);
  sigma = cckm / (1.0 + pertp);
  a = sigma*(ccmld + pertld);
  b = sigma*pertpd;
  dxhd = dxhd + (odouble) (a*sinlm) + (odouble) (b*coslm);
  dyhd = dyhd - (odouble) (a*coslm) + (odouble) (b*sinlm);
  dzhd =      - (odouble) (sigma*ccfdi*cos (forbel[2]));

  /* Barycentric motion of the earth. */
  dxbd = dxhd*dc1mme;
  dybd = dyhd*dc1mme;
  dzbd = dzhd*dc1mme;
  for (k= 1; k<=4; k++) { /* loop 60 */
    plon = forbel[k+2];
    pomg = sorbel[k];
    pecc = sorbel[k+8];
    tl = fmod (plon + 2.0*pecc*sin (plon - pomg), cc2pi);
    sinlp[k-1] = sin (tl);
    coslp[k-1] = cos (tl);
    dxbd = dxbd + (odouble) (ccpamv[k-1]*(sinlp[k-1] + pecc*sin (pomg)));
    dybd = dybd - (odouble) (ccpamv[k-1]*(coslp[k-1] + pecc*cos (pomg)));
    dzbd = dzbd - (odouble) (ccpamv[k-1]*sorbel[k+12]* cos (plon - sorbel[k+4]));
  } /* end loop  L60:  */

  /* Transition to mean equator of  date. */
  dcosep = cos (deps);
  dsinep = sin (deps);
  dyahd = dcosep*dyhd - dsinep*dzhd;
  dzahd = dsinep*dyhd + dcosep*dzhd;
  dyabd = dcosep*dybd - dsinep*dzbd;
  dzabd = dsinep*dybd + dcosep*dzbd;

  /* Heliocentric coordinates of the  earth. */
  dr = dpsi*d1pdro;
  flatm = ccim*sin (forbel[2]);
  a = sigma*cos (flatm);
  dxh = dr*dcosls - (odouble) (a*coslm);
  dyh = dr*dsinls - (odouble) (a*sinlm);
  dzh =           - (odouble) (sigma*sin (flatm));

  /* Barycentric coordinates of the  earth. */
  dxb = dxh*dc1mme;
  dyb = dyh*dc1mme;
  dzb = dzh*dc1mme;
  for (k= 1; k<=4; k++) { /* loop 70 */
    flat = sorbel[k+12]*sin (forbel[k+2] - sorbel[k+4]);
    a = ccpam[k-1]*(1.0 - sorbel[k+8]*cos (forbel[k+2] - sorbel[k]));
    b = a*cos (flat);
    dxb = dxb - (odouble) (b*coslp[k-1]);
    dyb = dyb - (odouble) (b*sinlp[k-1]);
    dzb = dzb - (odouble) (a*sin (flat));
  } /* end loop  L70:  */

  /* Transition to mean equator of  date. */
  dyah = dcosep*dyh - dsinep*dzh;
  dzah = dsinep*dyh + dcosep*dzh;
  dyab = dcosep*dyb - dsinep*dzb;
  dzab = dsinep*dyb + dcosep*dzb;

  /* Copy result components into  vectors, correcting for FK5  equinox. */
  depj = 2000.0 + ((JD - 2451545.0) / 365.25);
  deqcor = ds2r*(0.035 + 0.00085*(depj - 1950));
  dvh[0] = dxhd - deqcor*dyahd;
  dvh[1] = dyahd + deqcor*dxhd;
  dvh[2] = dzahd;
  dvb[0] = dxbd - deqcor*dyabd;
  dvb[1] = dyabd + deqcor*dxbd;
  dvb[2] = dzabd;
  dph[0] = dxh - deqcor*dyah;
  dph[1] = dyah + deqcor*dxh;
  dph[2] = dzah;
  dpb[0] = dxb - deqcor*dyab;
  dpb[1] = dyab + deqcor*dxb;
  dpb[2] = dzab;

  /* Was precession to another equinox requested? */
  if (ideq != 0) {
    /* yes: compute precession matrix from mjd date to julian epoch equin. */
    jprenu (-1, JD, equin, FALSE, dprema);

    /* rotate dvh */
    for (j= 0; j<3; j++) { /* loop 110 */
      w = 0;
      for (i= 0; i<3; i++) { /* loop 100 */
	w = w + dprema[i][j]*dvh[i];
      } /* end loop  L100: */;
      vw[j] = w;
    } /* end loop  L110: */

    dvh[0] = vw[0];
    dvh[1] = vw[1];
    dvh[2] = vw[2];

    /* Rotate dvb */
    for (j= 0; j<3; j++) { /* loop 140 */
      w = 0;
      for (i= 0; i<3; i++) { /* loop 130 */
	w = w + dprema[i][j]*dvb[i];
      } /* end loop  L130: */
      vw[j] = w;
    } /* end loop  L140: */

    dvb[0] = vw[0];
    dvb[1] = vw[1];
    dvb[2] = vw[2];

    /* Rotate dph */
    for (j= 0; j<3; j++) { /* loop 170 */
      w = 0;
      for (i= 0; i<3; i++) { /* loop 160 */
	/*	w = w + dprema[i][j-1]*dph[i]; ??? original */
	w = w + dprema[i][j]*dph[i];
      } /* end loop  L160: */
      vw[j] = w;
    } /* end loop  L170: */

    dph[0] = vw[0];
    dph[1] = vw[1];
    dph[2] = vw[2];

    /* Rotate dpb */
    for (j= 0; j<3; j++) { /* loop 200 */
      w = 0;
      for (i= 0; i<3; i++) { /* loop 190 */
	w = w + dprema[i][j]*dpb[i];
      } /* end loop  L190: */
      vw[j] = w;
    } /* end loop  L200: */

    dpb[0] = vw[0];
    dpb[1] = vw[1];
    dpb[2] = vw[2];
  } /* end need to precess */
} /* end EarthEphem */
