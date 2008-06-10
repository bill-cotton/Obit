# $Id: TableUtil.py,v 1.1 2006/06/13 18:32:38 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2006
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

# Python package for Obit Table Utilities 
import Obit, Table, OErr

def PCCMerge(inTab, outTab, err):
    """ Merge CC in the same cell

    inTab    = input Python Obit Table
    outTab   = extant output Python Obit Table
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    if not Table.PIsA(outTab):
        raise TypeError,"outTab MUST be a Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    Obit.TableCCUtilMerge(inTab.me, outTab.me, err.me)
    # end PCCMerge

def PConvolve(inFA1, inFA2, err):
    """ Return the convolution of two 2-D arrays

    inFA1 and inFA2 must be compatible size and should not contain magic value blanks.
    Convolved array returned.
    inFA1   = First array
    inFA2   = Second array
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FArray.PIsA(inFA1):
        raise TypeError,"inFA1 MUST be a Python Obit FArray"
    if not FArray.PIsA(inFA2):
        raise TypeError,"inFA2 MUST be a Python Obit FArray"
    # 
    outFA = FArray.FArray("None")
    outFA.me = Obit.FArrayUtilConvolve (inFA1.me, inFA2.me, err.me)
    return outFA
    # end PConvolve

