""" Python Obit multidimensional array of float class

ObitInfoList Linked list of labeled items class.
This facility allows storing arrays of values of the same (native) 
data type and retrieving them by name or order number in the list.
This container class is similar to a Python dictionary and is used
extensively in Obit to pass control information.
Most Obit objects contain an InfoList member.

Members:
Dict      - Python dictionary with contents of InfoList
"""
# $Id: InfoList.py,v 1.20 2008/02/15 19:41:55 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2008
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

# Python shadow class to ObitInfoList class
import Obit, OErr

class InfoListPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.InfoList_me_set(self.this,value)
            return
        if name=="Dict":
            return PSetDict(self, value)
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.InfoList_me_get(self.this)
        if name=="Dict":
            return PGetDict(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C InfoList instance>"
class InfoList(InfoListPtr):
    """ Python Obit InfoList class

    ObitInfoList Linked list of labeled items class.
    This facility allows storing arrays of values of the same (native) 
    data type and retrieving them by name or order number in the list.
    This container class is similar to a Python dictionary and is used
    extensively in Obit to pass control information.
    Most Obit objects contain an InfoList member.
    """
    def __init__(self) :
        self.this = Obit.new_InfoList()
    def __del__(self):
        if Obit!=None:
            Obit.delete_InfoList(self.this)
    
    def set (self, name, value, ttype=None):
        """ Save a value in an InfoList

        Set an entry in an InfoList, possibly redefining its type and dimension
        self     = input Python InfoList
        name     = name of desired entry
        value    = value to save, either a scalar integer, float, boolean or string
                   or a 1D array of one of these types
                   Type and dimensionality determined from value unless ttype is set
        ttype     = data type, "double", "long", None=>type of value
        """
        itype = type(value)
        dim = [1,1,1,1,1,]
        # Local value of the correct type
        tvalue = value
        if ttype:
            if itype==list:
                if ttype=="double":
                    tvalue = []
                    for x in value:
                        tvalue.append(float(x))
                if ttype=="long":
                    tvalue = []
                    for x in value:
                        tvalue.append(int(x))
            else:
                if ttype=="double":
                    tvalue = float(value)
                if ttype=="long":
                    tvalue = int(value)
        itype = type(tvalue)
        # List
        if itype==list:
            itype = type(tvalue[0])  # Type of array element
            dim[0] = len(tvalue)
            dim[1] = 1
            # float as double
            if itype==float:
               if ttype=="double":
                    Obit.InfoListAlwaysPutDouble(self.me, name, dim, tvalue)
               else:
                   Obit.InfoListAlwaysPutFloat(self.me, name, dim, tvalue)
               return
               # int
            elif itype==int:
                if ttype=="long":
                    Obit.InfoListAlwaysPutLong(self.me, name, dim, tvalue)
                else:
                    Obit.InfoListAlwaysPutInt(self.me, name, dim, tvalue)
                return
                # string
            elif itype==str:
                dim[0] = len(tvalue[0])
                for x in tvalue[0:]:
                    dim[0] = max(dim[0], len(x))
                dim[1] = len(value)
                if (dim[1]>1):    # array of strings
                    Obit.InfoListAlwaysPutString(self.me, name, dim, tvalue)
                else:             # single
                    Obit.InfoListAlwaysPutSString(self.me, name, dim, tvalue)
                return
                # boolean
            elif itype==bool:
                Obit.InfoListAlwaysPutBoolean(self.me, name, dim, tvalue)
                PAlwaysPutBoolean (self, name, dim, value)
                return
            # List
            if itype==list:
                dim[1] = len(tvalue)
                dim[0] = len(tvalue[0])
                temp = []
                for x in tvalue:
                    for y in x:
                        temp.append(y)
                jtype = type(tvalue[0][0])  # Type of array of arrays
                # float as double
                if jtype==float:
                    if ttype=="double":
                        Obit.InfoListAlwaysPutDouble(self.me, name, dim, temp)
                    else:
                        Obit.InfoListAlwaysPutFloat(self.me, name, dim, temp)
                    return
                # int
                elif jtype==int:
                    if ttype=="long":
                        Obit.InfoListAlwaysPutLong(self.me, name, dim, temp)
                    else:
                        Obit.InfoListAlwaysPutInt(self.me, name, dim, temp)
                    return
                # string
                elif jtype==str:
                    dim[2] = len(tvalue[0][0])
                    for x in value[0][0:]:
                        dim[2] = max(dim[2], len(x))
                    if (dim[2]>1):    # array of strings
                        Obit.InfoListAlwaysPutString(self.me, name, dim, temp)
                    else:             # single
                        Obit.InfoListAlwaysPutSString(self.me, name, dim, temp)
                    return
                # boolean
                elif jtype==bool:
                    Obit.InfoListAlwaysPutBoolean(self.me, name, dim, temp)
                    PAlwaysPutBoolean (self, name, dim, value)
                    return
                else:
                    raise RuntimeError,"Unknown type "+str(jtype)

            else:
                raise RuntimeError,"Unknown type "+str(itype)

            # Scalars
            # float 
        elif itype==float:
            if ttype=="double":
                Obit.InfoListAlwaysPutDouble(self.me, name, dim, [tvalue])
            else:
                Obit.InfoListAlwaysPutFloat(self.me, name, dim, [tvalue])
            return
            # int 
        elif itype==int:
            if ttype=="long":
                Obit.InfoListAlwaysPutLong(self.me, name, dim, [tvalue])
            else:
                Obit.InfoListAlwaysPutInt(self.me, name, dim, [tvalue])
            return
            # string
        elif itype==str:
            dim[0] = len(value)
            Obit.InfoListAlwaysPutSString(self.me, name, dim, [tvalue])
            return
            # boolean
        elif itype==bool:
            Obit.InfoListAlwaysPutBoolean(self.me, name, dim, [tvalue])
            return
        else :
            raise RuntimeError,"Unsupported type",itype
        return
    # end set
        
    def get (self, name):
        """ Retrieve a value from an InfoList
  
        returns python list containing data:
        0 - return code, 0=OK else failed
        1 - name
        2 - type
            int=1, oint=3, long=4, float=9, double=10, string=13, boolean=14
        3 - dimension array as list, e.g. [1,1,1,1,1] for scalar
        4 - data array
        self     = input Python InfoList
        name     = name of desired entry
        """
        # Ignore any exception when item not found
        try:
            out = Obit.InfoListGet(self.me, name)
            return out
        except:
            return out
        else:
            return out
    # end set
    # end InfoList Class
    
# Commonly used, dangerous variables
dim=[1,1,1,1,1]

def PCopy (inList):
    """ Copy list

    return copy of input InfoList
    inList    = input Python InfoList
    """
    ################################################################
    out = InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.InfoListCopy(inList.me)
    return out
# end PCopy 

def PCopyData (inList, outList):
    """ Copy all entries from inList to outList

    inList    = input Python InfoList
    outList   = output Python InfoList, previously exists
    """
    ################################################################
    Obit.InfoListCopyData(inList.me, outList.me)
    # end PCopyData

def PRemove (inList, name):
    """ Removes item name fro list

    inList    = input Python InfoList
    name      = name of desired entry
    """
    ################################################################
    Obit.InfoListRemove (inList.me, name)
    # end PRemove

def PItemResize (inList, name, type, dim):
    """ 

    inList   = input Python InfoList
    name     = name of desired entry
    type     = data type of object
               int=1, oint=3, long=4, float=9, double=10, string=13, boolean=14
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
    """
    ################################################################
    Obit.InfoListItemResize(inList, name, type, dim)
    # end PItemResize

def PIsA (inList):
    """ Tells if input really is InfoList

    returns true, false (1,0)
    inList    = input Python InfoList
    """
    ################################################################
     # Checks
    if inList.__class__ != InfoList:
        return 0
    return Obit.InfoListIsA(inList.me)
    # end PIsA

def PPutInt (inList, name, dim, data, err):
    """ Add an integer entry, error if conflict

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of integers
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if data[0].__class__ != int:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be int"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListPutInt(inList.me, name, dim, data, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error adding entry")
    # end PPutInt

def PAlwaysPutInt (inList, name, dim, data):
    """ Add an integer entry, changing type/dim of entry if needed

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of integers
    """
    ################################################################
    # Checks
    if data[0].__class__  != int:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be int"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListAlwaysPutInt(inList.me, name, dim, data)
    # end PAlwaysPutInt

def PPutLong (inList, name, dim, data, err):
    """ Add an long entry, error if conflict

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of integers
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if data[0].__class__ != int:
        raise TypeError,"data MUST be int"
    if len(dim) < 5:
        print "class is" , data[0].__class__
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListPutLong(inList.me, name, dim, data, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error adding entry")
    # end PPutLong

def PAlwaysPutLong (inList, name, dim, data):
    """ Add an long entry, changing type/dim of entry if needed

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of integers
    """
    ################################################################
    # Checks
    if data[0].__class__ != int:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be int"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListAlwaysPutLong(inList.me, name, dim, data)
    # end PAlwaysPutLong

def PPutFloat (inList, name, dim, data, err):
    """ Add an float entry, error if conflict

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of float
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if data[0].__class__ != float:
        raise TypeError,"data MUST be float"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
        print "data, size",prod,data
        raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListPutFloat(inList.me, name, dim, data, err.me)
    if err.isErr:
       OErr.printErrMsg(err, "Error adding entry")
    # end PPutFloat

def PAlwaysPutFloat (inList, name, dim, data):
    """ Add an float entry, changing type/dim of entry if needed

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of float
    """
    ################################################################
    # Checks
    if data[0].__class__ != float:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be float"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListAlwaysPutFloat(inList.me, name, dim, data)
    # end PAlwaysPutFloat

def PPutDouble (inList, name, dim, data, err):
    """ Add an double entry, error if conflict

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of double
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if data[0].__class__ != float:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be float"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListPutDouble(inList.me, name, dim, data, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error adding entry")
    # end  PPutDouble

def PAlwaysPutDouble (inList, name, dim, data):
    """ Add an integer entry, changing type/dim of entry if needed

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of double
    """
    ################################################################
    # Checks
    if data[0].__class__ != float:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be float"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListAlwaysPutDouble(inList.me, name, dim, data)
    # end PAlwaysPutDouble

def PPutBoolean (inList, name, dim, data, err):
    """ Add an boolean entry (1,0), error if conflict

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of boolean (1,0)
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if data[0].__class__ != bool:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be bool"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListPutBoolean(inList.me, name, dim, data, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error adding entry")
    # end PPutBoolean

def PAlwaysPutBoolean (inList, name, dim, data):
    """ Add an boolean entry, changing type/dim of entry if needed

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of boolean (1,0)
    """
    ################################################################
    # Checks
    if data[0].__class__ != bool:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be bool"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
        raise RuntimeError,"more data than defined in dim"
    if prod > len(data):
        raise RuntimeError,"less data than defined in dim"
    Obit.InfoListAlwaysPutBoolean(inList.me, name, dim, data)
    # end PAlwaysPutBoolean

def PPutString (inList, name, dim, data, err):
    """ Add an string entry, error if conflict

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of strings (rectangular char array)
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if data[0].__class__ != str:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be a string"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
         raise RuntimeError,"more data than defined in dim"
    # Call depends on whether single string, 1 or 2 D arrays
    if (dim[1]>1):    # array of strings
        Obit.InfoListPutString(inList.me, name, dim, data, err.me)
    else:             # single
        Obit.InfoListPutSString(inList.me, name, dim, data, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error adding entry")
    # end  PPutString

def PAlwaysPutString (inList, name, dim, data):
    """ Add an String entry, changing type/dim of entry if needed

    inList   = input Python InfoList
    name     = name of desired entry
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
               MUST have 5 entries
    data     = data as a 1-D array of (rectangular char array)
    """
    ################################################################
    # Checks
    if data[0].__class__ != str:
        print "class is" , data[0].__class__
        raise TypeError,"data MUST be a string"
    if len(dim) < 5:
        raise RuntimeError,"dim has fewer then 5 entries"
    prod = 1
    for x in dim:
        if (x>0):
            prod = prod*x
    if prod < len(data):
        raise RuntimeError,"more data than defined in dim"
    # Call depends on whether single string, 1 or 2 D arrays
    if (dim[1]>1):    # array of strings
        Obit.InfoListAlwaysPutString(inList.me, name, dim, data)
    else:             # single
        Obit.InfoListAlwaysPutSString(inList.me, name, dim, data)
     # end PAlwaysPutString

def PGet (inList, name):
    """ 
    returns python list containing data:
       0 - return code, 0=OK else failed
       1 - name
       2 - type
       3 -  dimension array
       4 - data array
    inList   = input Python InfoList
    type     = data type of object
               int=1, oint=3, long=4, float=9, double=10, string=13, boolean=14
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
    data     = data as a multi dimensional array of (rectangular char array)
    """
    ################################################################
    # Checks
    if not PIsA(inList):
        print "Actually ",inList.__class__
        raise TypeError,'inList MUST be a Python Obit InfoList'
    
    retVal = Obit.InfoListGet(inList.me, name)
    if retVal[0] != 0:
        raise RuntimeError,name+' not found in InfoList'
    return retVal
    # end PGet

def PGetDict (inList):
    """ 
    returns contents as python dict, each element is:
       0 - type
       1 - dimension array
       2 - data array
    type     = data type of object
               int=2, oint=3, long=4, float=9, double=10, string=13, boolean=14
    dim      = dimensionality of array as list, e.g. [1,1,1,1,1] for scalar
    data     = data as a multi dimensional array of (rectangular char array)
    inList    = input Python InfoList
    """
    ################################################################
    # Checks
    if not PIsA(inList):
        print "Actually ",inList.__class__
        raise TypeError,'inList MUST be a Python Obit InfoList'
    
    return Obit.InfoListGetDict(inList.me)
    # end PGetDict

def PSetDict (outList, inDict):
    """ 
    Adds entries in InfoList from elements of a python dict, each element is
    a list with the following elements:
       [0] - type code:int=2, oint=3, long=4, float=9, double=10, string=13, boolean=14
       [1] - dim, dimension array, an list of 5 int
       [2] - data list, a 1-D list all of the data type given by type

    outList  = output Python InfoList
    inDict   = input dict strictly in form generated by PGetDict
    """
    ################################################################
    # Checks
    if not PIsA(outList):
        print "Actually ",outList.__class__
        raise TypeError,'outList MUST be a Python Obit InfoList'
    if not type(inDict)==dict:
        raise TypeError,'inDict MUST be a Python dict'
    # Loop over dict
    for x in inDict:
        data = inDict[x]
        #print "DEBUG",x,data
        # Check that entry is an list
        if not type(data)==list:
            print "Bad entry:",data
            raise TypeError,"entry for "+x+" MUST be a list"
        # Check type (data[0])
        if not type(data[0])==int:
            raise TypeError,"type for "+x+" MUST be an int"
        # Check dim (data[1])
        if not type(data[1])==list:
            raise TypeError,"dim for "+x+" MUST be a list of 5 int"
        if not type(data[1][0])==int:
            raise TypeError,"dim for "+x+" MUST be a list of 5 int"
        # data (data[2]) must be as a list
        if not type(data[2])==list:
            raise TypeError,"data for "+x+" MUST be a list"
        #
        # Switch on type
        if data[0] in [1,2,3,4]:       # int
            PAlwaysPutInt (outList, x, data[1], data[2])
        elif data[0]==9:             # float
            PAlwaysPutFloat (outList, x, data[1], data[2])
        elif data[0]==10:            # double
            PAlwaysPutDouble (outList, x, data[1], data[2])
        elif data[0]==13:            # string
            PAlwaysPutString (outList, x, data[1], data[2])
        elif data[0]==14:            # boolean
            PAlwaysPutBoolean (outList, x, data[1], data[2])
        else:  # Unknown
            raise RuntimeError,"Unknown data type code "+str(data[0])+" for "+x
    # end loop over entries
    # end PSetDict

def PIsA (inList):
    """ Tells if input really a Python Obit InfoList

    return true, false (1,0)
    inList   = Python InfoList object
    """
    ################################################################
     # Checks
    if inList.__class__ != InfoList:
        return 0
    return Obit.InfoListIsA(inList.me)
    # end PIsA

def PUnref (inList):
    """ Decrement reference count

    Decrement reference count which will destroy object if it goes to zero
    Python object stays defined.
    inList   = Python InfoList object
    """
    ################################################################
     # Checks
    if not PIsA(inList):
        raise TypeError,"inList MUST be a Python Obit InfoList"

    inList.me = Obit.InfoListUnref(inList.me)
    # end PUnref

