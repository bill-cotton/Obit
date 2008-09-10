"""
This module provides the Task class.  It extends the MinimalMatch
class from the MinimalMatch module with type and range checking on its
attributes:

>>> class MyTask(Task):
...     indisk = 0
...     inseq  = 0
...     infile = ''
...     pixavg = 1.0
...     aparms = 10*[0.0]
...     def __init__(self):
...         Task.__init__(self)
...         self._min_dict = {'inseq': 0, 'aparms': 0}
...         self._max_dict = {'inseq': 4, 'aparms': 10}
...         self._strlen_dict = {'infile': 14}
...         self.__dict__['bparms'] = List(self, 'bparms', [None, 1, 2, 3])
...
>>> my_task = MyTask()

It still has the property that attribute names can be abbreviated:

>>> print my_task.ind
0
>>> my_task.ind = 1
>>> print my_task.ind
1

But an exception will be thrown if you try to assign a value that is
out of range:

>>> my_task.ins = 5
Traceback (most recent call last):
  ...
ValueError: value '5' is out of range for attribute 'inseq'

Or if you try to assign a value that has the wrong type, such
assigning a string to an integer attribute:

>>> my_task.ind = 'now'
Traceback (most recent call last):
  ...
TypeError: value 'now' has invalid type for attribute 'indisk'

Assigning strings to string attributes works fine of course:

>>> my_task.infile = 'short'

As long as there is no limit on the length of a string:

>>> my_task.infile = 'tremendouslylong'
Traceback (most recent call last):
  ...
ValueError: string 'tremendouslylong' is too long for attribute 'infile'

Assigning an integer value to a floating point attribute is perfectly
fine of course:

>>> my_task.pixavg = 2
>>> print my_task.pixavg
2.0

The same should happen for lists:

>>> my_task.aparms = 10*[1]
>>> print my_task.aparms
[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

For subscripting:

>>> my_task.aparms[0] = 0
>>> print my_task.aparms
[0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

And slice assignment:

>>> my_task.aparms[1:3] = [1, 2]
>>> print my_task.aparms
[0.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

You're not allowed to change the length of the list through slice
assignment though:

>>> my_task.aparms[3:6] = [3, 4, 5, 6]
Traceback (most recent call last):
  ...
TypeError: slice '3:6' changes the array size of attribute 'aparms'

To provide 1-based indexing used by several packages, you can set the
element at index zero of an array to 'None'.  This prevents setting that
element to anything other than 'None'

>>> my_task.bparms[0] = 0
Traceback (most recent call last):
  ...
ValueError: setting element '0' is prohibited

"""
# Copyright (C) 2005 Joint Institute for VLBI in Europe
# Copyright (C) 2006,2007 Associated Universities, Inc. Washington DC, USA.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from MinimalMatch import MinimalMatch

# Generic Python stuff
import pydoc, sys

class List(list):
    """ List class to convert lists to POPS/AIPS form """
    def __init__(self, task, attr, value):
        self._task = task
        self._attr = attr
        _value = []
        for item in value:
            if isinstance(item, list):
                _value.append(List(task, attr, item))
            else:
                _value.append(item)
                pass
            continue
        list.extend(self, _value)
        return
    
    def __setitem__(self, key, item):
        if item != None and self[key] == None:
            msg = "setting element '%d' is prohibited" % key
            raise ValueError, msg
        item = self._task._validateattr(self._attr, item, self[key])
        list.__setitem__(self, key, item)
        return

    def __setslice__(self, low, high, seq):
        high = min (high, len(self))
        if len(seq) > high - low or \
               (len(seq) < high - low and high < len(self)):
            msg = "slice '%d:%d' changes the array size of" \
                  " attribute '%s'" % (low, high, self._attr)
            raise TypeError, msg
        for key in xrange(low, high):
            if key - low < len(seq):
                self[key] = seq[key - low]
            else:
                self[key] = self._task._default_dict[self._attr][key]
                pass
            continue
        return


class Task(MinimalMatch):
    """ Basic (virtual) client side Task class """
    def __init__(self):
        """  Create AIPS task object
        
        Creates task object 
        Following is a list of class members:
        _default_dict   = Dictionary with default values of parameters
        _input_list     = List of input parameters in order
        _output_list    = List of output parameters in order
        _min_dict       = Parameter minimum values as a List
        _max_dict       = Parameter maximum values as a List
        _hlp_dict       = Parameter descriptions (list of strings)
        as a dictionary
        _strlen_dict    = String parameter lengths as dictionary
        _help_string    = Task Help documentation as list of strings
        _explain_string = Task Explain documentation as list of strings
        _short_help     = One line description of task
        """
        self._default_dict = {}
        self._min_dict = {}
        self._max_dict = {}
        self._strlen_dict = {}
        self._help_dict = {}
        self._short_help = ''
        self._help_string = ''
        self._explain_string = ''
        self._remainder = ""   # Partial message buffer

    def help(self):
        """Display help for this task."""

        if self._help_string:
            pydoc.pager(self._help_string)

    def explain(self):
        """Display help+explain for this task."""

        if self._help_string and self._explain_string:
            pydoc.pager(self._help_string+self._explain_string)
        else:
            print "No explanation available - see help"

    def _validateattr(self, attr, value, default):
        """Check whether VALUE is a valid valid for attribute ATTR."""

        # Do not check private attributes.
        if attr.startswith('_'):
            return value

        # Short circuit.
        if value == None and default == None:
            return value

        # Handle lists recursively.
        if isinstance(value, list) and isinstance(default, list):
            if len(value) > len(default):
                msg = "array '%s' is too big for attribute '%s'" \
                      % (value, attr)
                raise TypeError, msg
            validated_value = List(self, attr, default)
            for key in xrange(len(value)):
                validated_value[key] = value[key]
            return validated_value

        # Convert integers into floating point numbers if necessary.
        if type(value) == int and type(default) == float:
            value = float(value)

        # Check attribute type.
        if type(value) != type(default):
            msg = "value '%s' has invalid type for attribute '%s'" \
                  % (value, attr)
            raise TypeError, msg

        # Check range.
        if attr in self._min_dict:
            min = self._min_dict[attr]
            if not min <= value:
                msg = "value '%s' is out of range for attribute '%s'" \
                      % (value, attr)
                raise ValueError, msg
        if attr in self._max_dict:
            max = self._max_dict[attr]
            if not value <= max:
                msg = "value '%s' is out of range for attribute '%s'" \
                      % (value, attr)
                raise ValueError, msg

        # Check string length.
        if attr in self._strlen_dict:
            if len(value) > self._strlen_dict[attr]:
                msg = "string '%s' is too long for attribute '%s'" \
                      % (value, attr)
                raise ValueError, msg

        return value

    def __setattr__(self, name, value):
        attr = self._findattr(name)

        # Validate based on the value already present.
        if hasattr(self, attr):
            # Don't check proxies
            if attr != "proxy":
                value = self._validateattr(attr, value, getattr(self, attr))
        self.__dict__[attr] = value

    def __getattr__(self, name):
        """ Allow access to functions without parentheses """
        if name=="input":
            return self.inputs()
        if name=="i":
            return self.inputs()
        if name=="o":
            return self.outputs()
        if name=="h":
            return self.help()
        if name=="e":
            return self.explain()
        if name=="g":
            return self.go()
        if name=="a":
            return self.abort()
        if name=="w":
            return self.wait()

    def parseMessage(self, inMess):
        """
        Repackage messages to prevent line breaks
        
        Given the message buffer including line break characters, return
        an array of strings containing full lines (originally terminated
        by line break character).  Partial lines are saved in variable
        remainder and are prepended to the next input message buffer
        """
        outMess=[]
        #print "inMess=",inMess,"\n\n"
        for msg in inMess:
            prio = msg[0]
            mess = self._remainder+msg[1]
            #print "mess=",mess,"\n\n"
            self._remainder=""
            while len(mess)>0:
                # Look for newline or carrage return, either or both may be
                # end of line marker
                lf = mess.find("\n")
                cr = mess.find("\r")
                es = max (lf, cr)
                # Find one?
                if es<0:
                    # No - keep rest in remainder
                    self._remainder = mess
                    mess = ""
                    break
                # Copy line (minus end of line characters) to outMess
                if cr>0 and cr<lf:
                    es = cr
                omess = mess[0:es].strip("\n\r")
                outMess.append((prio,omess))
                # Drop line through end of line character(s)
                es = max (cr, lf) + 1
                es = max (es,1)
                mess = mess[es:]
                #print "save msg",mess[0:es],"\n\n"
       # End loop over entries in inMess
                
        #print "outMess=",outMess,"\n\n"
        return outMess
    # end parseMessage
# End class Task

# Tests.
if __name__ == '__main__':
    import doctest, sys
    results = doctest.testmod(sys.modules[__name__])
    sys.exit(results[0])


