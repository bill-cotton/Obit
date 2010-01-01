# Copyright (C) 2005 Joint Institute for VLBI in Europe
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

"""
Modify rlcompleter.Completer class such that it does not show hidden
attributes.
"""

import readline
import rlcompleter

class ObitTalkCompleter(rlcompleter.Completer):
    def attr_matches(self, text):
        words = rlcompleter.Completer.attr_matches(self, text)
        new_words = []
        for word in words:
            pos = word.rfind('.')
            if pos > 0:
                if word.find('._', pos) == pos:
                    continue;
            new_words.append(word)
        return new_words

readline.set_completer(ObitTalkCompleter().complete)
