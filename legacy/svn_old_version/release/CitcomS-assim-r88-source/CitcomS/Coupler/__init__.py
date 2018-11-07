#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#<LicenseText>
#
# CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
# Copyright (C) 2002-2005, California Institute of Technology.
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
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#</LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


def containingcoupler(name, facility):
    from ContainingCoupler import ContainingCoupler
    return ContainingCoupler(name, facility)


def embeddedcoupler(name, facility):
    from EmbeddedCoupler import EmbeddedCoupler
    return EmbeddedCoupler(name, facility)


def multicontainingcoupler(name, facility):
    from MultiC_Coupler import MultiC_Coupler
    return MultiC_Coupler(name, facility)


def multiembeddedcoupler(name, facility):
    from MultiE_Coupler import MultiE_Coupler
    return MultiE_Coupler(name, facility)


# version
__id__ = "$Id: __init__.py 7713 2007-07-19 18:37:13Z hlin $"

# End of file
