# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

'''
This package contains basic classes that are used throughout Quantas for::

    - reading input file;
    - set up calculations.

Each class used in Quantas for one of these purposes inherits from the classes
here present.

'''

from quantas.core.calculator import BasicCalculator
from quantas.core.reader import BasicReader

__all__ = ['BasicCalculator', 'BasicReader']
