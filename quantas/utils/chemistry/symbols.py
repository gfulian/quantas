# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

from quantas.utils.chemistry import chemical_symbols, atomic_numbers

def number2symbol(Z):
    return chemical_symbols[Z]

def symbol2number(s):
    return atomic_numbers[s]
