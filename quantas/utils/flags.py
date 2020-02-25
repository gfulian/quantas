# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

"""
This module contains classes that can be used as flags.
"""

class Flag(object):
    """
    Generic flag object.

    Parameters
    ----------

    activated: bool, optional
        When the object is instantiated, its status is set according to this
        value.

    """

    def __init__(self, activated=False):
        if activated:
            self.on()
        else:
            self.off()
        return

    def on(self):
        """
        Activate the flag.
        """
        self.__status = True
        return

    def off(self):
        """
        Deactivate the flag.
        """
        self.__status = False
        return

    def is_on(self):
        """
        Return the current status of the flag with  `bool` type.
        """
        return self.__status
