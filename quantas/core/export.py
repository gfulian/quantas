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
This module contains a basic export class that can be used as a base for
other more sophisticated ones.
"""


class BasicHDF5Export(object):
    """
    A barebone export class.

    Attributes
    ----------

    error: None or str
        Error encountered during data handling and storing.

    _data: dict
        Dictionary containing the data collected from the input.

        .. warning::

            If multiple instances of BasicReader (or any class that inherits
            from it) are used at the same time, the _data attribute **must**
            be initialized in the constructor method. Otherwise, all the
            instances will contain the *same* data (which is usually not
            desired).

    Methods
    -------

    load(filename)
        Read the file and store the data contained in it.

    _check()
        Check for errors in the stored data.
    """

    completed = False
    error = None
    _data = None

    def __init__(self, input_file=None):
        """
        Constructor for the basic reader.

        Each class that inherits from quantascli.core.calculator.BasicReader
        must employ the same construction method, namely the class immediately
        reads an input file, if provided during instantiation.

        Parameters
        ----------

        input_file: str, optional
            Path to the input file. If not None, the file is immediately read.

        """
        if not isinstance(input_file, type(None)):
            self.load(input_file)
        return

    def load(self, filename):
        """
        Dummy method employed to read the input file and store data. It could
        be called during class instantiation if an input file is provided.

        If the input file is entirely read, the BasicReader.completed flag
        should be set to True.

        If any error is encountered in parsing the file (e.g. missing data),
        the BasicReader.error attribute should be set to a string value that
        describes the error source.

        Nothing should be returned, everything must be kept inside this class.

        Parameters
        ----------

        filename: str
            Path to the input file

        """
        return

    def run(self):
        """
        Dummy method that starts an interactive procedure to export the
        collected data.
        """
        return
