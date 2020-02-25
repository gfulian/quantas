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
This module contains a basic reader classes that can be used as a base for
other more sophisticated ones.
"""


class BasicReader(object):
    """
    A barebone reader class.

    Attributes
    ----------

    completed: bool
        Flag that is True when the input file was completely read, without
        errors, otherwise it should be False.

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

    def _check(self):
        """
        Dummy, 'internal use' method employed to check if any possible error
        is present in the collected data.

        For example, you could check if the stored data in attribute 'x' is
        a float number, or if it is an array.

        If an error is encountered, a small description of the error should
        be put in the BasicReader.error attribute.

        Nothing should be returned.

        """
        return

    @property
    def data(self):
        """
        Return the information collected from the input file.

        For each information stored in BasicReader._data, a property that has
        direct (read-only) access to that information should be created.
        """
        return self._data
