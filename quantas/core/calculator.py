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
This module contains a basic class that can be used as a base for other more
sophisticated classes.
"""

from quantas.cmdline.utils.messages import echo
from quantas.cmdline.utils.messages import echo_warning
from quantas.cmdline.utils.messages import echo_error


class BasicCalculator(object):
    """
    A barebone calculator class.

    Attributes
    ----------

    completed: bool
        Flag that is True when the input file was completely read, without
        errors, otherwise it should be False.

    error: None or str
        Error encountered during data handling and storing.

    warnings: list
        List containing warnings encountered during the calculations, which
        did not cause an interruption of the workflow.

    descriptions: dict
        Dictionary containing descriptions of each data calculated by the
        class.

    _results: dict
        Dictionary containing the data collected from the input.

        .. warning::

            If multiple instances of BasicCalculator (or any class that
            inherits from it) are used at the same time, the _data attribute
            **must** be initialized in the constructor method.
            Otherwise, all the instances will contain the *same* data (which
            is usually not desired).

    Methods
    -------

    read_input(filename)
        Read the file and store the data contained in it. In Quantas, the
        input file is loaded by a reader class that inherits from
        quantascli.core.reader.BasicReader class.

    report_input()
        Report the input data collected.

    run()
        Start the expected calculations.

    report_results()
        Report the calculated data.

    set_array(n, m=None)
        Set a 1D- or 2D-array that is used during calculations.

    """

    completed = False

    error = None

    warnings = []

    debug = False

    log = None

    silent = False
    
    _results = {
        'var1': None,
        'var2': None,
        }

    descriptions = {
        'var1': 'var1 description',
        'var2': 'var2 description',
        }

    def echo(self, msg):
        echo(msg, self.log, silent=self.silent)
        return

    def echo_override(self, msg):
        echo(msg, self.log, silent=False)
        return

    def echo_warning(self, msg):
        echo_warning(msg, self.log, silent=self.silent)
        return

    def echo_debug(self, msg):
        if self.debug:
            echo('Debug: ' + msg, self.log, silent=self.silent)
        return

    def echo_time(self, t0, t1):
        elapse = t1 - t0
        echo('   Operation time {0:15.3f} msec\n'.format(elapse * 1000),
             self.log, silent=self.silent)
        return

    def __init__(self):
        """ Constructor for the basic calculator. """
        return

    def read_input(self, filename):
        """
        Dummy method that will be overwritten in other classes.

        When creating a class that inherits from
        quantascli.core.calculator.BasicCalculator, this method accepts
        an input file and must return None if no error was encountered during
        the reading process, otherwise a string containing the source of
        error in input file should be provided.

        Parameters
        ----------
        filename: str
            Path to the input file.

        Returns
        -------
        error: None or str
            If a string, it represents the encountered error when file was
            read.

        """
        return None

    def report_input_data(self):
        """
        Dummy method that will be overwritten in other classes.

        When creating a class that inherits from
        quantascli.core.calculator.BasicCalculator, no argument will be passed
        to this method and nothing should be returned.

        """
        return

    def run(self):
        """
        Dummy method that will be overwritten in other classes.

        When creating a class that inherits from
        quantascli.core.calculator.BasicCalculator, no argument will be passed
        to this method and nothing should be returned.

        The results of the calculations should be stored inside the class.
        """
        return

    def report_results(self):
        """
        Dummy method that could be overwritten in other classes.

        When creating a class that inherits from
        quantascli.core.calculator.BasicCalculator, no argument will be passed
        to this method and nothing should be returned.

        """
        return

    def export(self, filename):
        """
        Dummy method that could be overwritted in other classes.

        In Quantas, the output data are usually stored in HDF5 file format
        (binary).

        """
        return

    @property
    def results(self):
        """
        Return the calculated results with dict format.
        """
        return self._results
