# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import h5py

class QuantasHDF5Reader(object):

    def __init__(self, filename):
        self._settings = {}
        self._settings['input'] = filename
        return

    @property
    def settings(self):
        return self._settings

    def read(self):
        with h5py.File(self.settings['input'],'r') as f:
            self._info = f.attrs['info']
            self._info += '\nAvailable data:\n'
            for key in f.keys():
                self._info += '{: >15}: '.format(key)
                for attribute in f[key].attrs.keys():
                    if attribute == 'unit':
                        self._info += '({})'.format(f[key].attrs[attribute])
                    else: self._info += f[key].attrs[attribute] + ' '
                self._info += '\n'
        return

    @property
    def info(self):
        return self._info

    def get_data(self, key):
        if self.__has_dataset(key) != True:
            return self.__has_dataset(key)
        else:
            with h5py.File(self.settings['input'],'r') as f:
                return f[key][()]

    def has_data(self, key):
        return self.__has_dataset(key)

    def get_data_column(self, key, column):
        if self.__has_dataset(key) != True:
            return self.__has_dataset(key)
        else:
            try:
                with h5py.File(self.settings['input'],'r') as f:
                    return f[key][:,column]
            except ValueError:
                return 'Index {0} not available'.format(column)

    def get_data_row(self, key, row):
        if self.__has_dataset(key) != True:
            return self.__has_dataset(key)
        else:
            try:
                with h5py.File(self.settings['input'],'r') as f:
                    return f[key][row]
            except ValueError:
                return 'Index {0} not available'.format(row)

    def get_attribute(self, key, attribute):
        if not self.__has_dataset(key):
            return self.__has_dataset(key)
        elif not self.__has_attribute(key, attribute):
            return self.__has_attribute(key, attribute)
        else:
            with h5py.File(self.settings['input'], 'r') as f:
                return f[key].attrs[attribute]

    def __has_dataset(self, key):
        with h5py.File(self.settings['input'],'r') as f:
            if key not in f.keys():
                return "Dataset '{0}' not found".format(key)
            else:
                return True

    def __has_attribute(self,key,attribute):
        with h5py.File(self.settings['input'], 'r') as f:
            if attribute not in f[key].attrs.keys():
                return "Dataset '{0}' does not have '{1}' attribute".format(
                    key,attribute)
            else:
                return True
