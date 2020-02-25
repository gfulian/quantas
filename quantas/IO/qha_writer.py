# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import os
import numpy as np
from copy import deepcopy
from quantas.interfaces.crystal import CrystalPhononReader
from quantas.interfaces.crystal import CrystalQHAReader
from quantas.interfaces.phonopy import PhonopyReader


class QHAInputCreator(object):
    """
    This class creates input files for the (Q)HA framework implemented in
    Quantas.
    """
    
    interface_filter = {
        'crystal': CrystalPhononReader,
        'crystal-qha': CrystalQHAReader,
        'phonopy': PhonopyReader,
        }
    
    interface = None
    phondata = []
    reference = {}

    def __init__(self, interface=None):
        """
        """
        self.interface_flag = interface
        return

    def read(self, filename, is_list, reference=0, use_symm=False):
        """
        """
        if self.interface_flag == 'crystal-qha' and is_list:
            error = 'Only one CRYSTAL output file for --crystal-qha interface'
            return False, error

        c = 0

        if is_list:
            file_list = self.get_files_from_list(filename)
        else:
            file_list = [filename]
            
        for file in file_list:

            if not os.path.exists(file):
                error = 'File {} does not exists'.format(file)
                return False, error
            print('  - {} ...'.format(file), end='')

            interface = self.interface_filter[self.interface_flag]
            data = interface(file)

            if not data.completed:
                return False, data.error

            if self.interface_flag == 'crystal-qha':
                self.phondata = data
            else:
                self.phondata.append(data)

            if c == reference and self.interface_flag != 'crystal-qha':
                print(' Done! (Reference)')
            else:
                print(' Done!')

            c += 1
            
        if use_symm and self.interface_flag != 'crystal-qha':
            print()
            self.phonon_symmetry(reference, file_list[reference])
        return True, None

    def get_files_from_list(self, file_list):
        with open(file_list, "r") as f:
            files = f.readlines()
        for i in range(len(files)):
            files[i] = files[i].rstrip()
        return files

    def phonon_symmetry(self, reference=0, fname=''):
        """
        """
        refdata = self.phondata[reference]
        tmp_total_phonons = refdata.phonons.copy()
        reduced_phonons = {}
        reduced_qcoord = {}
        reduced_weights = {}
        c = 0
        print('Starting the reduction of phonon bands using the reference')
        print('file: {}'.format(fname))
        for key in refdata.phonons:
            print_progress(key, refdata.qpoints)
            k_phonons = refdata.phonons[key]
            check = False  # Flag to check if the phonon band was found
            for phon in reduced_phonons:
                if np.array_equal(k_phonons, reduced_phonons[phon]):
                    tmp_total_phonons.pop(key, None)
                    reduced_weights[phon] += 1
                    if not check:
                        check = True

            if not check:
                reduced_phonons[c] = refdata.phonons[key]
                reduced_weights[c] = 1
                reduced_qcoord[c] = refdata.qcoords[key] 
                tmp_total_phonons.pop(key, None)
                c += 1
        print()                
        print(' Done!')
        print()
        #
        # Store the reduced data in the reference
        self.phondata[reference].data['qcoords'] = reduced_qcoord
        self.phondata[reference].data['phonons'] = reduced_phonons
        self.phondata[reference].data['weights'] = reduced_weights
        self.phondata[reference].data['qpoints'] = len(reduced_qcoord)
        #
        # Search the same q-points in the other files
        if len(self.phondata) > 1:
            print()
            print('Apply changes to other data')
            for i in range(len(self.phondata)):
                if i == reference:
                    print(' - Skipping reference data')
                    print()
                    continue
                print(' - Processing data #{}'.format(i))
                tmp_coord = {}
                tmp_phons = {}
                tmp_weight = {}
                for key, ref_qcoord in refdata.qcoords.items():
                    print_progress(key, refdata.qpoints)
                    for q, item in self.phondata[i].qcoords.items():
                        if np.array_equal(ref_qcoord, item):
                            tmp_coord[key] = item
                            tmp_phons[key] = self.phondata[i].phonons[q]
                            tmp_weight[key] = refdata.weights[key]
                self.phondata[i].data['qcoords'] = tmp_coord
                self.phondata[i].data['phonons'] = tmp_phons
                self.phondata[i].data['weights'] = tmp_weight
                self.phondata[i].data['qpoints'] = len(refdata.qcoords)
                print()
                print()
            print('All done!')
        return

    def write(self, outfile, ref=0):
        """
        """
        msg = '\nPlease, enter a short description for the input file: '
        jobname = input(msg)

        if self.interface_flag == 'crystal-qha':
            self.write_from_single_file(outfile, jobname, ref)
        else:
            self.write_from_multiple_files(outfile, jobname, ref)
        return

    def write_from_multiple_files(self, outfile, jobname, ref=0):
        n = len(self.phondata)
        reference = self.phondata[ref]
        data = []
        
        data.append('job: {}'.format(jobname))
        data.append('natom:   %-7d' % reference.natom)
        data.append('supercell:')
        for i in range(3):
            vect_line = '- [ ' 
            for j in range(3):
                if j != 2:
                    vect_line += '%6d, '
                else:
                    vect_line += '%6d ]'
            data.append(vect_line % tuple(reference.dim[i]))
        data.append('qpoints: %i' % reference.qpoints)
        vline = 'volume: [ '
        eline = 'energy: [ '
        for i in range(n):
            if i != n - 1:
                vline += '{:16.8f},'.format(self.phondata[i].volume)
                eline += '{: 20.12E},'.format(self.phondata[i].energy)
            else:
                vline += '{:16.8f} ]'.format(self.phondata[i].volume)
                eline += '{: 20.12E} ]'.format(self.phondata[i].energy)
        data.append(vline)
        data.append(eline)
        data.append('phonon:')
        for i in range(reference.qpoints):
            data.append('- q-position: [ %12.7f, %12.7f, %12.7f ]'
                        % tuple(reference.qcoords[i]/reference.shrinkf))
            data.append('  weight: %-5d' % reference.weights[i])
            data.append('  band:')
            for j in range(reference.nphonon):
                data.append('  - # {0}'.format(j+1))
                band = []
                for item in self.phondata:
                    band.append(item.phonons[i][j])
                bline = '    frequency: [ '
                for k in range(n):
                    if k != n - 1:
                        bline += '{:16.10f}, '.format(band[k])
                    else:
                        bline += '{:16.10f} ]'.format(band[k])
                data.append(bline)
            data.append('')

        with open(outfile, 'w') as f:
            f.write('\n'.join(data))
        return

    def write_from_single_file(self, outfile, jobname, ref=0):
        n = self.phondata.points
        data = []
        data.append('job: {}'.format(jobname))
        data.append('natom:   %-7d' % self.phondata.natom)
        data.append('supercell:')
        for i in range(3):
            vect_line = '- [ ' 
            for j in range(3):
                if j != 2:
                    vect_line += '%6d, '
                else:
                    vect_line += '%6d ]'
            data.append(vect_line % tuple(self.phondata.dim[i]))
        data.append('qpoints: %i' % self.phondata.qpoints)
        vline = 'volume: [ '
        eline = 'energy: [ '
        for i in range(n):
            if i != n - 1:
                vline += '{:16.8f},'
                eline += '{: 20.12E},'
            else:
                vline += '{:16.8f} ]'
                eline += '{: 20.12E} ]'
        data.append(vline.format(*self.phondata.volume))
        data.append(eline.format(*self.phondata.energy))
        data.append('phonon:')
        for i in range(self.phondata.qpoints):
            data.append(
                '- q-position: [ %12.7f, %12.7f, %12.7f ]'
                % tuple(self.phondata.qcoords[i]/self.phondata.shrinkf))
            data.append('  weight: %-5d' % self.phondata.weights[i])
            data.append('  band:')
            for j in range(self.phondata.nphonon):
                data.append('  - # {0}'.format(j+1))
                band = self.phondata.phonons[i][j]
                bline = '    frequency: [ '
                for k in range(n):
                    if k != n - 1:
                        bline += '{:16.10f}, '
                    else:
                        bline += '{:16.10f} ]'
                data.append(bline.format(*band))
            data.append('')

        with open(outfile, 'w') as f:
            f.write('\n'.join(data))
        return


def print_progress(iteration, total, prefix = '', suffix = '',
                   decimals = 1, length = 50, fill = '█', printEnd = '\r'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. '\r', '\r\n') (Str)
    """
    percent = ('{0:.' + str(decimals) + 'f}').format(
        100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix),
          end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
