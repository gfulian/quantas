# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

# ---------------------------------------------------------------------------#
# Formatting constants
# ---------------------------------------------------------------------------#
sep = '{s:{c}<{w}}{s}'                          # Row separator format
head = '{s}{v1:{c}^{w1}}{s}{v2:{c}^{w2}}{s}'    # Header format
dw = 10                                         # Space for degree values
vw = 15                                         # Space for data values
nl = '\n'                                       # New line
fileHeader = '*- This file was created using Quantas -*\n'

# ---------------------------------------------------------------------------#
# Writing functions
# ---------------------------------------------------------------------------#
def calc_space(n):
    totw = (n * vw) + dw
    dataw = (n * vw)
    return totw, dataw


def separator(totw, sym, char):
    return sep.format(s=sym, w=totw, c=char)


def header(dataw, tag1, tag2):
    return head.format(s=' ', v1=tag1, w1=dw, v2=tag2, w2=dataw, c=' ')


def newrow():
    return endrow.format(s=' \n', w=0, c=' ')


# ---------------------------------------------------------------------------#
# Writing class
# ---------------------------------------------------------------------------#
class SOECTableFile():
    """
    File export class for SOEC results. It exports data in a table format
    in selected file.

    Parameters
    ----------

    outfile = str
        name of the file where data will be saved.

    degrees = ndarray
        First independent variable data.

    values = ndarray
        Results in a 2D-array format (1D or 2D).

    name = str
        Name of the result data.

    unit = str
        Measurement unit of result data.

    """
    types = {
        'E': {'desc': "Young's modulus", 'label': {
            '(xy)': [''],
            '(xz)': [''],
            '(yz)': [''],
            }},
        'LC': {'desc': 'Linear compressibility', 'label': {
            '(xy)': ['positive', 'negative'],
            '(xz)': ['positive', 'negative'],
            '(yz)': ['positive', 'negative'],
            }},
        'G': {'desc': 'Shear modulus', 'label': {
            '(xy)': ['minimum', 'maximum'],
            '(xz)': ['minimum', 'maximum'],
            '(yz)': ['minimum', 'maximum'],
            }},
        'Nu': {'desc': "Poisson's ratio", 'label': {
            '(xy)': ['negative', 'positive', 'measurement'],
            '(xz)': ['negative', 'positive', 'measurement'],
            '(yz)': ['negative', 'positive', 'measurement'],
            }},
        'waves': {'desc': 'Seismic wave velocities', 'label': {
            '(xy)': ['v_s1', 'vs_2', 'v_l'],
            '(xz)': ['v_s1', 'vs_2', 'v_l'],
            '(yz)': ['v_s1', 'vs_2', 'v_l'],
            }},
        }

    def __init__(self):
        return

    def write(self, outfile, degrees, value, name, unit):
        with open(outfile, 'wt', encoding='utf8') as f:
            totw, dataw = calc_space(value.shape[1])
            # Write the header
            f.write(fileHeader)
            f.write(nl)
            f.write(self.types[name]['desc'])
            f.write(nl)
            if unit == '':
                f.write('Adimensional data')
            else:
                f.write('Data in ' + unit + ' units')
            f.write(nl)
            f.write(nl)
            # Start constructing the table
            f.write(separator(totw, '+', '-'))
            f.write(nl)
            f.write(header(dataw, 'Degrees', self.types[name]['desc']))
            f.write(nl)
            subdir = '{: ^{d}}'.format('', d=dw)
            subval = '{: ^{d}}'.format('', d=dw)
            for direction in self.types[name]['label']:
                subdir += '{: ^{w}}'.format(
                    direction, w=int(value.shape[1]/3)*vw)
                for quantity in self.types[name]['label'][direction]:
                    subval += '{: ^{w}}'.format(quantity, w=vw)
            f.write(subdir)
            f.write(nl)
            f.write(subval)
            f.write(nl)
            f.write(separator(totw, '+', '='))
            f.write(nl)
            for i in range(value.shape[0]):
                line = '{: ^{d}.1f}'.format(degrees[i], d=dw)
                for j in range(value.shape[1]):
                    line += '{: ^{w}.4f}'.format(value[i, j], w=vw)
                f.write(line)
                f.write(nl)
            f.write(separator(totw, '+', '='))
        return
    
