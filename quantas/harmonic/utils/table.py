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
subhd = '{s}{v:{c}<{w}}{s}'                     # Sub-Header format
t_space = 10                                    # Space for T values
t_line = '{s}{t:{c}^{w}}'                       # T format
d_space = 20                                    # Space for data values
d_line = '{s}{d:{c}^{w}}'                       # Data format
endrow = '{s:{c}>{w}}'                          # End of row format
ttype = '%.2f'                                  # T numeric format
etype = '%.12E'                                 # Energy numeric format
vtype = '%.8f'                                  # Volume numeric format
ptype = '%.8f'                                  # Pressure numeric format
stype = '%.s'                                   # Text format
nl = '\n'                                       # New line
fileHeader = '*- This file was created using Quantas -*\n'


# ---------------------------------------------------------------------------#
# Writing functions
# ---------------------------------------------------------------------------#
def calc_space(n):
    totw = (n * d_space) + t_space + (n + 1)
    dataw = (n * d_space) + (n - 1)
    return totw, dataw


def separator(totw, sym, char):
    return sep.format(s=sym, w=totw, c=char)


def header(dataw, tag1, tag2):
    return head.format(s=' ', v1=tag1, w1=t_space, v2=tag2, w2=dataw, c=' ')


def subhead(totw, tag):
    return subhd.format(s=' ', v=tag, w=totw-1, c='')


def put_data(x, dtype):
    if x == '':
        value = t_line.format(s=' ', t='', w=t_space, c=' ')
    elif dtype == ttype:
        value = t_line.format(s=' ', t=str(dtype % x), w=t_space, c=' ')
    elif dtype == stype:
        value = d_line.format(s=' ', d=(x), w=d_space, c=' ')
    else:
        value = d_line.format(s=' ', d=str(dtype % x), w=d_space, c=' ')
    return value


def newrow():
    return endrow.format(s=' \n', w=0, c=' ')


# ---------------------------------------------------------------------------#
# Writing classes
# ---------------------------------------------------------------------------#
class QHATableFile():
    '''
    File export class for (Q)HA results. It exports data in a table format
    in selected file.

    Parameters
    ----------

    outfile = str
        name of the file where data will be saved.

    x = ndarray
        First independent variable data.

    xname = str
        Name of the variable of var1.

    xunit = str
        Measurement unit of var1.

    y = ndarray
        Second independent variable data.

    yname = str
        Name of the variable of var2.

    yunit = str
        Measurement unit of var1.

    z = ndarray
        Results in a 2D-array format (1D or 2D).

    zname = str
        Name of the result data.

    zunit = str
        Measurement unit of result data.

    '''
    types = {
        'T': {'desc': 'Temperature', 'label': 'T ({})', 'numstyle': ttype},
        'V': {'desc': 'Volume', 'label': 'V ({}^3)', 'numstyle': vtype},
        'U0': {'desc': 'Static energy', 'label': 'E0 ({})', 'numstyle': etype},
        'Uzp': {'desc': 'Zero-point energy', 'label': 'Ezp ({})',
                'numstyle': etype},
        'Uth': {'desc': 'Thermal energy', 'label': 'Eth ({})',
                'numstyle': etype},
        'Utot': {'desc': 'Internal energy', 'label': 'E ({})',
                 'numstyle': etype},
        'S': {'desc': 'Entropy', 'label': 'S ({})', 'numstyle': etype},
        'Cv': {'desc': 'Isochoric heat capacity', 'label': 'Cv ({})',
               'numstyle': etype},
        'Fvib': {'desc': 'Vibrational Helmholtz free energy',
                 'label': 'Fvib ({})', 'numstyle': etype},
        'F': {'desc': 'Helmholtz free energy', 'label': 'F ({})',
              'numstyle': etype},
        }
    def __init__(self):
        return

    def write(self, outfile,
              x, xname, xunit,
              y, yname, yunit,
              z, zname, zunit
              ):
        with open(outfile, 'wt', encoding='utf8') as f:
            totw, dataw = calc_space(len(y))
            # Write the header
            f.write(fileHeader)
            f.write(nl)
            f.write(self.types[zname]['desc'])
            f.write(nl)
            f.write('Data in ' + zunit + ' units')
            f.write(nl)
            f.write(nl)
            # Start constructing the table
            f.write(separator(totw, '+', '-'))
            f.write(nl)
            f.write(header(dataw, self.types[xname]['label'].format(xunit),
                           self.types[yname]['label'].format(yunit)))
            f.write(nl)
            f.write(put_data('', self.types[xname]['numstyle']))
            for i in range(len(y)):
                f.write(put_data(y[i], self.types[yname]['numstyle']))
            f.write(newrow())
            f.write(separator(totw, '+', '='))
            f.write(nl)
            for i in range(len(x)):
                f.write(put_data(x[i], self.types[xname]['numstyle']))
                for j in range(len(y)):
                    f.write(put_data(z[i][j], self.types[zname]['numstyle']))
                f.write(newrow())
            f.write(separator(totw, '+', '='))
        return
    
