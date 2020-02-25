# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

from quantas.cmdline.utils.messages import echo
from quantas.cmdline.utils.messages import echo_error
from quantas.cmdline.utils.messages import confirm
from quantas.cmdline.utils.messages import prompt

import numpy as np
from .pv_eos import Murnaghan
from .pv_eos import BirchMurnaghan
from .pv_eos import NaturalStrain
from .pv_eos import Vinet
from .pv_eos import Tait

EOS_FORMULATIONS = {
    '1': Murnaghan(),
    '2': BirchMurnaghan(),
    '3': NaturalStrain(),
    '4': Vinet(),
    '5': Tait(),
    }

EOS_FUNCTIONS = {
    '1': Murnaghan,
    '2': BirchMurnaghan,
    '3': NaturalStrain,
    '4': Vinet,
    '5': Tait,
    }

EOS_ORDERS = {
    '1': [3],
    '2': [2, 3, 4],
    '3': [2, 3, 4],
    '4': [2, 3],
    '5': [2, 3, 4],
    }

DATA_TYPES = {
    'v': 'Volume',
    'a': 'a-axis',
    'b': 'b-axis',
    'c': 'c-axis',
    }

VOLUME_LABELS = ['K0', "K'", "K''", 'V0']

AXIS_LABELS = ['M0({})', "M'({})", "M''({})", '{}0']

PLOT_LABES = {
    'v' : '$Volume (' + u'\u212B^3' + ')$',
    'a' : '$a lattice parameter (' + u'\u212B' + ')$',
    'b' : '$b lattice parameter (' + u'\u212B' + ')$',
    'c' : '$c lattice parameter (' + u'\u212B' + ')$'
    }


def is_available(data):
    """
    Check if the data is available for fitting procedure.

    Parameters
    ----------

    data: ndarray
        Array containing the data.

    Returns
    -------

    bool
        If the provided data is an array, return True, otherwhise False
    """
    if isinstance(data, np.ndarray):
        return True
    else:
        return False


def show_choices(choices):
    """
    Print the available choices for a specific selection by the user.

    Paramters
    ---------

    choices: dict or list
        Dictionary or list containing the available choices. If a dictionary,
        its keys will be used as bullet for the printed list.

    """
    if isinstance(choices, dict):
        for key, item in choices.items():
            echo('{: >3}. {}'.format(key, item))
    elif isinstance(choices, list):
        for choice in choices:
            echo(' - {}'.format(choice))
    echo('')
    return


def choose(msg, choices, default=None):
    """
    Prompt a choice to the user.

    Parameters
    ----------

    msg: str
        The question asked to the user.

    choices: dict or list
        List of the available choices to pick from.

    default: str, optional
        Default choice.

    """
    show_choices(choices)
    while True:
        selection = prompt(msg, default=default)
        if selection in choices:
            return selection
        else:
            echo_error('Error: invalid selection')


def set_eos(default='1'):
    """
    Ask the user to pick one of the available equation of state formulation.

    Parameters
    ----------

    default: str, optional
        The default EoS formulation.

    Returns
    -------

    selection: str
        The key for the selected EoS formulation stored in the
        EOS_FORMULATIONS dictionary.
        
    """
    echo('EOS FORMULATION')
    selection = choose('Select an EoS formulation', EOS_FORMULATIONS,
                       default=default)
    echo('')
    return selection


def set_eos_order(eos_key, default=2):
    """
    Ask the user to select the order of the picked equation of state
    formulation. The orders can be 2, 3 or 4, depending on the chosen EoS.

    Parameters
    ----------

    eos_key: str
        The key of the selected EoS Formulation.

    default: int, optional
        The default order of the EoS formulation.

    Returns
    -------

    selection: int
        Order of the EoS formulation.
        
    """
    choices = EOS_ORDERS[eos_key]
    if default not in choices:
        default = choices[0]
    echo('EOS ORDER')
    selection = np.int(choose('Select the order of the EoS', EOS_ORDERS[eos_key],
                       default=default))
    echo('')
    return selection


def select_eos(defaults=None):
    """
    Ask the user to select both equation of state formulation and its order.

    Parameters
    ----------

    default: dict, optional
        Dictionary containing default choices that will be prompted as default
        answers.

    Returns
    -------

    eos: cls
        Class of the selected equation of state.

    eos_order: int
        Order ot the selected equation of state.

    defaults: dict
        Updated (or new) dictionary containing the current EoS selections.

    """
    if isinstance(defaults, type(None)):
        defaults={'formulation': '1', 'order': 2}
    eos_sel = set_eos(defaults['formulation'])
    eos_order = set_eos_order(eos_sel, defaults['order'])
    defaults['formulation'] = eos_sel
    defaults['order'] = eos_order
    return EOS_FUNCTIONS[eos_sel](), eos_order, defaults


def select_data(data):
    """
    Ask the user to select among the available data provided as input.

    Parameters
    ----------

    data: dict
        Dictionary containing the input data.

    Returns
    -------

    selection: str
        Key of the selected data.

    data_name: str
        Name of the selected data.

    selected_data: ndarray
        Array of the selected data.

    linearized: bool
        Whether to treat the data as linear (for cell axes) or not
        (for volume).

    """
    available = {}
    for key in DATA_TYPES.keys():
        if is_available(data[key]):
            available[key] = DATA_TYPES[key]

    echo('DATA SELECTION')
    selection = choose('Select the values to fit', available, None)
    if selection != 'v':
        linearized = True
    else:
        linearized = False
    return selection, DATA_TYPES[selection], data[selection], linearized


def get_label(q):
    """
    Get the label of the EoS parameters, according to the selected data.

    Parameters
    ----------

    q: str
        Key of the selected data.

    Returns
    -------

    list
        List of the labels of the EoS parameters.

    """
    if q == 'v':
        return VOLUME_LABELS
    else:
        return [ label.format(q) for label in AXIS_LABELS ]


def set_parameters(parameters, labels):
    """
    Ask the user if any of the EoS parameters has to be modified. If yes,
    prompt the new values. Default values are those of the unmodified
    parameters.

    Parameters
    ----------

    parameters: ndarray
        Array of the EoS parameter values.

    labels: list
        List of the labels of the EoS parameters.

    Returns
    -------

    ndarray
        Array of the updated EoS parameters.

    """
    echo('EOS PARAMETERS')
    for i in range(len(parameters)):
        echo('  {: <8} = {: 12.5f}'.format(labels[i], parameters[i]))
    echo('')
    if confirm('Would you like to modify any of the EoS parameters?'):
        for i in range(len(parameters)):
            parameters[i] = np.float(prompt("  {: <8}".format(
                labels[i]), np.round(parameters[i], 5)))
    echo('')
    return parameters


def set_fixed_parameters(labels):
    """
    Ask the user if any of the EoS parameters has to be fixed during the
    fitting procedure.

    Parameters
    ----------

    labels: list
        List of the labels of the EoS parameters.

    Returns
    -------

    list
        list containing 0 (for fixed parameters) or 1 (for free parameters)
        values.

    """
    if confirm('Would you like to fix any of the EoS parameters?'):
        return [0 if confirm('  Fix {: <6}'.format(label)) else 1
                for label in labels]
    else:
        return [1, 1, 1, 1]


def set_weights(q, data):
    """
    Ask the user if weighting scheme has to be used for one or both the
    pressure and selected data values. The question will be asked only if
    at least one of pressure or selected data have experimental uncertainties.

    Parameters
    ----------

    q: str
        Key of the selected data.

    data: dict
        Dictionary of the input data.

    Returns
    -------

    list
        List containing either None (for unweighted data) or the experimental
        uncertainties related to the data.

    """
    weights = [None, None]
    if not is_available(data['sigmap']) and not is_available(data['sigma'+q]):
        return weights

    if confirm('Would you like to weight the data during fitting?',
               default=True):
        if is_available(data['sigmap']):
            if confirm(' Use weights for {: <8}?'.format('pressure'),
                       default=True):
                weights[0] = data['sigmap'].copy()

        if is_available(data['sigma'+q]):
            if confirm(' Use weights for {: <8}?'.format(DATA_TYPES[q]),
                       default=True):
                weights[1] = data['sigma'+q].copy()
    return weights
