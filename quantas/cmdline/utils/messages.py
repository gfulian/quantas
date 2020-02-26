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
import click
import quantas

FGCOLORS = {
    'success': 'green',
    'error': 'bright_red',
    'warning': 'bright_yellow',
    'highlight': 'bright_white',
    }

def quantas_title():
    """ Return the Quantas software title. """
    msg = r"""
________                       __                 
\_____  \  __ _______    _____/  |______    ______
 /  / \  \|  |  \__  \  /    \   __\__  \  /  ___/
/   \_/.  \  |  // __ \|   |  \  |  / __ \_\___ \ 
\_____\ \_/____/(____  /___|  /__| (____  /____  >
       \__>          \/     \/          \/     \/ """
    msg += '\n'
    msg += '{0:>45}{1}\n'.format('v',quantas.__version__)
    msg += '{0}: {1}\n'.format('Authors', quantas.__author__)
    msg += '{0}\n'.format(quantas.__copyright__)
    return msg


def quantas_error():
    """ Return a big error message. """
    msg = r"""
_______________________________ ________ __________._.
\_   _____/\______   \______   \\_____  \\______   \ |
 |    __)_  |       _/|       _/ /   |   \|       _/ |
 |        \ |    |   \|    |   \/    |    \    |   \\|
/_______  / |____|_  /|____|_  /\_______  /____|_  /__
        \/         \/        \/         \/       \/ \/"""
    msg += '\n'
    return msg


def quantas_warning():
    """ Return a big warning message. """
    msg = r"""
 __      __                     .__                
/  \    /  \_____ _______  ____ |__| ____    ____  
\   \/\/   /\__  \\_  __ \/    \|  |/    \  / ___\ 
 \        /  / __ \|  | \/   |  \  |   |  \/ /_/  >
  \__/\  /  (____  /__|  |___|  /__|___|  /\___  / 
       \/        \/           \/        \//_____/"""
    msg += '\n'
    return msg


def quantas_finish():
    """ Return a big software closing message. """
    msg = r"""
___________.__       .__       .__     
\_   _____/|__| ____ |__| _____|  |__  
 |    __)  |  |/    \|  |/  ___/  |  \ 
 |     \   |  |   |  \  |\___ \|   Y  \
 \___  /   |__|___|  /__/____  >___|  /
     \/            \/        \/     \/"""
    msg += '\nThank you for using this software!\n'
    return msg


def init_logfile(logfile=None):
    """
    Init the file used to log information on the current software run.

    Parameters
    ----------

    logfile: str
        Path to the log file.
    
    """
    if not isinstance(logfile, type(None)):
        print('', file=open(logfile, 'w'))
    return


def echo(msg, logfile=None, silent=False, bold=False):
    """
    Print a message via the click.echo function to stdout and also to a log
    file, it its path is provided.

    Parameters
    ----------

    msg: str
        Message to be printed.

    logfile: str, optional
        Path to the log file.

    silent: bool, optional
        If False, the message will not be printed to stdout (output strem).

    bold: bool, optional
        Whether to print the message using a bold style.

    """
    if not silent:
        click.secho(msg, bold=bold)
    if not isinstance(logfile, type(None)):
        print(msg, file=open(logfile, 'a'))
    return


def echo_error(msg, logfile=None, silent=False, bold=False):
    """
    Print an error message via the click.echo function to stdout and also to a
    log file, it its path is provided.

    Parameters
    ----------

    msg: str
        Message to be printed.

    logfile: str, optional
        Path to the log file.

    silent: bool, optional
        If False, the message will not be printed to stdout (output strem).

    bold: bool, optional
        Whether to print the message using a bold style.

    """
    click.secho(msg, bold=bold, fg=FGCOLORS['error'])
    if not isinstance(logfile, type(None)):
        print(msg, file=open(logfile, 'a'))
    return

def echo_warning(msg, logfile=None, silent=False, bold=False):
    """
    Print a warning message via the click.echo function to stdout and also to
    a log file, it its path is provided.

    Parameters
    ----------

    msg: str
        Message to be printed.

    logfile: str, optional
        Path to the log file.

    silent: bool, optional
        If False, the message will not be printed to stdout (output strem).

    bold: bool, optional
        Whether to print the message using a bold style.

    """
    click.secho(msg, bold=bold, fg=FGCOLORS['warning'])
    if not isinstance(logfile, type(None)):
        print(msg, file=open(logfile, 'a'))
    return


def echo_highlight(msg, logfile=None, silent=False, bold=True):
    """
    Print an highlighted message via the click.echo function to stdout and
    also to a log file, it its path is provided.

    Parameters
    ----------

    msg: str
        Message to be printed.

    logfile: str, optional
        Path to the log file.

    silent: bool, optional
        If False, the message will not be printed to stdout (output strem).

    bold: bool, optional
        Whether to print the message using a bold style.

    """
    if not silent:
        click.secho(msg, bold=bold, fg=FGCOLORS['highlight'])
    if not isinstance(logfile, type(None)):
        print(msg, file=open(logfile, 'a'))
    return


def confirm(question, default=False, suffix=': ', show_default=True):
    """
    Print a confirm message via the click.confirm function to stdout and
    prompt for a user answer.

    Parameters
    ----------

    question: str
        Question message to be printed.

    default: bool, optional
        The default answer that can be prompted by simply giving Return.

    suffix: str, optional
        The suffix appended to the message, before the user prompt.

    show_default: bool, optional
        Whether to show the default values in the prompt.

    Returns
    -------

    bool
        Boolean answer to the yes/no question
    """
    return click.confirm(
        click.style(question, fg=FGCOLORS['highlight']),
        default=default, prompt_suffix=suffix, show_default=show_default)


def prompt(question, default=None, type=None, suffix=': ', show_default=True):
    """
    Print a message via the click.prompt function to stdout and prompt for an
    answer from the user.

    Parameters
    ----------

    question: str
        Question message to be printed.

    default: bool, optional
        The default answer that can be prompted by simply giving Return.

    suffix: str, optional
        The suffix appended to the message, before the user prompt.

    show_default: bool, optional
        Whether to show the default values in the prompt.

    """
    return click.prompt(
        click.style(question, fg=FGCOLORS['highlight']), type=type,
        default=default, prompt_suffix=suffix, show_default=show_default)
