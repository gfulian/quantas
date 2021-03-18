# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

import click
import difflib

from .messages import quantas_title
from .messages import quantas_error
from .messages import echo_highlight
from .messages import echo_error


class MostSimilarCommandGroup(click.Group):
    """
    Overloads the get_command to display a list of possible command
    candidates if the command could not be found with an exact match.
    """

    def get_command(self, ctx, cmd_name):
        """
        Override the default click.Group get_command with one giving the
        user a selection of possible commands if the exact command name
        could not be found.
        """
        cmd = click.Group.get_command(self, ctx, cmd_name)

        if cmd is not None:
            return cmd

        matches = difflib.get_close_matches(
            cmd_name, self.list_commands(ctx), cutoff=0.5)

        if not matches:
            matches = [c for c in sorted(self.list_commands(ctx))
                       if c.startswith(cmd_name)][:3]

        if matches:
            echo_error(quantas_error())
            ctx.fail(
                "'{}' is not a Quantas command.\n\n"
                'The most similar commands are: \n'
                '{matches}'.format(
                    cmd_name,
                    matches='\n'.join('\t{}'.format(m) for m in sorted(matches))
                    )
            )
        else:
            echo_error(quantas_error())
            ctx.fail("'{}' is not a Quantas command.\n\n"
                     'No similar commands found.'.format(cmd_name))

        return None


def print_version(ctx, param, value):
    """ Check the context and print Quantas version. """
    if not value or ctx.resilient_parsing:
        return
    echo_highlight(quantas_title(), bold=True)
    ctx.exit()
