# -*- coding: utf-8 -*-

"""Compatibility import for shared grouped Click options.

The implementation is generic and lives in :mod:`quantas.cli.grouped_options`.
This historical EOS import remains stable for internal and external callers.
"""

from quantas.cli.grouped_options import GroupedCommand, GroupedOption, grouped_option

__all__ = ["GroupedCommand", "GroupedOption", "grouped_option"]
