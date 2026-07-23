# -*- coding: utf-8 -*-

"""Compatibility facade for pressure-volume EOS fitting adapters.

New internal code should import from :mod:`quantas.modules.eos.domains.pv`.
The historical module path remains stable for third-party callers.
"""

from .domains import pv as _domain
from .domains.pv import *  # noqa: F403

__all__ = _domain.__all__
