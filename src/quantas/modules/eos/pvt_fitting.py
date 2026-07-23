# -*- coding: utf-8 -*-

"""Compatibility facade for pressure-volume-temperature fit adapters.

New internal code should import from :mod:`quantas.modules.eos.domains.pvt`.
The historical module path remains stable for third-party callers.
"""

from .domains import pvt as _domain
from .domains.pvt import *  # noqa: F403

__all__ = _domain.__all__
