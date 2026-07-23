# -*- coding: utf-8 -*-

"""Compatibility facade for volume-temperature EOS fitting adapters.

New internal code should import from :mod:`quantas.modules.eos.domains.vt`.
The historical module path remains stable for third-party callers.
"""

from .domains import vt as _domain
from .domains.vt import *  # noqa: F403

__all__ = _domain.__all__
