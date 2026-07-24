Numerical precision and tolerances
==================================

Quantas calculations and native HDF5 values use ``float64`` for real quantities
and ``complex128`` for complex quantities. Display precision belongs to
renderers and does not modify stored data.

.. admonition:: Work in progress

   This page will explain tolerance selection, conditioning, uncertainty-aware
   comparisons, and the distinction between regression tolerances and physical
   acceptance criteria.
