"""Implements most of the numerical experiments in the paper.

Almost every script here implements a stand-alone experiment, which is run from
``Makefile`` with parameters passed via the command line from the ``make``
recipe. Passing ``--help`` to a script usually prints the available options.

The only notable exception is :py:mod:`UFLP_2_cav`, describing the j-UFLP
experiment, which is implemented as a stand-alone module.

(Note that the module contains some experiments that were left out from the
second revision of the paper.)
"""
