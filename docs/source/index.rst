.. Align-BDD documentation master file, created by
   sphinx-quickstart on Tue Jun 22 23:36:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:
   :maxdepth: 3

   self


Code documentation for Align-BDD project.
=========================================
These pages document the implementation of the numerical experiments for "Align-BDD" working paper by Bochkarev and Smith. These notes along with the accompanying source code and data are assumed to provide enough information to (1) inspect the raw data, (2) reproduce (and adjust, if nececssary) the figures from the paper, and (3) completely reproduce the paper results from scratch, with minimal efforts.

Questions can be addressed to `Alexey Bochkarev <https://www.bochkarev.io/contact>`_.

Problem instances and raw data.
---------------------------------
All the source data is available in the public repository. In particular,

* 📁 `instances` containts all the instance ("raw") data in a compressed format (`tar.gz`),
* 📁 `run_logs` contains all the numerical results ("run logs" as CSV files)

Reproducibility
---------------
Assuming the necessary software is present (see the next section), one can rebuild any figure from the existing log file (specific experiment results) by running a ``make`` command from the root project directory (``align-BDD``). For example,

.. code:: bash

   make figures/simpl_heuristics.eps  # to rebuild Figure 4, or
   make figures  # to rebuild all figures for the paper

GNU make will consider files' timestamps and re-run the experiments as necessary (e.g., if the instance generation procedure has changed). To force a complete run (generate all the data, make all experiments, and build all figures), just clean up the previous results first:

.. code:: bash

   make fresh  # to delete all instances and run logs
   make figures  # to rebuild everything

Note that depending on the parameters (which are set in ``Makefile``) and the hardware, it might take quite some time to recalculate everything. PBS scripts [#]_ that were used to run our experiments on the computation cluster are kept in ``📁 pbs`` .

All the recipies are implemented in a `Makefile` (that's the default makefile name), which is basically a list of commands and source files needed to build each target.  

.. [#] obviously, they are dependent on specific software and hardware on the computation system you use.

Computational infrastructure
----------------------------
To run the experiments we used machines running GNU/Linux system (laptops and remote computational resources provided by `Clemson Palmetto cluster <https://www.palmetto.clemson.edu/palmetto/about/>`_), running the following software: 

* experiments implementation: python3 (:download:`python-packages.txt`),
* figures: R (:download:`R-libraries.txt`),
* parallel computation: `GNU parallel <https://www.gnu.org/software/parallel/>`_,
* experiments pipeline: GNU make and PBS (for job scheduling on the cluster).
* (we also used a separate `tool <https://github.com/alex-bochkarev/tgs-curl>`_ to receive the results promptly via `Telegram <https://telegram.org>`_, but this can be safely removed)
* numerical experiments require `Gurobi <https://www.gurobi.com>`_ solver.

Code organization
-----------------
The code is organized into several blocks.

Simple instructions to that physically produce figures from CSV run logs are written in `R <https://www.r-project.org/>`_ (and the wonderful tidyverse/`ggplot2 <https://ggplot2.tidyverse.org>`_ library). The code is conceptually simple and, hopefully, self-explanatory -- see ``.R`` - files in ``📁 post_processing``.

The core data structures critical to the paper are presented in the following modules (click on the links for more implementation details):

.. autosummary:: 
   :toctree: _autosummary
   :recursive:

   BDD
   varseq
   heuristics
   BB_search

Description for specific problems serves the basis for further numerical experiments, which are implemented as separate python programs in ``📁 experiments``

.. autosummary::
   :toctree: _autosummary
   :recursive:

   gen_BDD_pair
   cUFL
   jUFL
   experiments

Key implementations are covered with tests using `pytest <https://pytest.org>`_ framework (the tests are in ``📁 tests``).

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
