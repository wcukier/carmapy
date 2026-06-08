Citing CARMApy
==============

The exact citation for CARMApy will depend on which parts of it you are using.
At a minimum, please cite both CARMApy itself and the underlying models that it
is built on. For example::

   This work uses CARMApy (Cukier et al. in prep), which is built on the
   ExoCARMA (Gao et al. 2018, Powell et al. 2018) version of the CARMA model
   (Turco et al. 1979, Toon et al. 1988, Bardeen et al. 2008).

If you use additional features of CARMApy, such as the 2D mode or the FastChem
chemistry interface, please also cite the relevant papers listed below.

Base Citation
-------------

CARMApy
~~~~~~~

.. code-block:: bibtex

   @article{CukierinprepCARMApyOpensourcePython,
     title  = {{CARMApy}: A Open-Source Python Wrapper for the Microphysical Cloud Code {CARMA}},
     author = {Cukier, Wolf and Powell, Diana and Zhang, Xi and Gao, Peter and Samra, Dominic and Nagpal, Vighnesh},
     year   = {in prep}
   }

ExoCARMA
~~~~~~~~

CARMApy wraps ExoCARMA, the version of CARMA adapted for exoplanetary atmospheres
by Gao et al. (2018) and Powell et al. (2018):

.. code-block:: bibtex

   @article{Gao2018SedimentationEfficiencyCondensation,
     title        = {Sedimentation Efficiency of Condensation Clouds in Substellar Atmospheres},
     author       = {Gao, Peter and Marley, Mark S. and Ackerman, Andrew S.},
     date         = {2018-03-10},
     journaltitle = {The Astrophysical Journal},
     shortjournal = {ApJ},
     volume       = {855},
     number       = {2},
     pages        = {86},
     doi          = {10.3847/1538-4357/aab0a1}
   }

   @article{Powell2018FormationSilicateTitanium,
     title        = {Formation of Silicate and Titanium Clouds on Hot Jupiters},
     author       = {Powell, Diana and Zhang, Xi and Gao, Peter and Parmentier, Vivien},
     date         = {2018-06},
     journaltitle = {The Astrophysical Journal},
     shortjournal = {ApJ},
     volume       = {860},
     number       = {1},
     pages        = {18},
     doi          = {10.3847/1538-4357/aac215}
   }

CARMA
~~~~~

The underlying CARMA model was originally developed in the works below:

.. code-block:: bibtex

   @article{Turco1979OneDimensionalModelDescribing,
     title        = {A One-Dimensional Model Describing Aerosol Formation and Evolution in the Stratosphere: I. Physical Processes and Mathematical Analogs},
     author       = {Turco, R. P. and Hamill, P. and Toon, O. B. and Whitten, R. C. and Kiang, C. S.},
     date         = {1979-04-01},
     journaltitle = {Journal of the Atmospheric Sciences},
     volume       = {36},
     pages        = {699--717},
     doi          = {10.1175/1520-0469(1979)036<0699:AODMDA>2.0.CO;2}
   }

   @article{Toon1988MultidimensionalModelAerosols,
     title        = {A Multidimensional Model for Aerosols: Description of Computational Analogs},
     author       = {Toon, O. B. and Turco, R. P. and Westphal, D. and Malone, R. and Liu, M.},
     date         = {1988-08-01},
     journaltitle = {Journal of the Atmospheric Sciences},
     volume       = {45},
     number       = {15},
     pages        = {2123--2143},
     doi          = {10.1175/1520-0469(1988)045<2123:AMMFAD>2.0.CO;2}
   }

   @article{Bardeen2008NumericalSimulationsThreedimensional,
     title        = {Numerical Simulations of the Three-Dimensional Distribution of Meteoric Dust in the Mesosphere and Upper Stratosphere},
     author       = {Bardeen, C. G. and Toon, O. B. and Jensen, E. J. and Marsh, D. R. and Harvey, V. L.},
     date         = {2008},
     journaltitle = {Journal of Geophysical Research: Atmospheres},
     volume       = {113},
     number       = {D17},
     doi          = {10.1029/2007JD009515}
   }

2D CARMApy
----------

If you use the two-dimensional mode of CARMApy, please also cite the model it is
based on, Powell & Zhang (2024):

.. code-block:: bibtex

   @article{Powell2024TwodimensionalModelsMicrophysical,
     title        = {Two-Dimensional Models of Microphysical Clouds on Hot Jupiters. I. Cloud Properties},
     author       = {Powell, Diana and Zhang, Xi},
     date         = {2024-07-01},
     journaltitle = {The Astrophysical Journal},
     volume       = {969},
     pages        = {5},
     doi          = {10.3847/1538-4357/ad3de4}
   }

FastChem
--------

If you use the FastChem chemistry interface, please cite Stock et al. (2022):

.. code-block:: bibtex

   @article{Stock2022FASTCHEM2Improved,
     title        = {FastChem 2: An Improved Computer Program to Determine the Gas-Phase Chemical Equilibrium Composition for Arbitrary Element Distributions},
     author       = {Stock, Joachim W. and Kitzmann, Daniel and Patzer, A. Beate C.},
     date         = {2022-12-01},
     journaltitle = {Monthly Notices of the Royal Astronomical Society},
     volume       = {517},
     pages        = {4070--4080},
     doi          = {10.1093/mnras/stac2623}
   }
