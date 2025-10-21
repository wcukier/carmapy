Installation
===============

Dependencies
------------
Most dependencies will be automatically installed for you but in order to 
install carmapy you must first either have the intel fortran compiler (preferred)
or gfortran installed.  To check that they are installed you can run the 
following commands::

   which ifort
   which gfortran

as long as at least one of those commands returns as file path rather than saying
"ifort not found" or "gfortran not found" then you should be good to proceed with
installation.  If you are installing on an apple silicon device, you may be able
to get away with installing using pip even if you don't have a fortran compiler

Installation with pip
---------------------
To install with pip simply write::
   
   pip install carmapy

Installation with git
---------------------
To install and build loacally::

   git clone --recursive https://github.com/wcukier/carmapy
   cd carmapy
   pip install .


Testing Installation
--------------------
To test that carmapy is installed correctly you can run the following python code

.. code-block:: python

   import carmapy.example

   carma = carmapy.example.example_carma("test")
   carma.run()

if everything is running correctly, a directory named "test" should have been created in the directory where you ran this code.  Additionally you should see that the model began running as it prints its current timestep—you may stop the run once you see it is running but if you let it run uninterrupted, the code should finish running in a few minutes.  Feel free to delete the "test" directory once you are sure everything is running correctly
