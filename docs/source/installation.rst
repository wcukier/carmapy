Installation
===============

Installation with pip
---------------------
To install with pip simply write::
   
   pip install carmapy

Assuming you are on a common operating system, this should work and all 
dependencies should automatically be installed for you.

Installation with git
---------------------
To install and build loacally::

   git clone --recursive https://github.com/wcukier/carmapy
   cd carmapy
   pip install .


Most dependencies will be automatically installed for you but in order to 
compile carmapy from source you must first either have the intel fortran compiler 
or gfortran installed.  To check that they are installed you can run the 
following commands::

   which ifort
   which gfortran

as long as at least one of those commands returns as file path rather than saying
"ifort not found" or "gfortran not found" then you should be good to proceed with
installation.

OpenMP
------
carmapy allows for running the code on multiple threads with OpenMP.  To install
with OpenMP you will need to compile the code on your computer with the environment
variable `CARMAPY_OPENMP=1`.  To do this with pip and bypassing github run::

   CARMAPY_OPENMP=1 pip install --no-binary=carmapy --force-reinstall carmapy

If installing from github replace the pip install line with::

   CARMAPY_OPENMP=1 pip install --upgrade .

Testing Installation
--------------------
To test that carmapy is installed correctly you can run the following python code

.. code-block:: python

   import carmapy.example

   carma = carmapy.example.example_carma("test")
   carma.run()

if everything is running correctly, a directory named "test" should have been created in the directory where you ran this code.  Additionally you should see that the model began running as it prints its current timestep—you may stop the run once you see it is running but if you let it run uninterrupted, the code should finish running in a few minutes.  Feel free to delete the "test" directory once you are sure everything is running correctly
