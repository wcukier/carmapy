Installation
===============

Installation with git
---------------------
Currently this is the only way to install carmapy::

   git clone https://github.com/wcukier/carmapy
   cd carmapy
   pip install -e .


Testing Installation
--------------------
To test that carmapy is installed correctly you can run the following python code

.. code-block:: python

   import carmapy.example

   carma = carmapy.example.example_carma("test")
   carma.run()

if everything is running correctly, a directory named "test" should have been created in the directory where you ran this code.  Additionally you should see that the model began running as it prints its current timestep—you may stop the run once you see it is running but if you let it run uninterrupted, the code should finish running in a few minutes.  Feel free to delete the "test" directory once you are sure everything is running correctly
