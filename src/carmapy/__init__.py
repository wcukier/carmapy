from .carmapy import Carma, Gas, Group, Element, Coag, Growth, Nuc



import carmapy.example

from .results import *
from .base import *
from .utils import *
from .constants import *

from importlib.metadata import version
__version__ = version('carmapy')

