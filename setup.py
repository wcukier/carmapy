from setuptools import setup
from setuptools.dist import Distribution
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py

# bdist_wheel was vendored into setuptools (>=70.1); fall back to the wheel
# package for older toolchains.
try:
    from setuptools.command.bdist_wheel import bdist_wheel
except ImportError:  # pragma: no cover
    from wheel.bdist_wheel import bdist_wheel

import subprocess
import os
import shutil
import sys

SRC = os.path.join(os.path.dirname(__file__), "src")


class BinaryDistribution(Distribution):
    """Force setuptools to treat carmapy as an impure (platform) distribution.

    carmapy ships a compiled ``carmapy.exe`` but declares no ``ext_modules``,
    so setuptools would otherwise classify the package as pure and install it
    into *purelib*. With the wheel root declared platlib (see ``BinaryWheel``),
    that mismatch routes every file into a ``*.data/purelib/`` scheme dir, and
    auditwheel/delocate reject the binary ("shared library in purelib folder").
    Reporting ext modules here installs the package to *platlib*, so files land
    at the wheel root where the repair tools expect them.
    """

    def has_ext_modules(self):
        return True


class BuildFortranBinary(build_ext):
    def run(self):
        print("Building Fortran binary...")

        # Try Intel compiler first
        env = os.environ.copy()
        # Pass NATIVE through to the Makefile. Default 1 (source installs build for
        # the local CPU). Wheel builds set CARMAPY_NATIVE=0 for portable binaries.
        native = os.environ.get("CARMAPY_NATIVE", "1")
        openmp = os.environ.get("CARMAPY_OPENMP", "0")
        compilers = [
            ("ifort", ["make", "all", "FORTRAN=ifort", f"NATIVE={native}", f"OPENMP={openmp}"]),
            ("gfortran", ["make", "all", "FORTRAN=gfortran", f"NATIVE={native}", f"OPENMP={openmp}"]),
        ]
        built = False

        os.makedirs(os.path.join(SRC, "CARMA", "build", "carma"), exist_ok=True)
        shutil.copyfile(os.path.join(SRC, "CARMA", "Makefile"), os.path.join(SRC, "CARMA", "build", "carma", "Makefile"))

        try:
            subprocess.call(["make", "clean"], cwd=os.path.join(SRC, "CARMA", "build", "carma"))
        except Exception:
            pass

        for name, cmd in compilers:
            try:
                print(f"Trying to build with {name}")


                subprocess.check_call(cmd, cwd=os.path.join(SRC, "CARMA", "build", "carma"), env=env)
                built = True
                break
            except subprocess.CalledProcessError as e:
                print(e)
                print(f"{name} failed")

        if not built:
            raise RuntimeError("Fortran build failed with both Intel and gfortran compilers.")
        # Move binary to package dir
        binary_path = os.path.join(SRC, "CARMA", "build", "carma", "carmapy.exe")  # adjust to your binary name
        print(self.build_lib)
        target_path = os.path.join(self.build_lib, "carmapy", "carmapy.exe")
        os.makedirs(os.path.dirname(target_path), exist_ok=True)
        shutil.copy2(binary_path, target_path)
        print(f"Binary copied to {target_path}")
        # raise


class CustomBuildPy(build_py):
    def run(self):
        self.run_command("build_ext")
        super().run()


class BinaryWheel(bdist_wheel):
    """Tag the wheel ``py3-none-<platform>``.

    carmapy bundles a standalone Fortran executable (carmapy.exe) that is
    invoked via subprocess and does NOT link libpython, so the wheel is
    platform-specific (right OS/arch) but Python-version agnostic. Marking it
    impure but py3/none yields one wheel per platform that works on any
    Python 3, instead of a separate cp310/cp311/cp312/cp313 wheel each.
    """

    def finalize_options(self):
        super().finalize_options()
        # Wheel carries a compiled binary: not a pure-Python ("any") wheel.
        self.root_is_pure = False

    def get_tag(self):
        # Keep the platform tag computed by setuptools; relabel interpreter/ABI
        # as py3/none since the binary is independent of the Python version.
        _python, _abi, plat = super().get_tag()
        return "py3", "none", plat


setup(
    distclass=BinaryDistribution,
    cmdclass={
        "build_ext": BuildFortranBinary,
        "build_py": CustomBuildPy,
        "bdist_wheel": BinaryWheel,
    },
)




