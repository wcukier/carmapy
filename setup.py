from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
import subprocess
import os
import shutil
import sys

SRC = os.path.join(os.path.dirname(__file__), "src")


class BuildFortranBinary(build_ext):
    def run(self):
        print("Building Fortran binary...")

        # Try Intel compiler first
        env = os.environ.copy()
        compilers = [
            ("ifort", ["make", "all", "FORTRAN=ifort"]),
            ("gfortran", ["make", "all", "FORTRAN=gfortran"]),
        ]
        built = False

        os.makedirs(os.path.join(SRC, "CARMA", "build", "carma"), exist_ok=True)
        shutil.copyfile(os.path.join(SRC, "CARMA", "Makefile"), os.path.join(SRC, "CARMA", "build", "carma", "Makefile"))

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

setup(    cmdclass={
        "build_ext": BuildFortranBinary,
        "build_py": CustomBuildPy,
    })