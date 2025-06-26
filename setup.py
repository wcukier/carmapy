from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
import subprocess
import os
import shutil
import sys

class BuildFortranBinary(build_ext):
    def run(self):
        print("Building Fortran binary...")

        # Try Intel compiler first
        env = os.environ.copy()
        compilers = [
            ("ifort", ["./make-carma.csh", "all", "ifort"]),
            ("gfortran", ["./make-carma.csh", "all", "gfortran"]),
        ]
        built = False
        for name, cmd in compilers:
            try:
                print(f"Trying to build with {name}")
                print(os.getcwd())
                subprocess.check_call(cmd, cwd=os.path.join("src", "CARMA"), env=env)
                built = True
                break
            except subprocess.CalledProcessError as e:
                print(e)
                print(f"{name} failed")

        if not built:
            raise RuntimeError("Fortran build failed with both Intel and gfortran compilers.")
        # Move binary to package dir
        # binary_path = os.path.join("CARMA", "carmapy.exe")  # adjust to your binary name
        # target_path = os.path.join(self.build_lib, "mypackage", "mybinary")
        # os.makedirs(os.path.dirname(target_path), exist_ok=True)
        # shutil.copy2(binary_path, target_path)
        # print(f"Binary copied to {target_path}")


class CustomBuildPy(build_py):
    def run(self):
        self.run_command("build_ext")
        super().run()

setup(    cmdclass={
        "build_ext": BuildFortranBinary,
        "build_py": CustomBuildPy,
    })