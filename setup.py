from setuptools import setup, Extension, find_packages
import os
from subprocess import check_output, CalledProcessError


VERSION = "0.2.1"

try:
    RATROOT = os.environ["RATROOT"]
except OSError:
    print(
        "could not detect RATROOT. Make sure you have source 'ratpac.sh' before installing"
    )

# Get RATPAC info
rat_incdir = os.path.join(RATROOT, "include")
rat_lib = "RATEvent"
rat_libdir = os.path.join(RATROOT, "lib")

# Get root info
try:
    root_args = check_output(["root-config", "--cflags", "--libs"]).decode().split()
except CalledProcessError:
    print(
        "Error calling 'root-config'. Make sure you have source 'thisroot.sh' before installing"
    )

extensions = [
    Extension(
        "fastrat",
        ["sibyl_cpp/fastrat.cpp"],
        include_dirs=[rat_incdir],
        libraries=[rat_lib],
        library_dirs=[rat_libdir],
        extra_compile_args=root_args,
    )
]


setup(
    name="sibyl",
    version=VERSION,
    author="Morgan Askins",
    author_email="maskins@berkeley.edu",
    scripts=["scripts/sibyl"],
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    install_requires=[
        "matplotlib",
        "numpy",
        "pyqt5",
        "pyqtgraph",
        "pyopengl",
        "markdown",
    ],
    install_package_data=True,
    package_data={"sibyl": ["assets/*"]},
    zip_safe=False,
    ext_package="sibyl_cpp",
    ext_modules=extensions,
)
