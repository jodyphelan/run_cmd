import setuptools
import glob

version = [l.strip() for l in open("run_cmd/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(

	name="run_cmd",
	version=version,
	packages=["run_cmd"],
	license="GPL3",
	long_description="Wrapper functions for running commands",
)
