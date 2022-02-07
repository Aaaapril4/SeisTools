import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Seistools",
    version="0.0.1",
    author="YJie",
    author_email="jieyaqi@msu.edu",
    description="Modules for daily use in seismology",
    url="https://github.com/Aaaapril4/SeisTools",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
)