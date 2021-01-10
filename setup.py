import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="epos", # Replace with your own username
    version="0.0.1",
    author="David Fuhrlaender",
    author_email="d.fuhrlaender@uni-bremen.de",
    description="EpoS - Electrolysis plant operation Simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        #"License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
