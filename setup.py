import setuptools

setuptools.setup(
    name="CineMol",
    version="0.0.1",
    author="David Meijer",
    author_email="david.meijer@wur.nl",
    install_requires=[],
    package_dir={"": "src"},
    packages=["cinemol"],
    python_requires=">=3.10",
    entry_points={"console_scripts": ["cinemol = main:main"]}
)