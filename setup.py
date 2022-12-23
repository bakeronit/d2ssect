from setuptools import Extension,setup

setup(
	ext_modules=[
        Extension(
            name="d2ssect.jellyfish",
            sources=["d2ssect/jellyfish.cpp"], # all sources are compiled into a single binary file
            include_dirs  = ['/Users/iracooke/.local/include/jellyfish-2.3.0'],
            library_dirs = ['/Users/iracooke/.local/lib'],
            libraries = ['jellyfish-2.0.2'],
            language='c++',
        ),
    ]
)