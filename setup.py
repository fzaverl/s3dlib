import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="s3dlib",
    version="1.0.0",
    author="Frank Zaverl, Jr.",
    author_email="fzaverl@s3dlib.org",
    keywords='python matplotlib 3D surface visualization',
    description="Python classes to create 3D surface objects rendered in Matplotlib",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        'Documentation': 'https://s3dlib.org/',
        'Source': 'https://github.com/fzaverl/s3dlib/',
    },
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Visualization",
        "Intended Audience :: Science/Research",
    ],
    python_requires='>=3.6',
)