import setuptools
setuptools.setup(
    name='PyGSI',
    version='1.0.0',
    description='Tools to read and filter through GSI diagnostic files.',
    author='NOAA-EMC',
    author_email='Kevin.Dougherty@noaa.gov',
    url='https://github.com/noaa-emc/PyGSI',
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where='src'),
    classifiers=[
        'Development Status :: 1 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Operating System :: OS Independent',
        'Typing :: Typed'],
    python_requires='>=3.6',
    install_requires=[
        'pyyaml>=6.0',
        'pycodestyle>=2.8.0',
        'netCDF4>=1.5.3',
        'matplotlib>=3.5.2',
        'cartopy>=0.20.2',
        'scikit-learn>=1.0.2',
        'xarray>=0.11.3',
        'emcpy @ git+https://github.com/NOAA-EMC/' +
        'emcpy@4f36baf1a2302fb0daa49bd8415bb7d2a65347bb#egg=emcpy'
    ]
)
