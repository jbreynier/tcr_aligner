from setuptools import setup

setup(
    name='tcraligner',
    version='0.1.0',    
    description='Align TCR data ',
    url='https://github.com/jbreynier/tcr_aligner',
    author='Johnny B',
    author_email='jb@jb.com',
    license='BSD 2-clause',
    packages=['tcraligner'],
    install_requires=['biopython',
                      'numpy',  
                      'pandas'                   
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux'
    ],
)