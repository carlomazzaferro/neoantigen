from setuptools import setup

setup(name='variantannotation',
      version='0.17',
      description='This package is aimed at providing a way of retrieving variant information using ANNOVAR '
                  'http://annovar.openbioinformatics.org/en/latest/ and myvariant.info. In particular, it is'
                  'suited for bioinformaticians interested in aggregating variant information into a single database.'
                  'The data structure used is a python dictionaty, so the data can be easily parsed to a MongoDB'
                  'instance.',
      url='https://github.com/Mazzafish/variantannotation',
      author='Carlo Mazzaferro',
      author_email='cmazzafe@ucsd.edu',
      install_requires = ['pymongo', 'pysam', 'myvariant', 'pyvcf'],
      license='MIT',
      packages=['variantannotation'],
      zip_safe=False)