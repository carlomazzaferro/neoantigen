from setuptools import setup

setup(name='VarP',
      version='0.3',
      description='This package is aimed at providing a way of retrieving epitope candidates from a vcf file'
                  'instance.',
      url='https://github.com/Mazzafish/neoantigen',
      author='Carlo Mazzaferro',
      author_email='cmazzafe@ucsd.edu',
      install_requires='varcode',
      license='MIT',
      packages=['VarP'],
      zip_safe=False)