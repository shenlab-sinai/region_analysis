from setuptools import setup

setup(name='regionanalysis',
      version='0.1.1',
      description='A utility to annotate genomic intervals.',
      url='http://github.com/shenlab-sinai/region_analysis',
      author='Ning-Yi SHAO',
      author_email='shaoningyi@gmail.com',
      license='GPLv3',
      packages=['regionanalysis'],
      scripts=['bin/region_analysis.py'],
      include_package_data=True,
      long_description=open("README.rst").read(),
      test_suite='nose.collector',
      tests_require=['nose'],
      install_requires=[
          'pybedtools',
      ],
      zip_safe=False)