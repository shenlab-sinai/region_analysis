from setuptools import setup

# Hack to prevent stupid TypeError: 'NoneType' object is not callable error on
# exit of python setup.py test # in multiprocessing/util.py _exit_function when
# running python setup.py test (see
# http://www.eby-sarna.com/pipermail/peak/2010-May/003357.html)
# copy from https://github.com/pypa/virtualenv/blob/develop/setup.py
try:
    import multiprocessing
except ImportError:
    pass

setup(name='regionanalysis',
      version='0.1.1',
      description='A utility to annotate genomic intervals.',
      url='http://github.com/shenlab-sinai/region_analysis',
      author='Ning-Yi SHAO',
      author_email='shaoningyi@gmail.com',
      license='GPLv3',
      packages=['regionanalysis'],
      scripts=['bin/region_analysis.py', 'bin/region_analysis_db.py'],
      include_package_data=True,
      long_description=open("README.rst").read(),
      test_suite='nose.collector',
      tests_require=['nose'],
      install_requires=[
          'pybedtools',
      ],
      zip_safe=False)