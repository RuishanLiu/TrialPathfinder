try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name = 'TrialPathfinder',
    version = '1.1.0',
    packages=['TrialPathfinder'],
    license='MIT',
    description = 'Python library for systematic evaluation of clinical trial eligibility criteria.',
    url = 'https://github.com/RuishanLiu/TrialPathfinder',
    author = 'Ruishan Liu',
    author_email = 'ruishan@stanford.edu',
    install_requires=[
        'numpy',
        'pandas',
        'lifelines',
        'sklearn'
    ],
)
