try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name = 'TrialPathfinder',
    version = '1.1.0',
    description = 'Python library for systematic evaluation of clinical trial eligibility criteria.',
    author = 'Ruishan Liu',
    author_email = 'ruishan@stanford.edu',
    url = 'https://github.com/RuishanLiu/TrialPathfinder',
    # download_url = '',
    packages=['TrialPathfinder'],
    keywords = ['Clinical Trial', ''],
    install_requires=[
        'numpy',
        'pandas',
        'lifelines',
        'sklearn'
    ],
)
