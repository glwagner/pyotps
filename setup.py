from setuptools import setup

setup(
    name = 'pyotps',
    version = '0.1',
    description = "A Python wrapper around the Oregon Tidal Prediction System..", 
    url = 'http://github.com/glwagner/pyotps',
    author = 'Gregory L. Wagner',
    author_email = 'wagner.greg@gmail.com',
    license = 'MIT',
    packages = ['pyotps'],
    install_requires = [
        'numpy', 
        'matplotlib', 
    ],
    zip_safe = False,
)
