from setuptools import setup, find_packages

setup(
    name='e3pred',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'e3pred': ['model/*.pkl', 'pfeature_comp.py', 'readme.txt','license.txt','Data/*'],
    },
    install_requires=[
        'scikit-learn',
        'tqdm',
        'pandas', 
        'openpyxl'
    ],
    entry_points={
        'console_scripts': [
            'e3pred=e3pred.module1:main',
        ],
    },
)

