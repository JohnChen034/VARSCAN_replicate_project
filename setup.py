from setuptools import setup
from setuptools import find_packages

setup(
    name='Vcall',
    version='0.1',
    description='Tool developed for snp variant calling',
    author='Kaifu Yang, Jiayu Chen',
    author_email='kay002@ucsd.edu, jic034@ucsd.edu',
    packages=find_packages(),
    install_requires=['Pandas', 'scipy'],
    url = "https://github.com/JohnChen034/CSE185_project",
    entry_points = {
        "console_scripts" : [
                "Vcallsnp = Vcall.main:main",
        ],
    },
)