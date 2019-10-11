from setuptools import setup, find_packages

setup(
    name = 'demo20191011',
    version = '0.3',
    author = 'abricot',
    author_email = 'zayintopology@gmail.com',
    url = 'https://github.com/kxxdhdn/hello-world',

    install_requires = ['numpy>=1.17.2'],
    
    entry_points={
        # installation test with command line
        'console_scripts': [
            'test_demo = src:test',
        ],
    },
    
    packages = find_packages(),
)
