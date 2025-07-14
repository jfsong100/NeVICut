from setuptools import setup, find_packages

setup(
    name='cutbayesflow',
    version='0.1.0',
    description='Cut Bayes Variational Inference with Normalizing Flows',
    author='Jiafang Song',
    packages=find_packages(),
    install_requires=[
        'torch>=2.0',
        'nflows>=0.14'
    ],
)
