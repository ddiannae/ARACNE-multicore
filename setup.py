from setuptools import setup, find_packages

setup(
        name="MultiAracne", 
        version='0.0.1',
        author="Diana García-Cortés",
        author_email="diana.gco@gmail.com",
        description="Python package to run a parallel version of Aracne",
        packages=find_packages(),
        install_requires=[
		'pandas', 
		'numpy'], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        keywords=[],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Operating System :: Unix",
			"Topic :: Utilities"
        ]
)
