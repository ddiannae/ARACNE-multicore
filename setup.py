from setuptools import setup, find_packages

setup(
        name="MultiAracne", 
        version='0.0.2',
        author="Diana García-Cortés",
        author_email="diana.gco@gmail.com",
        description="Python package to run a parallel version of Aracne",
        packages=["MultiAracne"],
        install_requires=[
		'pandas', 
		'numpy'],
        keywords=[],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Operating System :: Unix",
			"Topic :: Utilities"
        ]
)
