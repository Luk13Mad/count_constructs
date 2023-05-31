from setuptools import setup,find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read()
    
setup(
    name = 'count_constructs',
    version = '1.0',
    author = 'Lukas Madenach',
    author_email = 'Lukas.madenach@kitz-heidelberg.de',
    license = 'MIT',
    description = 'Tool producing count tables for 2D CRISPR screen.',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = '<github url where the tool code will remain>',
    py_modules = ['main', 'screen_2D','screen_3D'],
    packages = find_packages(),
    install_requires = [requirements],
    python_requires='>=3.10',
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    entry_points = '''
        [console_scripts]
        count_constructs=src.main:cli
    '''
)
