# Validation of Variance for Sobol' and CvM rank-based estimators

## Content

The repository is structured as follows. We only describe the most important files for a new user.
```bash
./
|-- ipynb: Contains Python notebooks which demonstrate how the code works
|-- README.md: This file
```

## Default installation and env setup

1. Create new virtual environment

```bash
$ python3 -m venv .venv_active_sampling_flow
```

(Do
sudo apt install python3-venv
if needed)

3. Activate virtual environment

```bash
$ source .venv_active_sampling_flow/bin/activate
```

4. Upgrade pip, wheel and setuptools 

```bash
$ pip install --upgrade pip
$ pip install --upgrade setuptools
$ pip install wheel
```

5. Install the `<name>` package.

```bash
python setup.py develop
```

6. (Optional) In order to use Jupyter with this virtual environment .venv
```bash
pip install ipykernel
python -m ipykernel install --user --name=.venv_active_sampling_flow
```
(see https://janakiev.com/blog/jupyter-virtual-envs/ for details)

## Specific to this package
Note that the dependencies have been left to a bare minimum in order to run the package. One can do

```bash
$ pip install tqdm numpy matplotlib sympy
```


## Configuration
Nothing to do

## Credits
The authors would like to thanks [...]
