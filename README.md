# pyDrizzlePac
(updated on 2020. 07. 02.)


## Description
A Python package for automatically drizzling the raw data of the Hubble Space Telescope


## Prerequisites
* You have to retrieve the HST raw images with ``*_flc.fits`` or ``*_flt.fits`` from [the Hubble MAST archive](http://archive.stsci.edu/hst/search.php).
* The following Python modules should be installed.
  * ``numpy <= 1.17.0``
  * ``matplotlib >= 3.2.2``
  * ``pandas >= 1.0.5``
  * ``astropy >= 4.0.0``
  * ``astroscrappy >= 1.0.5`` ([Reference link](https://astroscrappy.readthedocs.io/en/latest/))
  * ``sep >= 1.0.3`` ([Reference link](https://sep.readthedocs.io/en/v1.0.x/))
  * ``stwcs >= 1.5.3`` ([Reference link](https://stwcs.readthedocs.io/en/latest/hstwcs.html))
* The ``astroconda`` environement which contains ``DrizzlePac`` module should be installed. ([Reference link](https://astroconda.readthedocs.io/en/latest/getting_started.html#))
  * ``drizzlepac >= 3.1.6`` ([Reference link](https://drizzlepac.readthedocs.io/en/latest/))
* [init_param.py](https://github.com/joungh93/pyDrizzlePac/blob/master/init_param.py) is the initial configurations to run the drizzle tasks. (You can revise it!)


## Workflows
```
cd /your_working_directory/
git clone https://github.com/joungh93/pyDrizzlePac.git
```
You should activate the ``astroconda`` environment which contains ``DrizzlePac`` module.

```
conda activate astroconda
```

After revising ``init_param.py``, the following simple commands will work well.

(Still you have to check whether the resulting drizzled images are good or not.)

```
ipython
run start_drz.py
run opt_1.py
run mat_1.py
run twk_1.py
run drz_1.py
```

The above python codes can be also executed by a simple shell script command.

```
vi run.sh
  python *.py
sh run.sh
```


## Future works
* Current version sets the reference image to the shot with _**the first observation date**_.
This setting is generally good enough to get well-drizzled image, but it could get some troubles if the image with the first observation date has bad quality.
It can be possible in the future version to choose the reference image as one wants to use.
* Iterations of drizzling might be added.

To update this code is an on-going task... :crying_cat_face: :sweat_drops:
