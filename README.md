# pyDrizzlePac
(updated on 2020. 06. 30.)

## Prerequisites
* You have to retrieve the HST raw images with ``*_flc.fits`` or ``*_flt.fits`` from [the Hubble MAST archive](http://archive.stsci.edu/hst/search.php).
* [init_param.py](https://github.com/joungh93/pyDrizzlePac/blob/master/init_param.py) is the initial configurations to run the drizzle tasks. (You can revise it!)

## Workflow
```
cd /your_working_directory/
git clone https://github.com/joungh93/pyDrizzlePac.git
```
You should activate the ``astroconda`` environment which contains ``DrizzlePac`` module.
* Reference links
  * [``astroconda``](https://astroconda.readthedocs.io/en/latest/getting_started.html#)
  * [``DrizzlePac``](https://drizzlepac.readthedocs.io/en/latest/)

```
conda activate astroconda
python start_drz.py
python opt_1.py
python mat_1.py
python twk_1.py
python drz_1.py
```
