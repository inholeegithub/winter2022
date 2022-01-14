# XRD
X-ray diffraction calculations.

This branch is initally planned for my [Computational Physics class in 2017 fall](https://github.com/qzhu2017/2017-cmp)

Currently, there are four classes,
- Element
- crystal
- cif
- XRD

One could load the crystal from 
- dictionary
- POSCAR
- CIF file 

To perform XRD calculation, one needs to provide the following info
- crystal structure
- wavelength (default is Cu-Ka: 1.54184 \AA)
- maximum 2\theta value (defult: 180 degree)

The atomic scattering factor is calculated from 9-parameter equation by Don Cromer and J. Mann.

More detailed usage could be found in the [jupyter notebook](https://github.com/qzhu2017/XRD/blob/master/Demo.ipynb).

## usage
```
$ python XRD.py -h
Usage: XRD.py [options]

Options:
  -h, --help            show this help message and exit
  -m hkl index, --hkl=hkl index
                        show hkl_index info, e.g., [1,0,0]
  -a angle, --angle=angle
                        2theta angle range, default=180
  -t files, --transform=files
                        export file in different format
  -p plot, --plot=plot  plot pxrd, default: yes
  -w wave, --wavelength=wave
                        wavelength: 1.54184
  -c crystal, --crystal=crystal
                        crystal from file, cif or poscar, REQUIRED
  -f full, --full=full  show full hkl reflections
  -i intensity, --intensity=intensity
                        the minimum intensity to show, default 0.01
 ```
 ## execute 
 one just needs to run the followings,
```
$ python XRD.py -c NaCl.cif
      2theta    Intensity     d_hkl    h    k    l
--  --------  -----------  --------  ---  ---  ---
 0   27.489     0.0766015  3.24471     1    1    1
 1   31.8464    1          2.81        0    0   -2
 2   45.6587    0.698779   1.98697     0    2    2
 3   54.1242    0.0692494  1.69449    -1    3    1
 4   56.7429    0.241953   1.62235    -2   -2   -2
 5   66.5554    0.112971   1.405       0    0    4
 6   73.4435    0.0377815  1.28932     3   -1   -3
 7   75.6806    0.317616   1.25667     2   -4    0
 8   84.4455    0.246837   1.14718    -2    2    4
 9   90.9229    0.0368239  1.08157     1    5    1
10  101.787     0.0975993  0.993485    0    4   -4
11  108.492     0.0539193  0.949953   -5   -3   -1
12  110.782     0.245344   0.936667    2    4    4
13  120.354     0.213082   0.8886     -6    0   -2
14  128.188     0.0349423  0.857042    5    3    3
15  130.986     0.25169    0.847247    2    2   -6
16  143.747     0.112952   0.811177   -4   -4   -4
17  156.826     0.121099   0.786957   -5   -5    1
18  163.127     0.737297   0.779354   -6    4    0
```
It will also generate a png file with PXRD plot as follows
![NaCl](https://github.com/qzhu2017/XRD/blob/master/images/NaCl.cif.png)

