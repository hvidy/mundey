# mundey
[![Licence](http://img.shields.io/badge/license-GPLv3-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
![](https://github.com/hvidy/mundey/workflows/integration/badge.svg)
[![PyPI version](https://badge.fury.io/py/mundey.svg)](https://badge.fury.io/py/mundey)

Saving Sirius - Recalibration of TESS target pixel files for the brightest stars

## Contributors

Tim White, Benjamin Pope

## Installation

The easiest option is just from PyPI:

`pip install mundey`

To install from source: clone this git repo, enter the directory, and run

`python setup.py install`

## Smear calibration problems for bright stars observed by TESS

Bright stars saturate the TESS cameras, resulting in excess charge bleeding along CCD columns. For particularly bright stars, and those high in the science frame and/or the upper buffer rows, this excess charge can bleed into the upper serial register and result in an overestimated smear correction. This not only renders data for these stars unusable, but also for any star that appears in the same CCD column. The image below shows one frame from the Sector 6 target pixel file for Sirius. Note that the aspect ratio has been adjusted because the image is extremely long (1287 pixels x 11 pixels). The colour scale is logarithmic; pixles that appear white are where the flux is negative as a result of the overestimated smear correction.

![Before image of Sirius](/images/Sirius-before.png)

## Using mundey

This code, `mundey`, recalibrates TESS target pixel files, including a new smear correction, calculated from the target pixel file itself, that circumvents this problem. It is simple to run mundey:

```python
from mundey.mundey import mundey_tpf

tpf_filename = 'tess2018349182459-s0006-0000000322899250-0126-s_tp.fits'

tpf = mundey_tpf(tpf_filename)

tpf.calibrate()
```

The `mundey_tpf` class inherits the methods and attributes of the `lightkurve.TessTargetPixelFile` class.

Besides the target pixel file, `mundey` also requires several files for calibration, most of which are available from MAST (https://archive.stsci.edu/tess/all_products.html). In particular, `mundey` will require the Trailing Virtual Column and Smear Row collateral data, the Two-Dimensional Black Model, the Linearity Model, and the Flat Field Model. These are available at MAST. Additionally, a table of undershoot correction factors, `tess-undershoot.xml` is also required, and is included in this repo.

After running `mundey`, the target pixel file will have flux values updated with the new calibration. To save this updated tpf, use the `tpf.tofits([output_fn, overwrite])` method of the `TessTargetPixelFile` class.

```python
tpf.to_fits(filename)
```

The image below shows the Sirius target pixel file after correction with `mundey`.

![After image of Sirius](/images/Sirius-after.png)

## License

We invite anyone interested to use and modify this code under a GPL v3 license. 

## Jack Mundey

Jack Mundey (1929-2020) was an Australian trade unionist and environmental activist. Rising to prominance as a leader of the New South Wales Builders' Labourers Federation (BLF) during the late 1960s and early 1970s, Mundey saw his role as not only advocating for improvement to the pay and working conditions of his comrades, but to ensure that their labour was used for the good of their society. To this end, and with the support of community groups, the BLF instituted a series of **green bans** on developments that would destroy the natural and built environment of Sydney. This is the original use of the word "green" in reference to environmental politics, a label that was taken up by Petra Kelly, co-founder of the German Green Party, following her visit to Sydney in 1977. Among the areas saved by green bans is The Rocks district of Sydney, which contains numerous historic buildings and was home to a vibrant working-class community.

The Sirius Building is an apartment complex in The Rocks that was built as public housing for local community members as a result of the green bans, and was completed in 1980. In 2015, the NSW state government made a decision to sell the building, driving out the working class and elderly residents, and placing the building at risk of demolition by the new owners. This decision prompted community protest and a new green ban was placed on the site by the CFMEU. While this action has saved the Sirius building from destruction, the sale eventually proceeded and the Sirius will no longer be public housing. It is apparently intolerable to the government that poor people may live in areas of high real estate value.

![Save Our Sirius protest, 17 September 2016. Photo: Tim White](/images/SaveOurSirius-2016-09-17.jpg)
Save Our Sirius protest march outside the Sirius Building, 17 September 2016. Photo: Tim White

While trying to save observations of the star Sirius from an overestimated smear correction obviously pales in importance to trying to save the community of residents in the Sirius Building, we have chosen to name this code `mundey` in memory of Jack Mundey and the NSW BLF.