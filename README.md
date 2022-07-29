```
 __________.__             ___.     |__|  __
 \______   \__|_____   ____\_ |__   _||__/  |_  ______ (C) George Jackson-Mills 2020
  |     ___/  \____ \_/ __ \| __ \ /  _ \   __\/  ___/
  |    |   |  |  |_> >  ___/| \_\ (  O_O )  |  \___ \
  |____|   |__|   __/ \___  >___  /\____/|__| /____  >
              |__|        \/    \/                 \/
```

# itur - Implementation of several ITU-R Recommendations

## Overview

This module contains a Python implementation of engineering formulas in several ITU-R Recommendations that are relevant and applicable to Pipebots. In particular, we have got:

- P.525 - Free-space path loss, also known as Friis equation. Used to calculate baseline radio wave propagation losses at short distances in pipes and through soils.
- P.527 - Electromagnetic properties of arbitrary soils, as well as those of fresh and salt water. Necessary for estimating attenuation in pipes when considered as lossy waveguides.
- P.2040 - Electromagnetic properties of building materials and structures. A succint summary of the many papers on rectangular lossy waveguides, as well as propagation losses within dielectric materials such as concrete, wood, brick, and so on.

## Requirements

Requires some standard Python packages for scientific computing. As my current development environment is a bit polluted, I'll list the required packages here. Apologies.

- `Python>=3.6`
- `numpy`
- `scipy`

## Tests

To be added. Bad practice, I know.

## Installation

Use `pip install -e .` in the folder to which you clone or download this. This will install `rflib` as an "editable" package in your current environment, meaning you could just do a `git pull` in the future to get any updates.

## Contributing

Contributions are more than welcome and are in fact actively sought! Please contact Viktor at [eenvdo@leeds.ac.uk](mailto:eenvdo@leeds.ac.uk).

## Acknowledgements

This work is supported by the UK's Engineering and Physical Sciences Research Council (EPSRC) Programme Grant EP/S016813/1.
