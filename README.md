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

This module contains the Python implementation of several ITU-R Recommendations that are relevant and applicable to Pipebots. In particular, we've got:

- P525 - Free-space attenuation. Baseline losses at short distances in pipes and through soils
- P527 - Electromagnetic properties of soil. Necessary for estimating attenuation in pipes when viewed as lossy waveguides.
- P2040 - Building materials and structures. A succint summary of the many papers on rectangular lossy waveguides, as well as propagation losses within dielectric materials such as concrete, wood, brick, and so on.

## Requirements

Requires some standard Python packages for scientific computing. As my current development environment is a bit polluted, I'll list the required packages here. Apologies.

- `Python>=3.6`
- `numpy`
- `scipy`

## Tests

To be added. Bad practice, I know.

## Installation

Use `pip install -e .` in the folder to which you clone or download this. This will install `rflib` as an "editable" package in your current environment, meaning you should just do a `git pull` in the future to get any updates.

## Contributing

Contributions are more than welcome and are in fact actively sought! Please contact Viktor either at [eenvdo@leeds.ac.uk](mailto:eenvdo@leeds.ac.uk) or on the `Pipebots` Slack.