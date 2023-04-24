# PPII Predictor
@Shashank Pritam - (shashankpritam[at]gmail[dot]com)
This repository contains the code for predicting the binding sites of Poly-L-Proline Helix (PPII) peptides in protein structures. PPII helices are a class of peptides that play a crucial role in cellular functions such as protein-nucleic acid interaction and protein assembly. They are often present on various protein binding surfaces and act as common binding motifs responsible for various biological activities. This study aims to predict PPII binding sites and analyze their structural conservation across different protein families.
This program takes a pdb structure as input and return possible binding site of the PPII Helix.
Input - 4 letter PDB ID
Output - In case the given structure does not contain a PPII binding site - You'll get notified.
In case the given structure does contain a PPII binding site - You'll get the structure with PPII bound to it on the "best possible orientation" at the predicted binding site of the PPII.

For any help please contact through the provided email ID.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Citation](#citation)

## Overview

The PPII Predictor is designed to predict the binding sites of PPII peptides in protein structures. Using a template dataset of 37 protein structures belonging to 8 different structural families, we have devised a superimposition-based prediction method for PPII binding sites. The PPII Predictor can also potentially help create methodologies for predicting other conformation's binding sites.

## Installation

To use the PPII Predictor, follow these steps:

Clone the repository:

```bash

git clone https://github.com/shashankpritam/PPIIPredictor.git

```
## Usage

Run the PPII Predictor:

```rust

cargo build

cargo run "PDB_ID"

```

## License
The PPII Predictor is released under the Lesser GPL, version 2.1. (LGPL-2.1)

## Citation
Please cite the following article when using the PPII Predictor in your research:

Pritam, Shashank. “Prediction of Polyproline Type II Helix Receptors,” 2022. http://dr.iiserpune.ac.in:8080/xmlui/handle/123456789/6766.
