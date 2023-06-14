# VenomPred 2.0 API
**VenomPred 2.0** platform allows you, through an innovative machine learning consensus strategy, to evaluate  the **toxicological profile** of one or multiple small molecules 

## Requirements

You must have the conda package manager installed
Install the virtual environment 
```
$ conda env create -f enviroment.ym
```

## Usage

Clone the repository and run the script.
```
$ python VenomPred2_predictions.py -in test.csv -o res.csv 
```
You need to provide a csv file with a column named SMILES containing the SMILES structure of the compounds you want to predict

