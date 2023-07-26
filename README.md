# VenomPred 2.0 API
 ![](http://www.mmvsl.it/wp/wp-content/uploads/2022/08/mmvsl_2022_lr.png) ![](http://www.mmvsl.it/wp/wp-content/uploads/2022/02/logo_web.png) 
 
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
![Static Badge](https://img.shields.io/badge/MMVSL-lab-blue?style=for-the-badge&link=http%3A%2F%2Fwww.mmvsl.it%2Fwp%2F)


**VenomPred 2.0** platform allows you, through an innovative machine learning consensus strategy, to evaluate  the **toxicological profile** of one or multiple small molecules 

## Requirements

You must have the conda package manager installed.
Install the virtual environment.
```
$ conda env create -f enviroment.yml
```

## Usage

Clone the repository and run the script.
```
$ python VenomPred2_predictions.py -in test.csv -o res.csv 
```
You need to provide a csv file with a column named SMILES containing the SMILES structure of the compounds you want to predict

## License

For this project the default copyright laws apply. All rights of the source code cannot be reproduced, distributed, or used to create derivative works.
The users can use the code to obtain predictions for including them in their work. 
