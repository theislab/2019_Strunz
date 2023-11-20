#!/bin/bash

# downloads 9Gb dataset files
wget https://1drv.ms/f/s!AiwuS6Fz15R0gsgiXtjxJZqr-BDt0A?e=nCnnJ5 -O 2019_Strunz.tar
mkdir data; cd data
tar -xvf ../2019_Strunz.tar
