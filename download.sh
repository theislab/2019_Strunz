#!/bin/bash

# downloads 9Gb dataset files
wget https://hmgubox.helmholtz-muenchen.de/f/492f4319237a464f9a28/?dl=1 -O 2019_Strunz.tar
mkdir data; cd data
tar -xvf ../2019_Strunz.tar
