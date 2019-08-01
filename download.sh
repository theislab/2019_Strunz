#!/bin/bash

# downloads XXXGb dataset files
wget https://hmgubox.helmholtz-muenchen.de/XXXXXXX -O 2019_Strunz.tar
mkdir data; cd data
tar -xvf ../2019_Strunz.tar
