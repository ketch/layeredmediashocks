#!/bin/csh
cd clawpack
make
cd ../sharpclaw
make
cd ..
cp clawpack/xclaw sharpclaw/xsclaw .
