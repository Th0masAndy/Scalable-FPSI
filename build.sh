#!/usr/bin/env bash

apt update
apt install -y build-essential cmake git libtool iproute2 python3 sudo nasm libssl-dev libgmp-dev wget libfmt-dev

cd thirdparty
mkdir -p out/src
cd out/src

## Build volePSI (libOTe and other dependencies are built automatically.)
git clone https://github.com/ladnir/volepsi.git
cd volepsi
git checkout 59e06bca81a3287257522cd261bad71e37780642
python3 build.py --install=../../install --system -DVOLE_PSI_ENABLE_BOOST=true -DVOLE_PSI_ENABLE_OPPRF=true -DVOLE_PSI_ENABLE_BITPOLYMUL=false -DVOLE_PSI_SODIUM_MONTGOMERY=false -DFETCH_SPARSEHASH=true -DCMAKE_PREFIX_PATH=/usr/local/
cp ./out/build/linux/volePSI/config.h ../../install/include/volePSI
cd ..

## Build BLAKE3
git clone https://github.com/BLAKE3-team/BLAKE3.git
cd BLAKE3
git checkout c7f0d216e6fc834b742456b39546c9835baa1277
cmake -S c -B c/build -DCMAKE_INSTALL_PREFIX=../../install
cmake --build c/build --target install -j
cd ..
rm -rf BLAKE3

cd ../../..
