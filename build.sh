mkdir -p thirdparty
cd thirdparty

## Build volePSI (We have included the modified version of volePSI in the thirdparty folder, so you can skip this step.)

# git clone https://github.com/ladnir/volepsi.git
# cd volepsi
# git checkout 59e06bca81a3287257522cd261bad71e37780642
# sed -i 's|36cd7242e085eddba34feaa63733ec4c6ded66c7|d21bc4d7aae941e276b92615252fd1760c902890|g' thirdparty/getLibOTe.cmake
cd volepsi
python3 build.py --install=../out/install --system -DVOLE_PSI_ENABLE_BOOST=true -DVOLE_PSI_ENABLE_OPPRF=true -DVOLE_PSI_ENABLE_BITPOLYMUL=false -DVOLE_PSI_SODIUM_MONTGOMERY=false -DFETCH_SPARSEHASH=true -DCMAKE_PREFIX_PATH=/usr/local/
cp ./out/build/linux/volePSI/config.h ../out/install/include/volePSI
cd ..


## Build BLAKE3
git clone https://github.com/BLAKE3-team/BLAKE3.git
cd BLAKE3
git checkout c7f0d216e6fc834b742456b39546c9835baa1277
cmake -S c -B c/build -DCMAKE_INSTALL_PREFIX=../out/install
cmake --build c/build --target install -j
cd ..
rm -rf BLAKE3