mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make network
mv network $PREFIX/bin/afines
