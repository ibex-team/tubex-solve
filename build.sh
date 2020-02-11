# ==================================================================
#  tubex-solve - build script
# ==================================================================

#!/bin/bash

mkdir make -p
cd make
cmake -DBUILD_TESTS=ON ..
make
cd ..