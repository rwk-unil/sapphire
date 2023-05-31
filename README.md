# PP

## Build instructions

```shell
# Note all this is done in the Dockerfile, this is only for manual testing

git submodule update --init --recursive
cd xSqueezeit

# Clone and build htslib (if you already have htslib set Makefile accordingly and skip)
cd htslib
autoheader
autoconf
automake --add-missing 2>/dev/null
./configure
make
cd ..

# Clone and build zstd (if you already have zstd set Makefile accordingly and skip)
git clone https://github.com/facebook/zstd.git
cd zstd
make
make install
ldconfig
cd ..
cd ..
make
```