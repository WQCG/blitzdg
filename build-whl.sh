#!/bin/bash
set -e -x

python3 setup.py bdist_wheel sdist

export PLAT=manylinux2010_x86_64

# Bundle external shared libraries into the wheels
cp pyblitzdg/*.so .
for whl in dist/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w wheelhouse/
done
