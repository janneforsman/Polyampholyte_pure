#!/bin/sh
time ./t.run > h26 2>&1
cp fcdfil fcd.h26
./dh.run
time ./t.run > h25 2>&1
./dh.run
time ./t.run > h24 2>&1
./dh.run
time ./t.run > h23 2>&1
./dh.run
time ./t.run > h22 2>&1
./dh.run
time ./t.run > h21 2>&1
./dh.run
time ./t.run > h20 2>&1
./dh.run
time ./t.run > h19 2>&1
./dh.run
time ./t.run > h18 2>&1
./dh.run
time ./t.run > h17 2>&1
./dh.run
time ./t.run > h16 2>&1


