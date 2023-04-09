# CSIDH-GMP

The implementation uses [mcl](https://github.com/herumi/mcl).

If you are using Ubuntu, you need the following environment.
```
sudo apt install g++ make libgmp-dev
```

Build mcl.
```
make -j lib/libmcl.a
```

How to get profile.

```
sudo env MCL_PERF=1 perf record ./csidh
sudo perf report
```
