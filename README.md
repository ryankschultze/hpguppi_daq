# SETI ATA Hashpipe Data Acquisition

This uses this version of [hashpipe](https://github.com/MydonSolutions/hashpipe/tree/ae6c1541a1000f921ba1e0d5aafee6be6f8a8740).
## Compilation

The `$ ./configure` step determines which threads are compiled and indeed available for deployment.

- [`--with-sla-lib`](#slahttpsgithubcomscottransompyslalib) is required for MJD time calculations. 
Provide the directory containing the library.

- [`--with-blade`](#bladehttpsgithubcomluigifcruzblade) will enable BLADE related threads. 
For this `CXX=g++-10` should be set.
Provide the installation directory of BLADE (-Dprefix defaults to /usr/local).

- `--with-cuda-include` required to compile BLADE related C-API.

- [`--with-xgpu`](https://github.com/GPU-correlators/xGPU) will enable xGPU related threads.
Provide the compilation directory (i.e. xGPU/src).

- [`--with-uvh5`](https://github.com/MydonSolutions/uvh5c99) will enable UVH5 related threads.
Provide the compilation directory (i.e. uvh5c99/build).

Thereafter, `$ make` compiles the threads into `$ ./.libs/hpguppi_daq.so`.

## Instantiation

The `init_hpguppi.py` script assesses the `config_hpguppi.yaml` file to achieve the instance-system named as the first positional argument. The script is typically executed with super-user privileges.