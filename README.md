# SETI ATA Hashpipe Data Acquisition

## Compilation

The `$ ./configure` step determines which threads are compiled and indeed available for deployment.

- [`--with-libsla`](https://github.com/scottransom/pyslalib) is required for MJD time calculations.

- [`--with-libblade`](https://github.com/luigifcruz/blade) will enable Blade related threads. For this `CXX=g++-10` should be set.

- [`--with-libxgpu`](https://github.com/GPU-correlators/xGPU) will enable xGPU related threads.

- [`--with-libuvh5`](https://github.com/MydonSolutions/uvh5c99) will enable UVH5 related threads, if xGPU is enabled as that, currently, is the only scenario in which UVH5 is used.

Thereafter, `$ make` compiles the threads into `$ ./.libs/hpguppi_daq.so`.

## Instantiation

The `init_hpguppi.py` script assesses the `config_hpguppi.yaml` file to achieve the instance-system named as the first positional argument. The script is typically executed with super-user privileges.