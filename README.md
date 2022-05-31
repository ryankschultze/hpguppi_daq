# SETI ATA Hashpipe Data Acquisition

This uses this version of [hashpipe](https://github.com/david-macmahon/hashpipe/tree/44b432af3c88ff6ccd224f526ac2b06674af04a8).
## Compilation

The `$ ./configure` step determines which threads are compiled and indeed available for deployment.

- [`--with-libsla`](https://github.com/scottransom/pyslalib) is required for MJD time calculations.

- [`--with-libblade`](https://github.com/luigifcruz/blade) will enable BLADE related threads. For this `CXX=g++-10` should be set.

- `--with-cuda-include` required to compile BLADE related C-API.

- [`--with-libxgpu`](https://github.com/GPU-correlators/xGPU) will enable xGPU related threads.

- [`--with-libuvh5`](https://github.com/MydonSolutions/uvh5c99) will enable UVH5 related threads, if xGPU is enabled as that, currently, is the only scenario in which UVH5 is used.

Thereafter, `$ make` compiles the threads into `$ ./.libs/hpguppi_daq.so`.

## Instantiation

The `init_hpguppi.py` script assesses the `config_hpguppi.yaml` file to achieve the instance-system named as the first positional argument. The script is typically executed with super-user privileges.