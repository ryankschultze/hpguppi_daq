ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS     =

AM_CPPFLAGS += -I"@HASHPIPE_INCDIR@"

# AM_CFLAGS is used for all C compiles
AM_CFLAGS = -fopenmp -ggdb -fPIC -O3 -Wall -Werror -Wno-error=cpp -fno-strict-aliasing @X86_FEATURE_CFLAGS@
AM_CFLAGS += -Wno-address-of-packed-member 
AM_CFLAGS += -Wno-unused-const-variable
AM_CFLAGS += -Wno-stringop-truncation
AM_CFLAGS += -Wno-unused-result
AM_CFLAGS += -Wno-format-truncation
AM_CFLAGS += -Wno-unused-but-set-variable

# Convenience variables to group source files
hpguppi_databuf = hpguppi_databuf.h \
		  hpguppi_databuf.c

hpguppi_support = hpguppi_atasnap.h \
		  hpguppi_atasnap.c \
		  hpguppi_params.c \
		  hpguppi_mkfeng.h \
		  hpguppi_pksuwl.h \
		  hpguppi_rawspec.h \
		  hpguppi_rawspec.c \
		  hpguppi_time.h   \
		  hpguppi_time.c   \
		  hpguppi_util.h   \
		  hpguppi_util.c   \
		  hpguppi_udp.h    \
		  hpguppi_udp.c    \
		  hpguppi_vdif.h   \
		  polyco_struct.h  \
		  psrfits.h \
		  ata_antenna_weights_binary/include/antenna_weights.h \
		  mysigproc/mysigproc_utils.h

libsla_support = slalib.h sla.c f77.h

hpguppi_threads = hpguppi_ibverbs_pkt_thread.c \
		  hpguppi_net_thread.c \
		  hpguppi_mb1_net_thread.c \
		  hpguppi_mb128ch_net_thread.c \
		  hpguppi_ata_ibv_payload_order_thread.c \
		  hpguppi_atasnap_obs_rawdisk_thread.c \
		  hpguppi_ata_fildisk_thread.c \
		  hpguppi_rawdisk_thread.c \
		  hpguppi_rawdisk_only_thread.c \
		  hpguppi_fildisk_only_thread.c \
		  ata_antenna_weights_binary/src/antenna_weights.c \
			mysigproc/mysigproc_utils.c

# This is the hpguppi_daq plugin
lib_LTLIBRARIES = hpguppi_daq.la
hpguppi_daq_la_SOURCES  = $(hpguppi_databuf) $(hpguppi_support) $(libsla_support) $(hpguppi_threads)
hpguppi_daq_la_LIBADD   = -lsla -lrt -lrawspec -lhashpipe_ibverbs
hpguppi_daq_la_CFLAGS 	= $(AM_CFLAGS)
hpguppi_daq_la_LDFLAGS  = -avoid-version -module -shared -export-dynamic
hpguppi_daq_la_LDFLAGS += -L"@LIBSLADIR@" -Wl,-rpath,"@LIBSLADIR@"
hpguppi_daq_la_LDFLAGS += -L"@HASHPIPE_LIBDIR@" -Wl,-rpath,"@HASHPIPE_LIBDIR@"
#hpguppi_daq_la_CFLAGS  += -I./mysigproc

if BLADE_ENABLED
# bladeapi library
AM_CPPFLAGS += -I"@BLADEINCDIR@"
AM_CPPFLAGS += -I"ata_antenna_weights_binary/include"

hpguppi_blade_c_api_files = \
		hpguppi_blade_capi.h \
		hpguppi_blade_capi.cpp \
		hpguppi_blade_ata_structs.h \
		hpguppi_blade_ata_mode.h \
		hpguppi_blade_ata_mode_a_capi.cpp \
		hpguppi_blade_ata_mode_a_capi.h \
		hpguppi_blade_ata_mode_a_config.h \
		hpguppi_blade_ata_mode_b_capi.cpp \
		hpguppi_blade_ata_mode_b_capi.h \
		hpguppi_blade_ata_mode_b_config.h

noinst_LTLIBRARIES = libbladeapi.la
libbladeapi_la_SOURCES = $(hpguppi_blade_c_api_files)
libbladeapi_la_CXXFLAGS = -std=c++20 -I "@BLADEINCDIR@"
libbladeapi_la_LDFLAGS = -L"@LIBBLADEDIR@"
libbladeapi_la_LIBADD = -lblade
if CUDA_ENABLED
libbladeapi_la_CXXFLAGS  += -I "@CUDAINCDIR@"
endif

# hpguppi extension
hpguppi_support += \
		hpguppi_blade_ata_mode_b_capi.h \
		hpguppi_blade_ata_mode_b_config.h
hpguppi_databuf += \
		hpguppi_blade_databuf.c \
		hpguppi_blade_databuf.h

hpguppi_threads += \
		hpguppi_ata_blade_beamformer_thread.c \
		hpguppi_ata_blade_obs_rawdisk_thread.c

hpguppi_daq_la_LIBADD  += -lblade libbladeapi.la
hpguppi_daq_la_LDFLAGS += -L"@LIBBLADEDIR@" -Wl,-rpath,"@LIBBLADEDIR@"
hpguppi_daq_la_CFLAGS  += -I "@BLADEINCDIR@"
endif

if CUDA_ENABLED
hpguppi_daq_la_CFLAGS  += -I "@CUDAINCDIR@"
endif

if XGPU_ENABLED
hpguppi_databuf += \
		hpguppi_xgpu_databuf.c \
		hpguppi_xgpu_databuf.h

hpguppi_threads += \
		hpguppi_ata_xgpu_thread.c
if UVH5C99_ENABLED
hpguppi_threads += \
		hpguppi_ata_xintegration_thread.c
endif

hpguppi_daq_la_LDFLAGS += -L"@XGPU_LIBDIR@" -lxgpu -Xlinker -lrt
hpguppi_daq_la_CFLAGS  += -I "@XGPU_INCDIR@"
endif

if UVH5C99_ENABLED
hpguppi_threads += \
		hpguppi_ata_obs_uvh5disk_thread.c
hpguppi_daq_la_LDFLAGS += -L"@UVH5C99_LIBDIR@"
hpguppi_daq_la_LIBADD += -luvh5
hpguppi_daq_la_CFLAGS  += -I "@UVH5C99_INCDIR@"
hpguppi_daq_la_CFLAGS  += -I "@UVH5C99_INCDIR@/../external_src/tomlc99"
endif # UVH5C99_ENABLED 

if RADIOINTERFEROMETRYC99_ENABLED
hpguppi_daq_la_LDFLAGS += -L"@RADIOINTERFEROMETRYC99_LIBDIR@"
hpguppi_daq_la_LIBADD += -lradiointerferometryc99
hpguppi_daq_la_CFLAGS  += -I "@RADIOINTERFEROMETRYC99_INCDIR@"
endif # RADIOINTERFEROMETRYC99_ENABLED 

# Installed scripts
dist_bin_SCRIPTS = hpguppi_gbtstatus_loop.rb \
		   hpguppi_init.sh           \
		   init_atasnap_hashpipe.sh           \
		   kill_hpguppi.sh           \
		   start_hpguppi.sh

check_PROGRAMS = test_vdifhdr

test_vdifhdr_SOURCES = test_vdifhdr.c hpguppi_pksuwl.h hpguppi_vdif.h

# vi: set ts=8 noet :
