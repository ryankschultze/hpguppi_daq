#include <memory>

#include "blade/base.hh"
#include "blade/module.hh"
#include "blade/types.hh"

extern "C" {
#include "hpguppi_blade_capi.h"
}

using namespace Blade;

bool blade_use_device(int device_id) {
    return SetCudaDevice(device_id) == Result::SUCCESS;
}

bool blade_pin_memory(void* buffer, size_t size) {
    return Memory::PageLock(Vector<Device::CPU, I8>(buffer, size)) == Result::SUCCESS;
}