#include <stdint.h>
#include <stddef.h>

typedef void* module_t;

// Initializes all modules from pipeline.
//
// Arguments
// ---------
// batch_size : size_t
//      specifies the number of parallel workers (usually higher than two)
//
// Return
// ------
// module_t : pointer to the internal state
//
module_t blade_init(size_t batch_size);

// Destroys the module created by init().
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
//
void blade_deinit(module_t mod);

// Get the expected size of each input buffer.
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
//
size_t get_input_size(module_t mod);

// Get the expected NANTS dimension of each input buffer.
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
//
size_t get_input_dim_NANTS(module_t mod);

// Get the expected NCHANS dimension of each input buffer.
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
//
size_t get_input_dim_NCHANS(module_t mod);

// Get the expected NTIME dimension of each input buffer.
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
//
size_t get_input_dim_NTIME(module_t mod);

// Get the expected NPOLS dimension of each input buffer.
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
//
size_t get_input_dim_NPOLS(module_t mod);

// Get the expected size of each output buffer.
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
//
size_t get_output_size(module_t mod);

// Get the number of beams formed in the output buffer.
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
//
size_t get_output_NBEAMS(module_t mod);

// Process the data.
//
// Arguments
// ---------
// mod : module_t
//      pointer to the internal state
// input : void**
//      array of input buffers of size of batch_size (complex CI8)
// output : void**
//      array of output buffers of size of batch_size (complex CF16)
//
// Return
// ------
// int : error indicator (zero indicate success)
//
int blade_process(module_t mod, void** input, void** output);
