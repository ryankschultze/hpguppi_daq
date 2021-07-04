#!/home/sonata/miniconda3/bin/python
import argparse
import os, sys
import yaml
import subprocess
import re
import socket
import datetime

default_prefix_exec =	"/opt/mnt/bin/"
default_prefix_lib 	=	"/opt/mnt/lib/"
default_logdir			= None # Switches off logging

############ Hands off from here on
parser = argparse.ArgumentParser(description='Starts an instance of Hpguppi_daq')
parser.add_argument('system', type=str,
										help='The name of the system the instance will be configured for')
parser.add_argument('instance', type=int,
										help='The instance enumeration and ID')
parser.add_argument('additional_arguments', type=str, nargs='*',
										help='Arguments added to the instance command')
parser.add_argument('-o', '--options', type=str, nargs='*', default=[],
										help='Additional key=val pairs for the Hpguppi instance\'s status buffer')
parser.add_argument('-e', '--environment-keys', type=str, nargs='*', default=[],
										help='Additional key=val pairs for the environment of the Hpguppi instance')
parser.add_argument('-d', '--dry-run', action='store_true',
										help='Do not start the instance')
parser.add_argument('--configfile', type=str, default='config_hpguppi.yml',
										help='The file containing systems\' configurations.')
args = parser.parse_args()

# Load configuration file 
with open(args.configfile, 'r') as fio:
	config = yaml.load(fio, Loader=yaml.SafeLoader)

subsystem_split_index = args.system.find('_')

system_full_name = args.system
# Select configuration for the system
if subsystem_split_index > -1 and args.system not in config:
	subsystem = args.system[subsystem_split_index+1:]
	args.system = args.system[0:subsystem_split_index]
	assert args.system in config , 'Root-system {} not defined in {}'.format(args.system, args.configfile)
	print('Accessing root-system \'{}\''.format(args.system))
else:
	subsystem = ''
	assert args.system in config , 'System {} not defined in {}'.format(args.system, args.configfile)

system = config[args.system]

# Gather cpu core count
cores_per_cpu = subprocess.run('grep cpu.cores /proc/cpuinfo'.split(' '), capture_output=True).stdout
try:
	cores_per_cpu = cores_per_cpu.decode().strip()
	m = re.match(r'cpu cores\s*:\s*(\d+).*', cores_per_cpu)
	cores_per_cpu = int(m.group(1))
except:
	print('Error trying to get the cpu core count.')
	exit(0)

# Gather system configuration for the cpu_core_count
assert 'cpu_core_count_config' in system, '{} not defined for system {} in {}'.format('cpu_core_count_config', args.system, args.configfile)
if cores_per_cpu not in system['cpu_core_count_config']:
	print('{}[{}] not defined for system {} in {}'.format('cpu_core_count_config', cores_per_cpu, args.system, args.configfile))
	exit(1)
system_ccc_config = system['cpu_core_count_config'][cores_per_cpu] if isinstance(system['cpu_core_count_config'], dict) else True

# Gather numanode_bind specifications
instance_numanode_bind = args.instance
if 'instance_numanode_bind' in system:
	instance_numanode_bind = system['instance_numanode_bind'][args.instance]
else:
	print('{} not found for system {} in {}, numactl binding matches instance enumeration'.format('instance_numanode_bind', args.system, args.configfile))

# Set cpu_core to first core
cpu_core = cores_per_cpu*instance_numanode_bind
if isinstance(system_ccc_config, dict) and 'instance_cpu_core_0' in system_ccc_config:
	if args.instance >= len(system_ccc_config['instance_cpu_core_0']):
		print('{} only defines {} instances for system {} ({} core) in {}'.format('instance_cpu_core_0', len(system_ccc_config['instance_cpu_core_0']), args.system, cores_per_cpu, args.configfile))
		exit(1)
	
	cpu_core = system_ccc_config['instance_cpu_core_0'][args.instance]
print('Fallback thread-masks start from CPU core {}.'.format(cpu_core))

# Gather the list of threads for the instance
if subsystem != '':
	assert 'subsystem_threads' in system, '\'subsystem_threads\' not defined for root-system {} in {}'.format(args.system, args.configfile)
	assert subsystem in system['subsystem_threads'], '\'{}\' not defined in the subsystem_threads for root-system {} in {}'.format(subsystem, args.system, args.configfile)
	threads = system['subsystem_threads'][subsystem]
else:
	assert 'threads' in system, '\'threads\' not defined for system {} in {}'.format(args.system, args.configfile)
	threads = system['threads']

# Gather dictionaries of thread masks' lengths and user-specified thread masks
thread_mask_length_dict = system['thread_mask_lengths'] if 'thread_mask_lengths' in system else {}
thread_instance_mask_dict = None
thread_instance_mask_dict_name = None
if isinstance(system_ccc_config, dict):
	if subsystem == '' and 'instance_thread_masks' in system_ccc_config:
		thread_instance_mask_dict = system_ccc_config['instance_thread_masks']
		thread_instance_mask_dict_name = 'instance_thread_masks'
	elif 'instance_subsystem_thread_masks' in system_ccc_config and subsystem in system_ccc_config['instance_subsystem_thread_masks']:
		thread_instance_mask_dict = system_ccc_config['instance_subsystem_thread_masks'][subsystem]
		thread_instance_mask_dict_name = 'instance_subsystem_thread_masks'

# Create the thread-list segment of the instance command
hpguppi_threads_cmd_segment = []
for thread in threads:
	if thread_instance_mask_dict is not None and thread in thread_instance_mask_dict: # explicit core masks
		assert isinstance(thread_instance_mask_dict[thread], list), '{}[{}] must define a mask for each instance as a list for system {} ({} core) in {}.'.format(
																			thread_instance_mask_dict_name, thread, system_full_name, cores_per_cpu, args.configfile
																		)
		assert args.instance < len(thread_instance_mask_dict[thread]), '{}[{}] doesn\'t define a mask for instance {} for system {} ({} core) in {}.'.format(
																			thread_instance_mask_dict_name, thread, args.instance, system_full_name, cores_per_cpu, args.configfile
																		)
		thread_mask = thread_instance_mask_dict[thread][args.instance]
		if isinstance(thread_mask, int):
			hpguppi_threads_cmd_segment.append('-c {} {}'.format(thread_mask, thread))
			cpu_core = thread_mask + 1
		elif isinstance(thread_mask, list):
			mask_val = 0
			for core_idx in thread_mask:
				mask_val += 2**core_idx
				cpu_core = core_idx + 1
			hpguppi_threads_cmd_segment.append('-m {} {}'.format(mask_val, thread))
	else: # sequential core masks
		if thread_instance_mask_dict is not None:
			print('Fallback masking thread {} from core {}'.format(thread, cpu_core))
		if thread not in thread_mask_length_dict:
			hpguppi_threads_cmd_segment.append('-c {} {}'.format(cpu_core, thread))
			cpu_core += 1
		else:
			mask_val = 0
			for core_idx in range(cpu_core, cpu_core+thread_mask_length_dict[thread]):
				mask_val += 2**core_idx
			hpguppi_threads_cmd_segment.append('-m {} {}'.format(mask_val, thread))
			cpu_core += thread_mask_length_dict[thread]

# Gather instance-agnostic instantiation variables
prefix_exec = system['prefix_exec'] if 'prefix_exec' in system else default_prefix_exec
prefix_lib = system['prefix_lib'] if 'prefix_lib' in system else default_prefix_lib
logdir = system['logdir'] if 'logdir' in system else default_logdir
# Gather optional instance-agnostic instantiation variables
command_prefix = system['command_prefix'] if 'command_prefix' in system else []
hpguppi_plugin = system['hpguppi_plugin'] if 'hpguppi_plugin' in system else 'hpguppi_daq.so'

# Gather required instance-sensitive instantiation variables
assert 'instance_datadir' in system, '{} for system {} ({} core) in {}'.format('instance_datadir', args.system, cpu_core_count, args.configfile)
instance_datadir = system['instance_datadir'][args.instance]
assert os.path.exists(instance_datadir), '{} datadir path does not exist for instance {} of system {} ({} core) in {}'.format(
	instance_datadir, args.instance, args.system, cpu_core_count, args.configfile)

assert 'instance_bindhost' in system, '{} for system {} ({} core) in {}'.format('instance_bindhost', args.system, cpu_core_count, args.configfile)
instance_bindhost = system['instance_bindhost'][args.instance]

# Gather optional instance-sensitive instantiation variables
instance_port_bind = None
if 'instance_port_bind' in system:
	instance_port_bind = system['instance_port_bind'][args.instance]
else:
	print('{} not found for system {} in {}, no per-instance BINDPORT option set'.format('instance_port_bind', args.system, args.configfile))

# Build options (key=value pairs for the status buffer)
options = [
	'DATADIR={}'.format(instance_datadir),
	'BINDHOST={}'.format(instance_bindhost),
]
if instance_port_bind is not None:
	options.append('BINDPORT={}'.format(instance_port_bind))

if 'options' in system:
	options.extend(system['options'])

# Print empty line to conclude setup and assumption prints
print()

# Build hpguppi_daq command
cmd = [
	'numactl --cpunodebind={} --membind={}'.format(instance_numanode_bind, instance_numanode_bind),
	'{}'.format(' '.join(command_prefix)),
	'{}hashpipe -p {} -I {}'.format(prefix_exec, os.path.join(prefix_lib, hpguppi_plugin), args.instance),
	' '.join(['-o {}'.format(opt) for opt in options]),
	' '.join(args.additional_arguments),
	' '.join(hpguppi_threads_cmd_segment),
]

cmd = [seg for seg in cmd if seg != '']

# Kill previous instances
previous_instance_cmd_pattern = 'hashpipe -p .* -I {}'.format(args.instance)
previous_instance_pid = subprocess.run(['pgrep', '-f', previous_instance_cmd_pattern], capture_output=True).stdout
if previous_instance_pid:
	print('Killing previous {} instance #{} ({})'.format(args.system, args.instance, previous_instance_pid.decode().strip()))
	if not args.dry_run:
		print(subprocess.run(['pkill', '-ef', previous_instance_cmd_pattern], capture_output=True).stdout.decode())

# Handle logging true switch
out_logpath = None
err_logpath = None
if logdir is not None:
	# Generate log filepaths
	hostname = socket.gethostname()
	out_logpath = os.path.join(logdir, '{}.{}.out'.format(hostname, args.instance))
	err_logpath = os.path.join(logdir, '{}.{}.err'.format(hostname, args.instance))

	# Trim existing logs
	for logpath in [out_logpath, err_logpath]:
		if os.path.exists(logpath):
			print('Trimming log: {}'.format(logpath))

			with open('tmp.out', 'w') as tmpio:
				subprocess.run('tail -n 100000 {}'.format(logpath).split(' '), stdout=tmpio)
			
			subprocess.run('mv tmp.out {}'.format(logpath).split(' '))
			with open(logpath, 'a') as logio:
				logio.write('\n{}\nStartup {}\n{}\n{}\n'.format('-'*20, datetime.datetime.now(), ' '.join(cmd), 'v'*20))
	print()

# Setup environment
environment_keys = args.environment_keys
if 'hashpipe_keyfile' in system:
	environment_keys.append('HASHPIPE_KEYFILE={}'.format(system['hashpipe_keyfile']))
if 'environment' in system:
	environment_keys.extend(system['environment'])

hashpipe_env = os.environ
for env_kv in environment_keys:
	env_kv_parts = env_kv.split('=')
	key = env_kv_parts[0]
	val = env_kv_parts[1]

	if '$'+key in val:
		replacement = hashpipe_env[key] if key in hashpipe_env else ''
		val = val.replace('$'+key, replacement)
	
	hashpipe_env[key] = val

if 'setup_commands' in system:
	setup_command_variable_dict = {
		'BINDHOST': instance_bindhost,
		'INSTANCE': args.instance,
	}

	for setup_command in system['setup_commands']:
		for var,val in setup_command_variable_dict.items():
			setup_command = setup_command.replace('${}'.format(var), str(val))
		print('#', setup_command)
		if not args.dry_run:
			print(subprocess.run(setup_command.split(' '), env=hashpipe_env, capture_output=True).stdout.decode())
	print()

cmd = ' '.join(cmd)
print(cmd)

out_logio = None if out_logpath is None else open(out_logpath, 'a')
err_logio = None if err_logpath is None else open(err_logpath, 'a')
if not args.dry_run:
	subprocess.Popen(cmd.split(' '), env=hashpipe_env, stdout=out_logio, stderr=err_logio)
else:
	print('^^^ Dry run ^^^')
	err_logio.write('Dry run')
	err_logio.write('Dry run')
	out_logio.close()
	out_logio.close()
