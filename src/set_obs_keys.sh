# hashpipe_check_status -k 	DATADIR  	-s 	/mnt/buf0
# hashpipe_check_status -k 	DESTIP  	-s 	10.11.1.156
# hashpipe_check_status -k 	FENCHAN  	-i 	4096
# hashpipe_check_status -k 	NANTS  		-i 	3
# hashpipe_check_status -k 	NPOL  		-i 	2
# hashpipe_check_status -k 	NBITS  		-i 	4 
# hashpipe_check_status -k 	NSTRM  		-i 	6
# hashpipe_check_status -k 	PKTNTIME  -i 	16
# hashpipe_check_status -k 	PKTNCHAN  -i 	240
# hashpipe_check_status -k 	OBSSTART  -i 	0
# hashpipe_check_status -k 	OBSSTOP  	-i 	0
hashpipe_check_status -k 	SCHAN  		-i 	2048
# hashpipe_check_status -k 	CHAN_BW  	-f 	0.25
# hashpipe_check_status -k 	OBSBW  		-d 	128.0
hashpipe_check_status -k 	PROJID  	-s 	dmptest2
# hashpipe_check_status -k 	PKTFMT  	-s 	ATASNAPV
# hashpipe_check_status -k 	OBSFREQ  	-i 	1420

hashpipe_check_status -k OBSSTART -i 1978750000
hashpipe_check_status -k OBSSTOP -i 2078747680
