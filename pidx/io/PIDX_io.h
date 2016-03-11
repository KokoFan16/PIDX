#include "../PIDX_inc.h"
#include "./PIDX_idx_io/PIDX_idx_io.h"
#include "./PIDX_raw_io/PIDX_raw_io.h"
#include "./PIDX_partitioned_idx_io/PIDX_partitioned_idx_io.h"

#ifndef __PIDX_IO_H
#define __PIDX_IO_H

///
double PIDX_get_time();


///
void PIDX_init_timming_buffers1(PIDX_time time, int variable_count);


///
void PIDX_init_timming_buffers2(PIDX_time time, int variable_count, int layout_count);



///
void PIDX_delete_timming_buffers1(PIDX_time time);


///
void PIDX_delete_timming_buffers2(PIDX_time time, int variable_count);

#endif