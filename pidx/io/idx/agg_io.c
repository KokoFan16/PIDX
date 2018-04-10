/*
 * Copyright (c) 2010-2018 ViSUS L.L.C., 
 * Scientific Computing and Imaging Institute of the University of Utah
 * 
 * ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
 * University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
 *  
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * For additional information about this project contact: pascucci@acm.org
 * For support: support@visus.net
 * 
 */
#include "../../PIDX_inc.h"
static int lgi = 0;

PIDX_return_code data_io(PIDX_io file, int gi, int svi, int end_index, int mode)
{
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  int j;
  int ret = 0;
  //svi = svi - lvi;

  assert(var_grp->shared_start_layout_index == 0);
  for(j = var_grp->shared_start_layout_index; j < var_grp->agg_level; j++)
  {
    Agg_buffer temp_agg = idx->agg_buffer[svi][j];
    PIDX_block_layout temp_layout = var_grp->block_layout_by_level[j];

    file->io_id[svi][j] = PIDX_file_io_init(file->idx, file->idx_d, file->idx_c, svi, svi);

    if (file->idx_dbg->debug_do_io == 1)
    {
      if (mode == PIDX_WRITE)
        ret = PIDX_file_io_blocking_write(file->io_id[svi][j], temp_agg, temp_layout, file->idx->filename_template_partition);
      else
        ret = PIDX_file_io_blocking_read(file->io_id[svi][j], temp_agg, temp_layout, file->idx->filename_template_partition);

      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
    }
  }

  return PIDX_success;
}



PIDX_return_code data_aggregate(PIDX_io file, int gi, int svi, int evi, int agg_mode, int mode )
{
  lgi = gi;
  int ret = 0;
  PIDX_variable_group var_grp = file->idx->variable_grp[gi];
  idx_dataset_derived_metadata idx = file->idx_d;

  int j;
  PIDX_time time = file->idx_d->time;
  //svi = svi - lvi;
  //evi = evi - lvi;

  assert(var_grp->shared_start_layout_index == 0);
  for (j = var_grp->shared_start_layout_index; j < var_grp->agg_level; j++)
  {
    if (agg_mode == AGG_SETUP_AND_PERFORM || agg_mode == AGG_SETUP)
    {
      //fprintf(stderr, "svi lvi j %d %d %d\n", svi, lvi, j);
      time->agg_init_start[lgi][svi][j] = PIDX_get_time();

      //file->agg_id[svi][j] = PIDX_agg_init(file->idx, file->idx_d, file->idx_c, svi + lvi, evi + lvi, lvi);
      file->agg_id[svi][j] = PIDX_agg_init(file->idx, file->idx_d, file->idx_c, svi, evi);
      idx->agg_buffer[svi][j] = malloc(sizeof(*(idx->agg_buffer[svi][j])));
      memset(idx->agg_buffer[svi][j], 0, sizeof(*(idx->agg_buffer[svi][j])));

      idx->agg_buffer[svi][j]->file_number = -1;
      idx->agg_buffer[svi][j]->var_number = -1;
      idx->agg_buffer[svi][j]->sample_number = -1;

      idx->agg_buffer[svi][j]->no_of_aggregators = 0;
      idx->agg_buffer[svi][j]->aggregator_interval = 0;
      idx->agg_buffer[svi][j]->agg_f = 1;
      time->agg_init_end[lgi][svi][j] = PIDX_get_time();

      time->agg_meta_start[lgi][svi][j] = PIDX_get_time();
      ret = PIDX_agg_meta_data_create(file->agg_id[svi][j], idx->agg_buffer[svi][j], var_grp->block_layout_by_level[j]);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->agg_meta_end[lgi][svi][j] = PIDX_get_time();

      time->agg_buf_start[lgi][svi][j] = PIDX_get_time();

      //ret = PIDX_agg_create_randomized_aggregation(file->agg_id[svi][j], idx->agg_buffer[svi][j], block_layout_by_level[j], j, svi, file_status);

      //ret = PIDX_agg_create_global_partition_localized_aggregation_buffer(file->agg_id[svi][j], idx->agg_buffer[svi][j], var_grp->block_layout_by_level[j], j);

      ret = PIDX_agg_create_local_partition_localized_aggregation_buffer(file->agg_id[svi][j], idx->agg_buffer[svi][j], var_grp->block_layout_by_level[j], j);

      //ret = PIDX_agg_buf_create_local_uniform_dist(file->agg_id[svi][j], idx->agg_buffer[svi][j], block_layout_by_level[j]);

      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->agg_buf_end[lgi][svi][j] = PIDX_get_time();
    }

    if (agg_mode == AGG_SETUP_AND_PERFORM || agg_mode == AGG_PERFORM)
    {

      if (file->idx_dbg->debug_do_agg == 1)
      {
        time->agg_start[lgi][svi][j] = PIDX_get_time();
        ret = PIDX_agg_global_and_local(file->agg_id[svi][j], idx->agg_buffer[svi][j], j, var_grp->block_layout_by_level[j], mode);
        if (ret != PIDX_success)
        {
          fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
          return PIDX_err_rst;
        }
        time->agg_end[lgi][svi][j] = PIDX_get_time();
      }

      time->agg_meta_cleanup_start[lgi][svi][j] = PIDX_get_time();
      ret = PIDX_agg_meta_data_destroy(file->agg_id[svi][j], var_grp->block_layout_by_level[j]);
      if (ret != PIDX_success)
      {
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
      }
      time->agg_meta_cleanup_end[lgi][svi][j] = PIDX_get_time();
    }
  }


  return PIDX_success;
}
