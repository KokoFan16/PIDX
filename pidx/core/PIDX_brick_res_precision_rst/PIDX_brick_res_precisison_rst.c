/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2010-2019 ViSUS L.L.C., 
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

/**
 * \file PIDX_brick_res_precision_rst.c
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_brick_res_precision_rst.h
 *
 */

#include "../../PIDX_inc.h"


PIDX_brick_res_precision_rst_id PIDX_brick_res_precision_rst_init(idx_dataset idx_meta_data, idx_comm idx_c, idx_debug idx_dbg, PIDX_restructured_grid restructured_grid, int var_start_index, int var_end_index)
{
  PIDX_brick_res_precision_rst_id brick_res_precision_rst_id;
  brick_res_precision_rst_id = (PIDX_brick_res_precision_rst_id)malloc(sizeof (*brick_res_precision_rst_id));
  memset(brick_res_precision_rst_id, 0, sizeof (*brick_res_precision_rst_id));

  brick_res_precision_rst_id->idx = idx_meta_data;
  brick_res_precision_rst_id->idx_c = idx_c;
  brick_res_precision_rst_id->idx_dbg = idx_dbg;

  brick_res_precision_rst_id->first_index = var_start_index;
  brick_res_precision_rst_id->last_index = var_end_index;

  brick_res_precision_rst_id->maximum_neighbor_count = 256;

  brick_res_precision_rst_id->reg_brick_res_precision_grp_count = 0;
  brick_res_precision_rst_id->sim_max_brick_res_precision_io_restructured_super_patch_count = 0;

  brick_res_precision_rst_id->restructured_grid = restructured_grid;

  brick_res_precision_rst_id->reg_patch_size[0] = restructured_grid->patch_size[0];
  brick_res_precision_rst_id->reg_patch_size[1] = restructured_grid->patch_size[1];
  brick_res_precision_rst_id->reg_patch_size[2] = restructured_grid->patch_size[2];

  return (brick_res_precision_rst_id);
}


PIDX_return_code PIDX_brick_res_precision_rst_finalize(PIDX_brick_res_precision_rst_id brick_res_precision_rst_id)
{
  brick_res_precision_rst_id->reg_brick_res_precision_grp_count = 0;
  brick_res_precision_rst_id->sim_max_brick_res_precision_io_restructured_super_patch_count = 0;
  free(brick_res_precision_rst_id);
  brick_res_precision_rst_id = 0;

  return PIDX_success;
}