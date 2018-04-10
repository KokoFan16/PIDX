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

/**
 * \file PIDX_rst.c
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_raw_rst.h
 *
 */

#include "../../PIDX_inc.h"

PIDX_return_code PIDX_raw_rst_buf_create(PIDX_raw_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  int j = 0, v = 0;
  int cnt = 0, i = 0;

  for (v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    cnt = 0;
    for (i = 0; i < rst_id->reg_raw_grp_count; i++)
    {
      if (rst_id->idx_c->grank == rst_id->reg_raw_grp[i]->max_patch_rank)
      {
        PIDX_super_patch patch_group = var->rst_patch_group[cnt]; // here use patch_group

        for(j = 0; j < rst_id->reg_raw_grp[i]->patch_count; j++)
        {
          patch_group->patch[j]->buffer = malloc(patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8);

          if (patch_group->patch[j]->buffer == NULL)
            return PIDX_err_rst;

          memset(patch_group->patch[j]->buffer, 0, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->vps * var->bpv/8));
        }
        cnt++;
      }
    }
    if (cnt != var->patch_group_count)
      return PIDX_err_rst;
  }

  return PIDX_success;
}



PIDX_return_code PIDX_raw_rst_buf_destroy(PIDX_raw_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  int i, j, v;

  for(v = rst_id->first_index; v <= rst_id->last_index; v++)
  {
    PIDX_variable var = var_grp->variable[v];
    for(i = 0; i < var_grp->variable[v]->patch_group_count; i++)
    {
      for(j = 0; j < var_grp->variable[v]->rst_patch_group[i]->patch_count; j++)
      {
        free(var->rst_patch_group[i]->patch[j]->buffer);
        var->rst_patch_group[i]->patch[j]->buffer = 0;
      }
    }
  }
  return PIDX_success;
}


PIDX_return_code PIDX_raw_rst_aggregate_buf_create(PIDX_raw_rst_id rst_id)
{
  int v = 0;

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    //int bytes_per_value = var->bpv / 8;

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      // copy the size and offset to output
      PIDX_patch out_patch = var->rst_patch_group[g]->restructured_patch;

      size_t nx = out_patch->size[0];
      size_t ny = out_patch->size[1];
      size_t nz = out_patch->size[2];

      var->rst_patch_group[g]->restructured_patch->buffer = malloc(nx * ny * nz * (var->bpv/8) * var->vps);
      memset(var->rst_patch_group[g]->restructured_patch->buffer, 0, nx * ny * nz * (var->bpv/8) * var->vps);

      if (var->rst_patch_group[g]->restructured_patch->buffer == NULL)
        return PIDX_err_chunk;
    }
  }

  return PIDX_success;
}



PIDX_return_code PIDX_raw_rst_aggregate_buf_destroy(PIDX_raw_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  int v = 0;
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
      free(var->rst_patch_group[g]->restructured_patch->buffer);
  }

  return PIDX_success;
}



PIDX_return_code PIDX_raw_rst_buf_aggregate(PIDX_raw_rst_id rst_id, int mode)
{
  int v = 0;

  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];
  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = var_grp->variable[v];
    //int bytes_per_value = var->bpv / 8;

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->patch_group_count; ++g)
    {
      // copy the size and offset to output
      PIDX_super_patch patch_group = var->rst_patch_group[g];
      PIDX_patch out_patch = var->rst_patch_group[g]->restructured_patch;

      size_t k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;

      size_t nx = out_patch->size[0];
      size_t ny = out_patch->size[1];

      for (r = 0; r < var->rst_patch_group[g]->patch_count; r++)
      {
        for (k1 = patch_group->patch[r]->offset[2]; k1 < patch_group->patch[r]->offset[2] + patch_group->patch[r]->size[2]; k1++)
        {
          for (j1 = patch_group->patch[r]->offset[1]; j1 < patch_group->patch[r]->offset[1] + patch_group->patch[r]->size[1]; j1++)
          {
            for (i1 = patch_group->patch[r]->offset[0]; i1 < patch_group->patch[r]->offset[0] + patch_group->patch[r]->size[0]; i1 = i1 + patch_group->patch[r]->size[0])
            {
              index = ((patch_group->patch[r]->size[0])* (patch_group->patch[r]->size[1]) * (k1 - patch_group->patch[r]->offset[2])) + ((patch_group->patch[r]->size[0]) * (j1 - patch_group->patch[r]->offset[1])) + (i1 - patch_group->patch[r]->offset[0]);
              send_o = index * var->vps * (var->bpv/8);
              send_c = (patch_group->patch[r]->size[0]);
              recv_o = (nx * ny * (k1 - out_patch->offset[2])) + (nx * (j1 - out_patch->offset[1])) + (i1 - out_patch->offset[0]);

#if !SIMULATE_IO
              if (mode == PIDX_WRITE)
                memcpy(out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), var->rst_patch_group[g]->patch[r]->buffer + send_o, send_c * var->vps * (var->bpv/8));
              else
                memcpy(var->rst_patch_group[g]->patch[r]->buffer + send_o, out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
#endif
            }
          }
        }
      }
    }
  }

  return PIDX_success;
}
