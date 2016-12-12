/*****************************************************
 **  PIDX Parallel I/O Library                      **
 **  Copyright (c) 2010-2014 University of Utah     **
 **  Scientific Computing and Imaging Institute     **
 **  72 S Central Campus Drive, Room 3750           **
 **  Salt Lake City, UT 84112                       **
 **                                                 **
 **  PIDX is licensed under the Creative Commons    **
 **  Attribution-NonCommercial-NoDerivatives 4.0    **
 **  International License. See LICENSE.md.         **
 **                                                 **
 **  For information about this project see:        **
 **  http://www.cedmav.com/pidx                     **
 **  or contact: pascucci@sci.utah.edu              **
 **  For support: PIDX-support@visus.net            **
 **                                                 **
 *****************************************************/

/**
 * \file PIDX_multi_patch_rst.c
 *
 * \author Steve Petruzza
 * \author Sidharth Kumar
 * \date   10/09/14
 *
 * Implementation of all the functions
 * declared in PIDX_multi_patch_rst.h
 *
 */

#include "../../PIDX_inc.h"
static PIDX_return_code dump_debug_data_init(PIDX_multi_patch_rst_id rst_id);
static PIDX_return_code dump_debug_data_finalie (PIDX_multi_patch_rst_id id);

PIDX_return_code PIDX_multi_patch_rst_staged_write(PIDX_multi_patch_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  if (dump_debug_data_init(rst_id) != PIDX_success)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return PIDX_err_mpi;
  }

  if (rst_id->idx->enable_rst != 1)
  {
    int v = 0, j = 0, p = 0;
    for (v = rst_id->first_index; v <= rst_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      for (p = 0; p < var->patch_group_count; p++)
      {
        Ndim_patch_group patch_group = var->rst_patch_group[p];
        for(j = 0; j < patch_group->count; j++)
          memcpy(patch_group->patch[j]->buffer, var->sim_patch[p]->buffer, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->bpv/8 * var->vps));
      }
    }
    return PIDX_success;
  }

#if PIDX_HAVE_MPI
  unsigned long long k1 = 0, i1 = 0, j1 = 0;
  unsigned long long i, j, v, index, count1 = 0, req_count = 0;
  int *send_count, *send_offset;
  unsigned long long send_c = 0, send_o = 0, counter = 0, req_counter = 0, chunk_counter = 0;
  int ret = 0;
  int pipe_length = 0;

  MPI_Request *req;
  MPI_Status *status;
  MPI_Datatype *chunk_data_type;

  //creating ample requests and statuses
  for (i = 0; i < rst_id->reg_multi_patch_grp_count; i++)
    for(j = 0; j < rst_id->reg_multi_patch_grp[i]->count; j++)
      req_count++;

  int end_index = 0;
  int start_index = 0;
  for (start_index = rst_id->first_index; start_index < (rst_id->last_index + 1); start_index = start_index + pipe_length + 1)
  {
    send_c = 0, send_o = 0, counter = 0, req_counter = 0, chunk_counter = 0;
    end_index = ((start_index + pipe_length) >= (rst_id->last_index + 1)) ? (rst_id->last_index) : (start_index + pipe_length);

    req = malloc(sizeof (*req) * req_count * 2 * (end_index - start_index + 1));
    if (!req)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return (-1);
    }
    memset(req, 0, sizeof (*req) * req_count * 2 * (end_index - start_index + 1));

    status = malloc(sizeof (*status) * req_count * 2 * (end_index - start_index + 1));
    if (!status)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return (-1);
    }
    memset(status, 0, sizeof (*status) * req_count * 2 * (end_index - start_index + 1));

    chunk_data_type =  malloc(sizeof (*chunk_data_type) * req_count  * (end_index - start_index + 1));
    if (!chunk_data_type)
    {
      fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
      return (-1);
    }
    memset(chunk_data_type, 0, sizeof (*chunk_data_type) * req_count  * (end_index - start_index + 1));

    for (i = 0; i < rst_id->reg_multi_patch_grp_count; i++)
    {
      if (rst_id->idx_c->grank == rst_id->reg_multi_patch_grp[i]->max_patch_rank)
      {
        for(j = 0; j < rst_id->reg_multi_patch_grp[i]->count; j++)
        {
          unsigned long long *reg_patch_offset = rst_id->reg_multi_patch_grp[i]->patch[j]->offset;
          unsigned long long *reg_patch_count  = rst_id->reg_multi_patch_grp[i]->patch[j]->size;

          if(rst_id->idx_c->grank == rst_id->reg_multi_patch_grp[i]->source_patch[j].rank)
          {
            count1 = 0;
            int p_index = rst_id->reg_multi_patch_grp[i]->source_patch[j].index;
            unsigned long long *sim_patch_offset = var_grp->variable[start_index]->sim_patch[p_index]->offset;
            unsigned long long *sim_patch_count = var_grp->variable[start_index]->sim_patch[p_index]->size;

            for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
              for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                {
                  index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                          (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                          (i1 - sim_patch_offset[0]);

                  for(v = start_index; v <= end_index; v++)
                  {
                    PIDX_variable var = var_grp->variable[v];
                    send_o = index * var->vps;
                    send_c = reg_patch_count[0] * var->vps;

                    if (rst_id->idx_dbg->simulate_rst == 0)
                      memcpy(var->rst_patch_group[counter]->patch[j]->buffer + (count1 * send_c * var->bpv/8), var->sim_patch[p_index]->buffer + send_o * var->bpv/8, send_c * var->bpv/8);


                    //if (rst_id->idx_dbg->dump_rst_info == 1)
                    //{
                    //  fprintf(rst_id->idx_dbg->rst_dump_fp, "[M] [%lld] Dest offset %lld Dest size %lld Source offset %lld Source size %lld\n", v, (unsigned long long)(count1 * send_c * var->bpv/8), (unsigned long long)(send_c * var->bpv/8), (unsigned long long)(send_o * var->bpv/8), (unsigned long long)(send_c * var->bpv/8));
                    //  fflush(rst_id->idx_dbg->rst_dump_fp);
                    //}

                  }
                  count1++;
                }
          }
          else
          {
            for(v = start_index; v <= end_index; v++)
            {
              PIDX_variable var = var_grp->variable[v];
              int length = (reg_patch_count[0] * reg_patch_count[1] * reg_patch_count[2]) * var->vps * var->bpv/8;
              if (rst_id->idx_dbg->simulate_rst == 0)
              {
                ret = MPI_Irecv(var->rst_patch_group[counter]->patch[j]->buffer, length, MPI_BYTE, rst_id->reg_multi_patch_grp[i]->source_patch[j].rank, 123, rst_id->idx_c->global_comm, &req[req_counter]);
                if (ret != MPI_SUCCESS)
                {
                  fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                  return PIDX_err_mpi;
                }
              }

              if (rst_id->idx_dbg->dump_rst_info == 1)
              {
                fprintf(rst_id->idx_dbg->rst_dump_fp, "[N REC] [%lld] Dest offset 0 Dest size %d My rank %d Source rank %d\n", v, length, rst_id->idx_c->grank,  rst_id->reg_multi_patch_grp[i]->source_patch[j].rank);
                fflush(rst_id->idx_dbg->rst_dump_fp);
              }

              req_counter++;
            }
          }
        }
        counter++;
      }
      else
      {
        for(j = 0; j < rst_id->reg_multi_patch_grp[i]->count; j++)
        {
          if(rst_id->idx_c->grank == rst_id->reg_multi_patch_grp[i]->source_patch[j].rank)
          {
            for(v = start_index; v <= end_index; v++)
            {
              PIDX_variable var = var_grp->variable[v];

              unsigned long long *reg_patch_count = rst_id->reg_multi_patch_grp[i]->patch[j]->size;
              unsigned long long *reg_patch_offset = rst_id->reg_multi_patch_grp[i]->patch[j]->offset;

              send_offset = malloc(sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));
              if (!send_offset)
              {
                fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                return PIDX_err_mpi;
              }
              memset(send_offset, 0, sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));

              send_count = malloc(sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));
              if (!send_count)
              {
                fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                return PIDX_err_mpi;
              }
              memset(send_count, 0, sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));

              count1 = 0;
              int p_index =  rst_id->reg_multi_patch_grp[i]->source_patch[j].index;
              unsigned long long *sim_patch_count  = var_grp->variable[start_index]->sim_patch[p_index]->size;
              unsigned long long *sim_patch_offset = var_grp->variable[start_index]->sim_patch[p_index]->offset;

              int total_send_count = 0;
              for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
                for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                  for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                  {
                    index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                            (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                            (i1 - sim_patch_offset[0]);
                    send_offset[count1] = index * var->vps * var->bpv/8;
                    send_count[count1] = reg_patch_count[0] * var->vps * var->bpv/8;
                    total_send_count = total_send_count + send_count[count1];

                    count1++;
                  }

              MPI_Type_indexed(count1, send_count, send_offset, MPI_BYTE, &chunk_data_type[chunk_counter]);
              MPI_Type_commit(&chunk_data_type[chunk_counter]);

              if (rst_id->idx_dbg->simulate_rst == 0)
              {
                ret = MPI_Isend(var->sim_patch[p_index]->buffer, 1, chunk_data_type[chunk_counter], rst_id->reg_multi_patch_grp[i]->max_patch_rank, 123, rst_id->idx_c->global_comm, &req[req_counter]);
                if (ret != MPI_SUCCESS)
                {
                  fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
                  return PIDX_err_mpi;
                }
              }

              if (rst_id->idx_dbg->dump_rst_info == 1)
              {
                fprintf(rst_id->idx_dbg->rst_dump_fp, "[N SND] [%lld] Source offset 0 Source size %d My rank %d Dest rank %d\n", v, total_send_count, rst_id->idx_c->grank,  rst_id->reg_multi_patch_grp[i]->max_patch_rank);
                fflush(rst_id->idx_dbg->rst_dump_fp);
              }

              req_counter++;
              chunk_counter++;

              free(send_offset);
              free(send_count);
            }
          }
        }
      }
    }

    if (rst_id->idx_dbg->simulate_rst == 0)
    {
      ret = MPI_Waitall(req_counter, req, status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
        return (-1);
      }
    }

    for (i = 0; i < chunk_counter; i++)
      MPI_Type_free(&chunk_data_type[i]);
    free(chunk_data_type);
    chunk_data_type = 0;

    free(req);
    req = 0;
    free(status);
    status = 0;
    req_counter = 0;
  }

  dump_debug_data_finalie(rst_id);

  return PIDX_success;
#else
  if (rst_id->idx->enable_rst == 1)
    return PIDX_err_rst;
  else
    return PIDX_success;
#endif
}

PIDX_return_code PIDX_multi_patch_rst_write(PIDX_multi_patch_rst_id rst_id)
{
  PIDX_variable_group var_grp = rst_id->idx->variable_grp[rst_id->group_index];

  if (rst_id->idx->enable_rst != 1)
  {
    int v = 0, j = 0, p = 0;
    for (v = rst_id->first_index; v <= rst_id->last_index; v++)
    {
      PIDX_variable var = var_grp->variable[v];
      for (p = 0; p < var->patch_group_count; p++)
      {
        Ndim_patch_group patch_group = var->rst_patch_group[p]; // used patch_group
        for(j = 0; j < patch_group->count; j++)
          memcpy(patch_group->patch[j]->buffer, var->sim_patch[p]->buffer, (patch_group->patch[j]->size[0] * patch_group->patch[j]->size[1] * patch_group->patch[j]->size[2] * var->bpv/8 * var->vps));
      }
    }
    return PIDX_success;
  }

#if PIDX_HAVE_MPI
  unsigned long long k1 = 0, i1 = 0, j1 = 0;
  unsigned long long i, j, v, index, count1 = 0, req_count = 0;
  int *send_count, *send_offset;
  unsigned long long send_c = 0, send_o = 0, counter = 0, req_counter = 0;
  int ret = 0;

  MPI_Request *req;
  MPI_Status *status;

  //printf("rst_id->reg_multi_patch_grp_count = %d\n", rst_id->reg_multi_patch_grp_count);
  for (i = 0; i < rst_id->reg_multi_patch_grp_count; i++)
    for(j = 0; j < rst_id->reg_multi_patch_grp[i]->count; j++)
      req_count++;

  //creating ample requests and statuses
  req = (MPI_Request*) malloc(sizeof (*req) * req_count * 2 * (rst_id->last_index - rst_id->first_index + 1));
  if (!req)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }
  memset(req, 0, sizeof (*req) * req_count * 2 * (rst_id->last_index - rst_id->first_index + 1));

  status = (MPI_Status*) malloc(sizeof (*status) * req_count * 2 * (rst_id->last_index - rst_id->first_index + 1));
  if (!status)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }
  memset(status, 0, sizeof (*status) * req_count * 2 * (rst_id->last_index - rst_id->first_index + 1));


  for (i = 0; i < rst_id->reg_multi_patch_grp_count; i++)
  {
    if (rst_id->idx_c->grank == rst_id->reg_multi_patch_grp[i]->max_patch_rank)
    {
      //  printf("%d: reg_multi_patch_grp %d count %d\n", rank, i, rst_id->reg_multi_patch_grp[i]->count);
      for(j = 0; j < rst_id->reg_multi_patch_grp[i]->count; j++)
      {
        unsigned long long *reg_patch_offset = rst_id->reg_multi_patch_grp[i]->patch[j]->offset;
        unsigned long long *reg_patch_count  = rst_id->reg_multi_patch_grp[i]->patch[j]->size;

        int patch_count = 0;
        //printf("rank %d is source patch %d??\n", rank, rst_id->reg_multi_patch_grp[i]->source_patch[j].rank);
        if(rst_id->idx_c->grank == rst_id->reg_multi_patch_grp[i]->source_patch[j].rank)
        {
          count1 = 0;
          int p_index = rst_id->reg_multi_patch_grp[i]->source_patch[j].index;

          unsigned long long *sim_patch_offset = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->offset;
          unsigned long long *sim_patch_count = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->size;
          //unsigned long long *sim_patch_offset = var_grp->variable[rst_id->first_index]->sim_patch[patch_count]->offset;
          //unsigned long long *sim_patch_count = var_grp->variable[rst_id->first_index]->sim_patch[patch_count]->size;

          //printf("%d: copy local in %d:%d [%lld %lld %lld : %lld %lld %lld]\n", rank, p_index,j, sim_patch_offset[0],sim_patch_offset[1], sim_patch_offset[2],  sim_patch_count[0], sim_patch_count[1],sim_patch_count[2]);

          // PIDX_variable var = var_grp->variable[0];
          // memset(var->rst_patch_group[counter]->patch[p_index]->buffer,100, sim_patch_count[0]*sim_patch_count[1]*sim_patch_count[2] * var->bpv/8);

          for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
            for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
              for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
              {
                // int p_index = rst_id->reg_multi_patch_grp[i]->source_patch[j].index;

                // unsigned long long *sim_patch_offset = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->offset;
                // unsigned long long *sim_patch_count = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->size;

                // printf("%d: copy local %lld %lld %lld : %lld %lld %lld\n", rank, sim_patch_offset[0],sim_patch_offset[1], sim_patch_offset[2],  sim_patch_count[0], sim_patch_count[1],sim_patch_count[2]);

                index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                    (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                    (i1 - sim_patch_offset[0]);


                for(v = rst_id->first_index; v <= rst_id->last_index; v++)
                {
                  PIDX_variable var = var_grp->variable[v];
                  send_o = index * var->vps;
                  send_c = reg_patch_count[0] * var->vps;
                  memcpy(var->rst_patch_group[counter]->patch[j]->buffer + (count1 * send_c * var->bpv/8), var->sim_patch[p_index]->buffer + send_o * var->bpv/8, send_c * var->bpv/8);
                }
                count1++;
              }
          patch_count++;
        }
        else
        {
          for(v = rst_id->first_index; v <= rst_id->last_index; v++)
          {
            PIDX_variable var = var_grp->variable[v];
            int length = (reg_patch_count[0] * reg_patch_count[1] * reg_patch_count[2]) * var->vps * var->bpv/8;
            ret = MPI_Irecv(var->rst_patch_group[counter]->patch[j]->buffer, length, MPI_BYTE, rst_id->reg_multi_patch_grp[i]->source_patch[j].rank, 123, rst_id->idx_c->global_comm, &req[req_counter]);
            if (ret != MPI_SUCCESS)
            {
              fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
              return PIDX_err_mpi;
            }
            req_counter++;
          }
        }
      }
      counter++;
    }
    else
    {
      for(j = 0; j < rst_id->reg_multi_patch_grp[i]->count; j++)
      {
        int patch_count = 0;
        if(rst_id->idx_c->grank == rst_id->reg_multi_patch_grp[i]->source_patch[j].rank)
        {
          for(v = rst_id->first_index; v <= rst_id->last_index; v++)
          {
            PIDX_variable var = var_grp->variable[v];

            unsigned long long *reg_patch_count = rst_id->reg_multi_patch_grp[i]->patch[j]->size;
            unsigned long long *reg_patch_offset = rst_id->reg_multi_patch_grp[i]->patch[j]->offset;

            send_offset = malloc(sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));
            if (!send_offset)
            {
              fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
              return PIDX_err_mpi;
            }
            memset(send_offset, 0, sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));

            send_count = malloc(sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));
            if (!send_count)
            {
              fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
              return PIDX_err_mpi;
            }
            memset(send_count, 0, sizeof (int) * (reg_patch_count[1] * reg_patch_count[2]));

            count1 = 0;
            int p_index =  rst_id->reg_multi_patch_grp[i]->source_patch[j].index;

            unsigned long long *sim_patch_count  = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->size;
            unsigned long long *sim_patch_offset = var_grp->variable[rst_id->first_index]->sim_patch[p_index]->offset;

            for (k1 = reg_patch_offset[2]; k1 < reg_patch_offset[2] + reg_patch_count[2]; k1++)
              for (j1 = reg_patch_offset[1]; j1 < reg_patch_offset[1] + reg_patch_count[1]; j1++)
                for (i1 = reg_patch_offset[0]; i1 < reg_patch_offset[0] + reg_patch_count[0]; i1 = i1 + reg_patch_count[0])
                {
                  index = (sim_patch_count[0] * sim_patch_count[1] * (k1 - sim_patch_offset[2])) +
                      (sim_patch_count[0] * (j1 - sim_patch_offset[1])) +
                      (i1 - sim_patch_offset[0]);
                  send_offset[count1] = index * var->vps * var->bpv/8;
                  send_count[count1] = reg_patch_count[0] * var->vps * var->bpv/8;

                  count1++;
                }

            MPI_Datatype chunk_data_type;
            MPI_Type_indexed(count1, send_count, send_offset, MPI_BYTE, &chunk_data_type);
            MPI_Type_commit(&chunk_data_type);
            ret = MPI_Isend(var->sim_patch[p_index]->buffer, 1, chunk_data_type, rst_id->reg_multi_patch_grp[i]->max_patch_rank, 123, rst_id->idx_c->global_comm, &req[req_counter]);

            if (ret != MPI_SUCCESS)
            {
              fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
              return PIDX_err_mpi;
            }
            req_counter++;

            MPI_Type_free(&chunk_data_type);
            free(send_offset);
            free(send_count);
          }
        }
        patch_count++;
      }
    }
  }

#if !SIMULATE_IO

  ret = MPI_Waitall(req_counter, req, status);
  if (ret != MPI_SUCCESS)
  {
    fprintf(stderr, "Error: File [%s] Line [%d]\n", __FILE__, __LINE__);
    return (-1);
  }

#endif

  free(req);
  req = 0;
  free(status);
  status = 0;

  return PIDX_success;
#else
  if (rst_id->idx->enable_rst == 1)
    return PIDX_err_rst;
  else
    return PIDX_success;
#endif

}


static PIDX_return_code dump_debug_data_init (PIDX_multi_patch_rst_id id)
{
  int ret = 0;
  if (id->idx_dbg->dump_rst_info == 1)
  {
    char rst_file_name[1024];
    ret = mkdir(id->idx_dbg->rst_dump_dir_name, S_IRWXU | S_IRWXG | S_IRWXO);
    if (ret != 0 && errno != EEXIST)
    {
      perror("mkdir");
      fprintf(stderr, " Error in rstregate_write_read Line %d File %s folder name %s\n", __LINE__, __FILE__, id->idx_dbg->rst_dump_dir_name);
      return PIDX_err_rst;
    }

    MPI_Barrier(id->idx_c->global_comm);

    sprintf(rst_file_name, "%s/rank_%d", id->idx_dbg->rst_dump_dir_name, id->idx_c->grank);
    id->idx_dbg->rst_dump_fp = fopen(rst_file_name, "a+");
    if (!id->idx_dbg->rst_dump_fp)
    {
      fprintf(stderr, " [%s] [%d] rst_dump_fp filename = %s is corrupt.\n", __FILE__, __LINE__, rst_file_name);
      return PIDX_err_rst;
    }
  }

  return PIDX_success;
}


static PIDX_return_code dump_debug_data_finalie (PIDX_multi_patch_rst_id id)
{

  if (id->idx_dbg->dump_rst_info == 1)
  {
    fprintf(id->idx_dbg->rst_dump_fp, "\n");
    fclose(id->idx_dbg->rst_dump_fp);
  }

  return PIDX_success;
}
