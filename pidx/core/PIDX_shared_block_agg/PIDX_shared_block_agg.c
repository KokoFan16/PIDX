#include "../../PIDX_inc.h"
#include <zfp.h>

#define PIDX_MIN(a,b) (((a)<(b))?(a):(b))

#define PIDX_ACTIVE_TARGET 1

static PIDX_return_code one_sided_data_com(PIDX_shared_block_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int mode);
static int intersectNDChunk(Ndim_patch A, Ndim_patch B);

static int shared_block_count = 0;

static PIDX_return_code create_shared_block_window(PIDX_shared_block_agg_id id, Agg_buffer ab);

struct PIDX_shared_block_agg_struct
{
  MPI_Win shared_block_win;

  idx_comm idx_c;

  /// Contains all relevant IDX file info
  /// Blocks per file, samples per block, bitmask, patch, file name template and more
  idx_dataset idx;

  /// Contains all derieved IDX file info
  /// number of files, files that are ging to be populated
  idx_dataset_derived_metadata idx_d;

  int gi;
  int fi;
  int li;

  int ***agg_r;
};



PIDX_shared_block_agg_id PIDX_shared_block_agg_init(idx_dataset idx_meta_data, idx_dataset_derived_metadata idx_d, idx_comm idx_c, int fi, int li)
{
  PIDX_shared_block_agg_id id;

  id = malloc(sizeof (*id));
  memset(id, 0, sizeof (*id));

  id->idx = idx_meta_data;
  id->idx_d = idx_d;
  id->idx_c = idx_c;

  id->gi = 0;
  id->fi = fi;
  id->li = li;

  return id;
}



PIDX_return_code PIDX_shared_block_agg_global_and_local(PIDX_shared_block_agg_id id, Agg_buffer ab, PIDX_block_layout lbl,  int MODE)
{
  if (create_shared_block_window(id, ab) != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

#ifdef PIDX_ACTIVE_TARGET
  if (MPI_Win_fence(0, id->shared_block_win) != MPI_SUCCESS)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif

  if (one_sided_data_com(id, ab, lbl, MODE) != PIDX_success)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }

#if PIDX_HAVE_MPI
#ifdef PIDX_ACTIVE_TARGET
  if (MPI_Win_fence(0, id->shared_block_win) != MPI_SUCCESS)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif
  if (MPI_Win_free(&(id->shared_block_win)) != MPI_SUCCESS)
  {
    fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
    return PIDX_err_agg;
  }
#endif

  return PIDX_success;
}


PIDX_return_code PIDX_shared_block_agg_buf_create(PIDX_shared_block_agg_id agg_id, Agg_buffer ab)
{
  PIDX_variable_group var_grp = agg_id->idx->variable_grp[agg_id->gi];
  PIDX_variable var = var_grp->variable[ab->var_number];

  if (ab->buffer_size != 0 && ab->file_number == 0)
  {
    int tcs = agg_id->idx->chunk_size[0] * agg_id->idx->chunk_size[1] * agg_id->idx->chunk_size[2];
    int bpdt = tcs * (var->bpv/8) / (agg_id->idx->compression_factor);

    shared_block_count = pow(2, agg_id->idx_d->shared_block_level - 1) / agg_id->idx_d->samples_per_block;
    agg_id->idx_d->shared_block_agg_buffer = malloc(shared_block_count * agg_id->idx_d->samples_per_block * bpdt);
    memset(agg_id->idx_d->shared_block_agg_buffer, 0, shared_block_count * agg_id->idx_d->samples_per_block * bpdt);
  }
  return PIDX_success;
}



PIDX_return_code PIDX_shared_block_agg_buf_destroy(PIDX_shared_block_agg_id agg_id, Agg_buffer ab)
{
  if (ab->buffer_size != 0 && ab->file_number == 0)
    free(agg_id->idx_d->shared_block_agg_buffer);

  return PIDX_success;
}


PIDX_return_code PIDX_shared_block_agg_finalize(PIDX_shared_block_agg_id id)
{
  free(id);
  id = 0;

  return PIDX_success;
}



static PIDX_return_code create_shared_block_window(PIDX_shared_block_agg_id id, Agg_buffer ab)
{
  int ret = 0;
  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  PIDX_variable var = var_grp->variable[ab->var_number];

  if (ab->buffer_size != 0 && ab->file_number == 0)
  {
    int tcs = id->idx->chunk_size[0] * id->idx->chunk_size[1] * id->idx->chunk_size[2];
    int bpdt = tcs * (var->bpv/8) / (id->idx->compression_factor);

    ret = MPI_Win_create(id->idx_d->shared_block_agg_buffer, (shared_block_count * id->idx_d->samples_per_block * bpdt), bpdt, MPI_INFO_NULL, id->idx_c->global_comm, &(id->shared_block_win));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }
  else
  {
    ret = MPI_Win_create(0, 0, 1, MPI_INFO_NULL, id->idx_c->global_comm, &(id->shared_block_win));
    if (ret != MPI_SUCCESS)
    {
      fprintf(stdout,"File %s Line %d\n", __FILE__, __LINE__);
      return PIDX_err_agg;
    }
  }

  return PIDX_success;
}




static PIDX_return_code one_sided_data_com(PIDX_shared_block_agg_id id, Agg_buffer ab, PIDX_block_layout lbl, int mode)
{
  int i, p, v, ret = 0;
  unsigned long long index = 0, count = 0;

  PIDX_variable_group var_grp = id->idx->variable_grp[id->gi];
  PIDX_variable var0 = var_grp->variable[id->fi];

  for(v = id->fi; v <= id->li; v++)
  {
    PIDX_variable var = var_grp->variable[v];
  }

  return PIDX_success;
}