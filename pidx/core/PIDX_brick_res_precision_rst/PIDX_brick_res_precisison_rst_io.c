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
 * \file PIDX_rst.c
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


PIDX_return_code PIDX_brick_res_precision_rst_buf_aggregate_and_write(PIDX_brick_res_precision_rst_id rst_id)
{
  int v;
  char *directory_path;
  int brick_res_precision_io_pipe_length = 0;

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  char time_template[512];
  sprintf(time_template, "%%s/%s/%%d_%%d", rst_id->idx->filename_time_template);

  int g = 0;
  PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index];
  for (g = 0; g < var0->brick_res_precision_io_restructured_super_patch_count; ++g)
  {
    //int bytes_per_value = var->bpv / 8;
    // loop through all groups
    char *file_name;
    file_name = malloc(PATH_MAX * sizeof(*file_name));
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));

    sprintf(file_name, time_template, directory_path, rst_id->idx->current_time_step, rst_id->idx_c->simulation_rank, g);
    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    int v_start = 0, v_end = 0;
    int start_var_index = rst_id->first_index;
    int end_var_index = rst_id->last_index + 1;
    for (v_start = start_var_index; v_start < end_var_index; v_start = v_start + (brick_res_precision_io_pipe_length + 1))
    {
      v_end = ((v_start + brick_res_precision_io_pipe_length) >= (end_var_index)) ? (end_var_index - 1) : (v_start + brick_res_precision_io_pipe_length);

      // copy the size and offset to output
      PIDX_variable var_start = rst_id->idx->variable[v_start];
      PIDX_super_patch patch_group = var_start->brick_res_precision_io_restructured_super_patch[g];
      PIDX_patch out_patch = var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch;

      uint64_t nx = out_patch->size[0];
      uint64_t ny = out_patch->size[1];
      uint64_t nz = out_patch->size[2];

      int bits = 0;
      for (v = v_start; v <= v_end; v++)
      {
        PIDX_variable var = rst_id->idx->variable[v];
        bits = bits + (var->bpv/8) * var->vps;
      }

      //PIDX_variable var = rst_id->idx->variable[v];
      unsigned char* reg_patch_buffer = malloc(nx * ny * nz * bits);
      memset(reg_patch_buffer, 0, nx * ny * nz * bits);
      if (reg_patch_buffer == NULL)
        return PIDX_err_chunk;

      uint64_t k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
      for (r = 0; r < var_start->brick_res_precision_io_restructured_super_patch[g]->patch_count; r++)
      {
        for (k1 = patch_group->patch[r]->offset[2]; k1 < patch_group->patch[r]->offset[2] + patch_group->patch[r]->size[2]; k1++)
        {
          for (j1 = patch_group->patch[r]->offset[1]; j1 < patch_group->patch[r]->offset[1] + patch_group->patch[r]->size[1]; j1++)
          {
            for (i1 = patch_group->patch[r]->offset[0]; i1 < patch_group->patch[r]->offset[0] + patch_group->patch[r]->size[0]; i1 = i1 + patch_group->patch[r]->size[0])
            {
              index = ((patch_group->patch[r]->size[0])* (patch_group->patch[r]->size[1]) * (k1 - patch_group->patch[r]->offset[2])) + ((patch_group->patch[r]->size[0]) * (j1 - patch_group->patch[r]->offset[1])) + (i1 - patch_group->patch[r]->offset[0]);
              send_o = index;
              send_c = (patch_group->patch[r]->size[0]);
              recv_o = (nx * ny * (k1 - out_patch->offset[2])) + (nx * (j1 - out_patch->offset[1])) + (i1 - out_patch->offset[0]);

              for (v = v_start; v <= v_end; v++)
              {
                int v1 = 0;
                uint64_t data_offset = 0;
                for (v1 = v_start; v1 < v; v1++)
                {
                  data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (rst_id->idx->variable[v1]->vps * (rst_id->idx->variable[v1]->bpv/8)));
                }
                PIDX_variable var = rst_id->idx->variable[v];
                memcpy(reg_patch_buffer + data_offset + (recv_o * var->vps * (var->bpv/8)), var->brick_res_precision_io_restructured_super_patch[g]->patch[r]->buffer + send_o * var->vps * (var->bpv/8), send_c * var->vps * (var->bpv/8));
              }
            }
          }
        }
      }

      uint64_t data_offset = 0, v1 = 0;
      for (v1 = 0; v1 < v_start; v1++)
        data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (rst_id->idx->variable[v1]->vps * (rst_id->idx->variable[v1]->bpv/8)));

      uint64_t buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
      uint64_t write_count = pwrite(fp, reg_patch_buffer, buffer_size, data_offset);
      if (write_count != buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }

      free(reg_patch_buffer);
      reg_patch_buffer = 0;
    }
    close(fp);
    free(file_name);
  }
  free(directory_path);

  return PIDX_success;
}



PIDX_return_code PIDX_brick_res_precision_rst_buf_read_and_aggregate(PIDX_brick_res_precision_rst_id rst_id)
{
  int v;
  MPI_File fh;
  char *directory_path;
  char *data_set_path;

  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);

  data_set_path = malloc(sizeof(*data_set_path) * PATH_MAX);
  memset(data_set_path, 0, sizeof(*data_set_path) * PATH_MAX);

  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);
  sprintf(data_set_path, "%s/time%09d/", directory_path, rst_id->idx->current_time_step);

  char time_template[512];
  sprintf(time_template, "%%s/%s/%%d_%%d", rst_id->idx->filename_time_template);

  for (v = rst_id->first_index; v <= rst_id->last_index; ++v)
  {
    PIDX_variable var = rst_id->idx->variable[v];

    // loop through all groups
    int g = 0;
    for (g = 0; g < var->brick_res_precision_io_restructured_super_patch_count; ++g)
    {
      // copy the size and offset to output
      PIDX_super_patch patch_group = var->brick_res_precision_io_restructured_super_patch[g];
      PIDX_patch out_patch = var->brick_res_precision_io_restructured_super_patch[g]->restructured_patch;

      uint64_t nx = out_patch->size[0];
      uint64_t ny = out_patch->size[1];
      uint64_t nz = out_patch->size[2];

      var->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer = malloc(nx * ny * nz * (var->bpv/8) * var->vps);
      memset(var->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer, 0, nx * ny * nz * (var->bpv/8) * var->vps);

      if (var->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer == NULL)
        return PIDX_err_chunk;

      uint64_t data_offset = 0, v1 = 0;
      for (v1 = 0; v1 < v; v1++)
        data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (rst_id->idx->variable[v1]->vps * (rst_id->idx->variable[v1]->bpv/8)));

      uint64_t buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (var->vps * (var->bpv/8));

      char *file_name;
      file_name = malloc(PATH_MAX * sizeof(*file_name));
      memset(file_name, 0, PATH_MAX * sizeof(*file_name));

      sprintf(file_name, time_template, directory_path, rst_id->idx->current_time_step, rst_id->idx_c->simulation_rank, g);

      MPI_Status status;
      int ret = 0;
      ret = MPI_File_open(MPI_COMM_SELF, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Line %d File %s File opening %s\n", __LINE__, __FILE__, file_name);
        return PIDX_err_rst;
      }

      ret = MPI_File_read_at(fh, data_offset, out_patch->buffer, (buffer_size), MPI_BYTE, &status);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }

      ret = MPI_File_close(&fh);
      if (ret != MPI_SUCCESS)
      {
        fprintf(stderr, "Line %d File %s\n", __LINE__, __FILE__);
        return PIDX_err_rst;
      }

      uint64_t k1, j1, i1, r, index = 0, recv_o = 0, send_o = 0, send_c = 0;
      for (r = 0; r < var->brick_res_precision_io_restructured_super_patch[g]->patch_count; r++)
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

              memcpy(var->brick_res_precision_io_restructured_super_patch[g]->patch[r]->buffer + send_o, out_patch->buffer + (recv_o * var->vps * (var->bpv/8)), send_c * var->vps * (var->bpv/8));
            }
          }
        }
      }

      free(var->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer);
    }
  }

  return PIDX_success;
}


// A wavelet helper
void PIDX_wavelet_helper(unsigned char* buf, int step, int ng_step, int flag, int bits, uint64_t x, uint64_t y, uint64_t z, char* type_name)
{
  int si = ng_step; int sj = ng_step; int sk = ng_step;

  // Define the start position based on orientations (x, y, z)
  if (flag == 0)
      sj = step;
  if (flag == 1)
      si = step;
  if (flag == 2)
      sk = step;

  int neighbor_ind = 1;

	// data types
  unsigned char c_data = 0;
  unsigned char c_neigb = 0;
  short s_data = 0;
  short s_neigb = 0;
  float f_data = 0;
  float f_neigb = 0;
  double d_data = 0;
  double d_neigb = 0;
  int i_data = 0;
  int i_neigb = 0;
  uint64_t u64i_data = 0;
  uint64_t u64i_neigb = 0;
  int64_t i64_data = 0;
  int64_t i64_neigb = 0;

  for (int k = 0; k < z; k+=sk)
  {
	for (int i = 0; i < y; i+=si)
	{
	  for (int j = 0; j < x; j+=sj)
	  {
	    int index = k * x * y + i * x + j;
	    // Define the neighbor position based on orientations (x, y, z)
	    if (flag == 0)
	      neighbor_ind = index + ng_step;
	    if (flag == 1)
	      neighbor_ind = index + ng_step * x;
	    if (flag == 2)
	      neighbor_ind = index + ng_step * y * x;


	    if (strcmp(type_name, PIDX_DType.UINT8) == 0 || strcmp(type_name, PIDX_DType.UINT8_GA) == 0 || strcmp(type_name, PIDX_DType.UINT8_RGB) == 0)
	    {
	      c_data = buf[index];
	      c_neigb = buf[neighbor_ind];
	      // Calculate wavelet coefficients and replace in the buffer
	      buf[index] = (c_data + c_neigb)/2;
	      buf[neighbor_ind] = buf[index] - c_neigb;
	    }
	    if (strcmp(type_name, PIDX_DType.INT16) == 0 || strcmp(type_name, PIDX_DType.INT16_GA) == 0 || strcmp(type_name, PIDX_DType.INT16_RGB) == 0)
	    {
	      // Covert unsigned char to short
	      memcpy(&s_data, &buf[index * sizeof(short)], sizeof(short));
	      memcpy(&s_neigb, &buf[neighbor_ind * sizeof(float)], sizeof(short));
	      // Calculate wavelet coefficients
	      short avg = (s_data + s_neigb) / 2;
	      short dif = avg - s_neigb;
	      // Replace buffer
	      memcpy(&buf[index * sizeof(short)], &avg, sizeof(short));
	      memcpy(&buf[neighbor_ind * sizeof(short)], &dif, sizeof(short));
	    }
	    if (strcmp(type_name, PIDX_DType.INT32) == 0 || strcmp(type_name, PIDX_DType.INT32_GA) == 0 || strcmp(type_name, PIDX_DType.INT32_RGB) == 0)
	    {
	      // Covert unsigned char to int
	      memcpy(&i_data, &buf[index * sizeof(int)], sizeof(int));
	      memcpy(&i_neigb, &buf[neighbor_ind * sizeof(int)], sizeof(int));
	      // Calculate wavelet coefficients
	      int avg = (i_data + i_neigb) / 2;
	      int dif = avg - i_neigb;
	      // Replace buffer
	      memcpy(&buf[index * sizeof(int)], &avg, sizeof(int));
	      memcpy(&buf[neighbor_ind * sizeof(int)], &dif, sizeof(int));
	    }
	    else if (strcmp(type_name, PIDX_DType.FLOAT32) == 0 || strcmp(type_name, PIDX_DType.FLOAT32_GA) == 0 || strcmp(type_name, PIDX_DType.FLOAT32_RGB) == 0)
	    {
	      // Covert unsigned char to float
	      memcpy(&f_data, &buf[index * sizeof(float)], sizeof(float));
	      memcpy(&f_neigb, &buf[neighbor_ind * sizeof(float)], sizeof(float));
	      // Calculate wavelet coefficients
	      float avg = (f_data + f_neigb) / 2.0;
	      float dif = avg - f_neigb;
		// Replace buffer
	      memcpy(&buf[index * sizeof(float)], &avg, sizeof(float));
	      memcpy(&buf[neighbor_ind * sizeof(float)], &dif, sizeof(float));
	    }
	    else if (strcmp(type_name, PIDX_DType.FLOAT64) == 0 || strcmp(type_name, PIDX_DType.FLOAT64_GA) == 0 || strcmp(type_name, PIDX_DType.FLOAT64_RGB) == 0)
	    {
	      // Covert unsigned char to double
	      memcpy(&d_data, &buf[index * sizeof(double)], sizeof(double));
	      memcpy(&d_neigb, &buf[neighbor_ind * sizeof(double)], sizeof(double));
	      // Calculate wavelet coefficients
	      double avg = (d_data + d_neigb) / 2.0;
	      double dif = avg - d_neigb;
	      // Replace buffer
	      memcpy(&buf[index * sizeof(double)], &avg, sizeof(double));
	      memcpy(&buf[neighbor_ind * sizeof(double)], &dif, sizeof(double));
	    }
	    else if (strcmp(type_name, PIDX_DType.INT64) == 0 || strcmp(type_name, PIDX_DType.INT64_GA) == 0 || strcmp(type_name, PIDX_DType.INT64_RGB) == 0)
	    {
	      // Covert unsigned char to int64_t
	      memcpy(&i64_data, &buf[index * sizeof(int64_t)], sizeof(int64_t));
	      memcpy(&i64_neigb, &buf[neighbor_ind * sizeof(int64_t)], sizeof(int64_t));
	      // Calculate wavelet coefficients
	      int64_t avg = (i64_data + i64_neigb) / 2.0;
	      int64_t dif = avg - i64_neigb;
	      // Replace buffer
	      memcpy(&buf[index * sizeof(int64_t)], &avg, sizeof(int64_t));
	      memcpy(&buf[neighbor_ind * sizeof(int64_t)], &dif, sizeof(int64_t));
	    }
	    else if (strcmp(type_name, PIDX_DType.UINT64) == 0 || strcmp(type_name, PIDX_DType.UINT64_GA) == 0 || strcmp(type_name, PIDX_DType.UINT64_RGB) == 0)
	    {
	      // Covert unsigned char to uint64_t
	      memcpy(&u64i_data, &buf[index * sizeof(uint64_t)], sizeof(uint64_t));
	      memcpy(&u64i_neigb, &buf[neighbor_ind * sizeof(uint64_t)], sizeof(uint64_t));
	      // Calculate wavelet coefficients
	      uint64_t avg = (u64i_data + u64i_neigb) / 2.0;
	      uint64_t dif = avg - u64i_neigb;
	      // Replace buffer
	      memcpy(&buf[index * sizeof(uint64_t)], &avg, sizeof(uint64_t));
	      memcpy(&buf[neighbor_ind * sizeof(uint64_t)], &dif, sizeof(uint64_t));
	    }
	  }
	}
  }
}

// Wavelet transform
void PIDX_wavelet_transform(unsigned char* buffer, uint64_t x, uint64_t y, uint64_t z, int bits, char* type_name, int wavelet_level)
{
  for (int level = 1; level <= wavelet_level; level++)
  {
    int step = pow(2, level);
    int ng_step = step/2;

    // Calculate x-dir
    PIDX_wavelet_helper(buffer, step, ng_step, 0, bits, x, y, z, type_name);
    // Calculate y-dir
    PIDX_wavelet_helper(buffer, step, ng_step, 1, bits, x, y, z, type_name);
    // Calculate z-dir
    PIDX_wavelet_helper(buffer, step, ng_step, 2, bits, x, y, z, type_name);
  }
}


// Structure of ZFP pointer
struct PIDX_zfp_compress_pointer
{
  unsigned char *p;
  int compress_size;
};


// ZFP compression
struct PIDX_zfp_compress_pointer PIDX_compress_3D_float(unsigned char* buf, int dim_x, int dim_y, int dim_z, int flag, float param, char* type_name)
{
  // ZFP data type according to PIDX data type
  zfp_type type = zfp_type_none;
  if (strcmp(type_name, PIDX_DType.INT32) == 0 || strcmp(type_name, PIDX_DType.INT32_GA) == 0 || strcmp(type_name, PIDX_DType.INT32_RGB) == 0)
	type = zfp_type_int32;
  else if (strcmp(type_name, PIDX_DType.FLOAT32) == 0 || strcmp(type_name, PIDX_DType.FLOAT32_GA) == 0 || strcmp(type_name, PIDX_DType.FLOAT32_RGB) == 0)
	type = zfp_type_float;
  else if (strcmp(type_name, PIDX_DType.INT64) == 0 || strcmp(type_name, PIDX_DType.INT64_GA) == 0 || strcmp(type_name, PIDX_DType.INT64_RGB) == 0)
	type = zfp_type_int64;
  else if (strcmp(type_name, PIDX_DType.FLOAT64) == 0 || strcmp(type_name, PIDX_DType.FLOAT64_GA) == 0 || strcmp(type_name, PIDX_DType.FLOAT64_RGB) == 0)
	type = zfp_type_double;
  else
	printf("ERROR: ZFP cannot handle type %s\n", type_name);

  zfp_field* field = zfp_field_3d(buf, type, dim_x, dim_y, dim_z);
  zfp_stream* zfp = zfp_stream_open(NULL);
  // Two compression modes
  if (flag == 0)
	zfp_stream_set_accuracy(zfp, param);
  else if (flag == 1)
	zfp_stream_set_precision(zfp, param);
  else
  {
	printf("ERROR: O means accuracy, and 1 means precision\n");
  }
  size_t max_compressed_bytes = zfp_stream_maximum_size(zfp, field);
  // ZFP pointer structure
  struct PIDX_zfp_compress_pointer output;
  output.p = (unsigned char*) malloc(max_compressed_bytes);
  bitstream* stream = stream_open(output.p, max_compressed_bytes);
  zfp_stream_set_bit_stream(zfp, stream);
  size_t compressed_bytes = zfp_compress(zfp, field);
  output.compress_size = compressed_bytes; // Data size after compression
  if (compressed_bytes == 0)
	puts("ERROR: Something wrong happened during compression\n");
  zfp_field_free(field);
  zfp_stream_close(zfp);
  stream_close(stream);
  return output;
}


// Calculate wavelet level dimensions
void PIDX_calculate_level_dimension(uint64_t* size, uint64_t* brick_size, int level)
{
  for (int i = 0; i < 3; i++)
  {
	size[i] = brick_size[i]/pow(2, level);
  }
}


// A reorganisation hepler
void PIDX_reorg_helper(unsigned char* buf, unsigned char* level_buf, int step, int *index, int sk, int si, int sj, uint64_t x, uint64_t y, uint64_t z, int bits)
{
  for (int k = sk; k < z; k+=step)
  {
	for (int i = si; i < y; i+=step)
	{
	  for (int j = sj; j < x; j+=step)
	  {
		int position = k * y * x + i * x + j;
		memcpy(&level_buf[(*index) * bits], &buf[position * bits], bits);
		*index += 1;
	  }
	}
  }
}

// Combine subbands of the top wavalet level when dc component is 8
void PIDX_compress_top_buffer(unsigned char* comp_buf, unsigned char* buf, uint64_t* comp_size, int bits, int wavelet_level, uint64_t* patch_size, char* type_name, int end_level, int dc_size, int comp_mode, float param)
{
  unsigned char* level_buf = NULL;

  if (dc_size < 64)
	dc_size = 64;
  level_buf = (unsigned char*) malloc(dc_size * bits);

  // Read DC component
  int step = pow(2, wavelet_level);
  int index = 0;
  PIDX_reorg_helper(buf, level_buf, step, &index, 0, 0, 0, patch_size[0], patch_size[1], patch_size[2], bits);

  // Read each subbands for top two level
  for (int level = wavelet_level; level > end_level; level--)
  {
	step = pow(2, level);
	int n_step[2] = {0, step/2};
	for(int k = 0; k < 2; k++)
	{
	  for(int i = 0; i < 2; i++)
	  {
		for(int j = 0; j < 2; j++)
		{
		  if (j == 0 && i == 0 && k == 0)
			continue;
		  PIDX_reorg_helper(buf, level_buf, step, &index, n_step[k], n_step[i], n_step[j], patch_size[0], patch_size[1], patch_size[2], bits);
		}
	  }
	}
  }

  uint64_t level_size[3] = {4, 4, 4};
  if (dc_size > 64)
	PIDX_calculate_level_dimension(level_size, patch_size, wavelet_level);

  // ZFP compress (0 means the accuracy, and followed 0 means the tolerance (need to be changed))
  struct PIDX_zfp_compress_pointer output = PIDX_compress_3D_float(level_buf, level_size[0], level_size[1], level_size[2], comp_mode, param, type_name);
  free(level_buf);
  memcpy(&comp_buf[*comp_size], output.p, output.compress_size); // Combine buffer
  *comp_size += output.compress_size;
}

// ZFP compression of each subband per level
void PIDX_compressed_subbands(unsigned char* comp_buf, unsigned char* buf, uint64_t* comp_size, int start_level, uint64_t* patch_size, int bits, char* type_name, int* comp_blocks_sizes, int comp_mode, float param)
{
  int count = 0;
  uint64_t level_size[3];

  for (int level = start_level; level > 0; level--)
  {
	// Calculate dimention per level
	PIDX_calculate_level_dimension(level_size, patch_size, level);
	int size = level_size[0] * level_size[1] * level_size[2];
	unsigned char* level_buf = (unsigned char *)malloc(size * bits);

	int step = pow(2, level);
	int n_step[2] = {0, step/2};

	// Define the start points for HHL, HLL, HLH ...
	for(int k = 0; k < 2; k++)
	{
	  for(int i = 0; i < 2; i++)
	  {
		for(int j = 0; j < 2; j++)
		{
		  if (j == 0 && i == 0 && k == 0)
			continue;
		  else
		  {
		    int index = 0;
		    // Read subbands per level
		    PIDX_reorg_helper(buf, level_buf, step, &index, n_step[k], n_step[i], n_step[j], patch_size[0], patch_size[1], patch_size[2], bits);
		    // ZFP compression per subbands of each level
		    struct PIDX_zfp_compress_pointer output = PIDX_compress_3D_float(level_buf, level_size[0], level_size[1], level_size[2], comp_mode, param, type_name);
		    comp_blocks_sizes[count] = output.compress_size;
		    memcpy(&comp_buf[*comp_size], output.p, output.compress_size); // Combine buffer
		    *comp_size += output.compress_size;
		    count++;
		  }
		}
	  }
	}
	free(level_buf);
  }
}


void PIDX_wavelet_compression(unsigned char* comp_buf, unsigned char* buf, int bits, uint64_t dc_size, uint64_t* patch_size, char* type_name, PIDX_patch out_patch, int comp_mode, float param)
{
  uint64_t comp_size = 0;
  int comp_count = 0;
  int wavelet_level = out_patch->wavelet_level;

  if (dc_size == 1)
  {
	comp_count = (wavelet_level-2) * 7 + 1;
	out_patch->compressed_blocks_sizes = (int*) malloc((comp_count - 1) * sizeof(int));
	PIDX_compress_top_buffer(comp_buf, buf, &comp_size, bits, wavelet_level, patch_size, type_name, wavelet_level-2, dc_size, comp_mode, param);
	out_patch->first_compressed_size = comp_size;
	PIDX_compressed_subbands(comp_buf, buf, &comp_size, wavelet_level-2, patch_size, bits, type_name, out_patch->compressed_blocks_sizes, comp_mode, param);
	out_patch->total_compress_size = comp_size; // Store the total compressed size for each brick
	out_patch->num_compress_blocks = comp_count; // Store the number of compressed blocks for each brick
  }
  else if (dc_size == 8)
  {
	comp_count = (wavelet_level-1) * 7 + 1;
	out_patch->compressed_blocks_sizes = (int*) malloc((comp_count - 1) * sizeof(int));
	PIDX_compress_top_buffer(comp_buf, buf, &comp_size, bits, wavelet_level, patch_size, type_name, wavelet_level-1, dc_size, comp_mode, param);
	out_patch->first_compressed_size = comp_size;
	PIDX_compressed_subbands(comp_buf, buf, &comp_size, wavelet_level-1, patch_size, bits, type_name, out_patch->compressed_blocks_sizes, comp_mode, param);
	out_patch->total_compress_size = comp_size; // Store the total compressed size for each brick
	out_patch->num_compress_blocks = comp_count; // Store the number of compressed blocks for each brick
  }
  else
  {
	comp_count = wavelet_level * 7 + 1;
	out_patch->compressed_blocks_sizes = (int*) malloc((comp_count - 1) * sizeof(int));
	PIDX_compress_top_buffer(comp_buf, buf, &comp_size, bits, wavelet_level, patch_size, type_name, wavelet_level, dc_size, comp_mode, param);
	out_patch->first_compressed_size = comp_size;
	PIDX_compressed_subbands(comp_buf, buf, &comp_size, wavelet_level, patch_size, bits, type_name, out_patch->compressed_blocks_sizes, comp_mode, param);
	out_patch->total_compress_size = comp_size; // Store the total compressed size for each brick
	out_patch->num_compress_blocks = comp_count; // Store the number of compressed blocks for each brick
  }
}


PIDX_return_code PIDX_wavelet_perform(PIDX_brick_res_precision_rst_id rst_id)
{
	PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index]; // first variable
	int patch_count = var0->brick_res_precision_io_restructured_super_patch_count; // local number of bricks

    int svi = rst_id->first_index;
    int evi = rst_id->last_index + 1;

    // Patch size (restructured size)
    uint64_t patch_x = rst_id->restructured_grid->patch_size[0];
    uint64_t patch_y = rst_id->restructured_grid->patch_size[1];
    uint64_t patch_z = rst_id->restructured_grid->patch_size[2];

    int min =  rst_id->restructured_grid->patch_size[0];
      for (int i = 1; i < 3; i++)
      {
    	if (rst_id->restructured_grid->patch_size[i] < min)
    	  min =  rst_id->restructured_grid->patch_size[i];
      }
      int max_wavelet_level = log2(min) - 2; // maximum wavelet level
      rst_id->restructured_grid->max_wavelet_level = max_wavelet_level; // Store the parameter into restructured_grid structure

	for (int g = 0; g < patch_count; ++g)
	{
		for (int v_start = svi; v_start < evi; v_start = v_start + 1)
		{
			// copy the size and offset to output
			PIDX_variable var_start = rst_id->idx->variable[v_start];
			PIDX_patch out_patch = var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch;

			// Calculate the bits per sample
			int bits = 0;
			PIDX_variable var = rst_id->idx->variable[v_start];
			bits = (var->bpv/8) * var->vps;

			int wavelet_level = max_wavelet_level;
			out_patch->wavelet_level = wavelet_level; // Store this parameter into PIDX_patch structure

			unsigned char* buf = var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer;
			uint64_t buffer_size = out_patch->size[0] * out_patch->size[1] * out_patch->size[2];
			uint64_t size = rst_id->restructured_grid->patch_size[0] * rst_id->restructured_grid->patch_size[1] * rst_id->restructured_grid->patch_size[2];
			unsigned char* res_buf = NULL;

			 // If the patch size is less than the brick size (e.g., 32x24x32), this patch should be padding with 0 to be 32x32x32.
			if (buffer_size < size)
			{
				res_buf = calloc(size * bits, sizeof(unsigned char));
				int index1 = 0; int index2 = 0;
				for (int i = 0; i < out_patch->size[2]; i++)
				{
					for (int j = 0; j < out_patch->size[1]; j++)
					{
						index1 = i * out_patch->size[1] * out_patch->size[0] + j * out_patch->size[0];
						index2 = i * patch_y * patch_x + j * patch_x;
						memcpy(&res_buf[index2 * bits], &buf[index1 * bits], out_patch->size[0] * bits);
					}
				}
				memcpy(buf, res_buf, size * bits);
			}
			free(res_buf);

			PIDX_wavelet_transform(buf, patch_x, patch_y, patch_z, bits, var->type_name, wavelet_level);
		}
	}
	return PIDX_success;
}



PIDX_return_code PIDX_zfp_compression_perform(PIDX_brick_res_precision_rst_id rst_id, unsigned long long* process_comp_size, int* max_patch_size)
{
	PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index]; // first variable
	int patch_count = var0->brick_res_precision_io_restructured_super_patch_count; // local number of bricks

    int svi = rst_id->first_index;
    int evi = rst_id->last_index + 1;

	uint64_t patch_x = rst_id->restructured_grid->patch_size[0];
	uint64_t patch_y = rst_id->restructured_grid->patch_size[1];
	uint64_t patch_z = rst_id->restructured_grid->patch_size[2];

	uint64_t patch_size[3] = {patch_x, patch_y, patch_z};

	// ZFP compression parameters
	rst_id->idx->comp_mode = 1;// The compression mode
	rst_id->idx->comp_param = 16; // The compression parameter

	rst_id->idx->procs_comp_buffer = (unsigned char*) malloc(patch_x * patch_y * patch_z * patch_count * 64);

	rst_id->compressed_sizes = (int *) malloc(patch_count * sizeof(int)); // local compressed bricks size array
	rst_id->patches_global_id = (int *) malloc(patch_count * sizeof(int));// local global id array
	rst_id->patches_rank = (int *) malloc(patch_count * sizeof(int)); // local brick rank array (which brick belongs to which process)

	for (int g = 0; g < patch_count; ++g)
	{
		int vars_comp_size = 0; // The size of a brick which contains all the variables

		for (int v_start = svi; v_start < evi; v_start = v_start + 1)
		{
			// copy the size and offset to output
			PIDX_variable var_start = rst_id->idx->variable[v_start];
			PIDX_patch out_patch = var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch;

			int bits = 0;
			PIDX_variable var = rst_id->idx->variable[v_start];
			bits = (var->bpv/8) * var->vps;

			unsigned char* buf = var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer;
			uint64_t size = rst_id->restructured_grid->patch_size[0] * rst_id->restructured_grid->patch_size[1] * rst_id->restructured_grid->patch_size[2];

		    uint64_t dc_dimension[3];
		    PIDX_calculate_level_dimension(dc_dimension, patch_size, rst_id->restructured_grid->max_wavelet_level);
		    uint64_t dc_size = dc_dimension[0] * dc_dimension[1] * dc_dimension[2];

		    // Wavelet compression
		    out_patch->compressed_buffer = (unsigned char*) malloc(size * bits * sizeof(unsigned char));
		    PIDX_wavelet_compression(out_patch->compressed_buffer, buf, bits, dc_size, patch_size, var->type_name, out_patch, rst_id->idx->comp_mode, rst_id->idx->comp_param);
		    out_patch->compressed_buffer = (unsigned char*) realloc(out_patch->compressed_buffer, out_patch->total_compress_size);
		    // copy all the brick buffers to one big buffer per process
		    memcpy(&rst_id->idx->procs_comp_buffer[(*process_comp_size)], out_patch->compressed_buffer, out_patch->total_compress_size);
		    (*process_comp_size) += out_patch->total_compress_size; // total compressed size per process
		    vars_comp_size += out_patch->total_compress_size; // total compressed size of all the variables

		}
		if (vars_comp_size > (*max_patch_size))
			(*max_patch_size) = vars_comp_size;

		rst_id->compressed_sizes[g] = vars_comp_size; // Store the compressed size per brick
		rst_id->patches_global_id[g] = var0->brick_res_precision_io_restructured_super_patch[g]->global_id; // store the global id per brick
		rst_id->patches_rank[g] = rst_id->idx_c->simulation_rank;
	}
	rst_id->idx->procs_comp_buffer = (unsigned char*)realloc(rst_id->idx->procs_comp_buffer, (*process_comp_size));
	return PIDX_success;
}


PIDX_return_code PIDX_brick_res_precision_rst_buf_aggregated_write(PIDX_brick_res_precision_rst_id rst_id)
{
	rst_id->idx->wave_start = MPI_Wtime();
	if (PIDX_wavelet_perform(rst_id) != PIDX_success)
	{
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
	}
	rst_id->idx->wave_end = MPI_Wtime();

	int max_patch_size = 0;
	int rank = rst_id->idx_c->simulation_rank; // The rank of processes
	unsigned long long process_comp_size = 0;   // The total size of a process after compression

	rst_id->idx->zfp_start = MPI_Wtime();
	if (PIDX_zfp_compression_perform(rst_id, &process_comp_size, &max_patch_size) != PIDX_success)
	{
        fprintf(stderr,"File %s Line %d\n", __FILE__, __LINE__);
        return PIDX_err_rst;
	}
	rst_id->idx->zfp_end = MPI_Wtime();

	rst_id->idx->agg_start = MPI_Wtime();
	/*********************** Aggregation **********************/
  PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index]; // first variable
  int patch_count = var0->brick_res_precision_io_restructured_super_patch_count; // local number of bricks

  unsigned long long max_file_size = rst_id->idx->max_file_size;   // Required file size
  int required_num_brick = rst_id->idx->required_num_brick; // Required number of bricks
  int process_count = rst_id->idx_c->simulation_nprocs; // Number of processes

  int total_patches_count = 0; // Total patches count
  MPI_Allreduce(&patch_count, &total_patches_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  rst_id->total_num_bricks = total_patches_count;

  // Total patches size over all the processes
  unsigned long long total_patches_size = 0;
  MPI_Allreduce(&process_comp_size, &total_patches_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  rst_id->idx->total_size = total_patches_size;

  int num_files = 0;
  if (required_num_brick > 0)
  {
	num_files = total_patches_count / required_num_brick; // The number of files
	num_files = (total_patches_count % required_num_brick == 0)? num_files: (num_files + 1);
	int mim_num_brick = total_patches_count / process_count; // The bottle line for the number of bricks
	if (num_files > process_count)
	{
	  printf("ERROR: The required number of bricks should be larger than %d\n", mim_num_brick);
	  return PIDX_err_io;
	}
  }
  if (max_file_size > 0)
  {
	// Maximun patch size over all the patches
	int max_pros_patch_size = 0;
	MPI_Allreduce(&max_patch_size, &max_pros_patch_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	// The required file size should larger than the maximum patch size
	if (max_file_size < (unsigned long long)max_pros_patch_size)
	{
	  printf("ERROR: The required file size should be larger than %d\n", max_pros_patch_size);
	  return PIDX_err_io;
	}
	// Total number of out files
	num_files = total_patches_size/max_file_size;
	unsigned long long remaining_size = total_patches_size % max_file_size;
	if (remaining_size != 0)
	{
	  unsigned long long tolerance = max_file_size * 0.5;
	  if (remaining_size > tolerance)
		num_files += 1;
	  else
		max_file_size += remaining_size/num_files;
	}
	// The total number of files should smaller than the number of processes
	unsigned long long min_file_size = total_patches_size/process_count;
	if (process_count < num_files)
	{
	  printf("ERROR: The required file size should be larger than %llu\n", min_file_size);
	  return PIDX_err_io;
	}
  }
  rst_id->idx->agg_counts = num_files; // The number of aggregation ranks
  // Patch counts array
  int patches_count_array[process_count];
  MPI_Allgather(&patch_count, 1, MPI_INT, patches_count_array, 1, MPI_INT, MPI_COMM_WORLD);
  // The displacement at which to place the incoming data from process i
  int displace_array[process_count];
  int value = 0;
  for (int i = 0; i < process_count; i++)
  {
	displace_array[i] = value;
	value += patches_count_array[i];
  }
  // Patches compressed size array
  int patch_size_array[total_patches_count];
  MPI_Allgatherv(rst_id->compressed_sizes, patch_count, MPI_INT,
		patch_size_array, patches_count_array, displace_array,
		MPI_INT, MPI_COMM_WORLD);
  // Patches global Id array
  int patch_global_id_array[total_patches_count];
  MPI_Allgatherv(rst_id->patches_global_id, patch_count, MPI_INT,
		patch_global_id_array, patches_count_array, displace_array,
		MPI_INT, MPI_COMM_WORLD);
  // Patches rank array
  int patch_rank_array[total_patches_count];
  MPI_Allgatherv(rst_id->patches_rank, patch_count, MPI_INT,
		patch_rank_array, patches_count_array, displace_array,
		MPI_INT, MPI_COMM_WORLD);

  rst_id->idx->is_aggregator = 0;

  // Decide aggregation ranks
  int aggregate_rank[num_files];
  int div = process_count/num_files;
  int agg_index = 0;
  for (int i = 0; i < process_count; i++)
  {
	if (i % div == 0 && agg_index < num_files)
	{
	  aggregate_rank[agg_index] = i;
	  agg_index++;
	  if (rank == i)
		  rst_id->idx->is_aggregator = 1;
	}
  }

  // Record which brick should send to which aggregate
  int aggregate_record[total_patches_count];
  memset(aggregate_record, -1, total_patches_count * sizeof(int));
  unsigned long long agg_size[num_files]; // The file size per aggregate

  if (required_num_brick > 0)
  {
	int i = 0;
	for (i = 0 ; i < (num_files - 1); i++)
	{
	  agg_size[i] = 0;
	  for (int j = 0; j < required_num_brick; j++)
	  {
		int id = i*required_num_brick + j;
		aggregate_record[id] = aggregate_rank[i];
		agg_size[i] += patch_size_array[id];
	  }
	}
	int index = i*required_num_brick;
	agg_size[i] = 0;
	while (index < total_patches_count)
	{
	  aggregate_record[index] = aggregate_rank[i];
	  agg_size[i] += patch_size_array[index];
	  index++;
	}
  }

  if (max_file_size > 0)
  {
	// Traverse all the bricks and assign them to aggregates based on the required file size
	for (int i = 0 ; i < num_files; i++)
	{
	  agg_size[i] = 0;
	  for (int j = 0; j < total_patches_count; j++)
	  {
		if (aggregate_record[j] == -1 && (agg_size[i] + patch_size_array[j]) < max_file_size)
		{
		  aggregate_record[j] = aggregate_rank[i];
		  agg_size[i] += patch_size_array[j];
		}
	  }
	}
	// Assign remaining bricks to the aggregate which holds the minimum size
	for (int i = (total_patches_count - 1); i > -1; i--)
	{
	  if (aggregate_record[i] == -1)
	  {
		unsigned long long min_size = max_file_size;
		int min_rank = -1;
		int min_index = -1;
		for (int j = 0; j < num_files; j++)
		{
		  if (agg_size[j] < min_size)
		  {
			min_size = agg_size[j];
			min_rank = aggregate_rank[j];
			min_index = j;
		  }
		}
		aggregate_record[i] = min_rank;
		agg_size[min_index] += patch_size_array[i];
	  }
	}
  }

  rst_id->idx->agg_size = 0; // The size per aggregate
  unsigned char* aggregate_buffer = NULL;
  for (int i = 0; i < num_files; i++)
  {
	if (rank == aggregate_rank[i])
	{
	  aggregate_buffer = (unsigned char*) malloc(agg_size[i]); // aggregate buffer
	  rst_id->idx->agg_size = agg_size[i];
	}
  }

  unsigned long long local_brick_disp[patch_count]; // The displacement for each brick per process
  int start_index = displace_array[rank]; // The index of first brick for each process
  int disp = 0;
  for (int i  = 0; i < patch_count; i++)
  {
	int id = start_index + i;
	local_brick_disp[i] = disp;
	disp += patch_size_array[id];
  }

  int local_recv_array[total_patches_count];
  int local_id_array[total_patches_count];
  int recv_ind = 0;
  for (int i = 0; i < total_patches_count; i++)
  {
	  if (aggregate_record[i] == rank && patch_rank_array[i] != rank)
	  {
		  local_recv_array[recv_ind] = patch_rank_array[i];
		  local_id_array[recv_ind] = i;
		  recv_ind++;
	  }
  }

  unsigned long long agg_cur_size = 0; // Current aggregate size
  int index = 0;
  rst_id->idx->agg_patch_array = (int*) malloc(total_patches_count * sizeof(int));
  rst_id->idx->agg_patches_size_array = (int*) malloc(total_patches_count * sizeof(int));
  int owned_patch_count = 0;
  MPI_Request req[patch_count + recv_ind];
  MPI_Status stat[patch_count + recv_ind];
  for (int i = 0; i < patch_count; i++)
  {
	  int id = start_index + i;
	  if (rank == aggregate_record[id])
	  {
		  memcpy(&aggregate_buffer[agg_cur_size], &rst_id->idx->procs_comp_buffer[local_brick_disp[i]], patch_size_array[id]);
		  agg_cur_size += patch_size_array[id];
		  rst_id->idx->agg_patch_array[owned_patch_count] = patch_global_id_array[id];
		  owned_patch_count++;
	  }
	  else
	  {
		  MPI_Isend(&rst_id->idx->procs_comp_buffer[local_brick_disp[i]], patch_size_array[id], MPI_UNSIGNED_CHAR, aggregate_record[id], 0,
		  			MPI_COMM_WORLD, &req[index]);
		  MPI_Wait(&req[index], &stat[index]);
		  index++;
	  }
  }

  for (int i = 0; i < recv_ind; i++)
  {
	  MPI_Irecv(&aggregate_buffer[agg_cur_size], patch_size_array[local_id_array[i]], MPI_UNSIGNED_CHAR, local_recv_array[i],
	  			0, MPI_COMM_WORLD, &req[index]);
	  MPI_Wait(&req[index], &stat[index]);
	  agg_cur_size += patch_size_array[local_id_array[i]];
	  rst_id->idx->agg_patch_array[owned_patch_count] = patch_global_id_array[local_id_array[i]];
	  rst_id->idx->agg_patches_size_array[owned_patch_count] = patch_size_array[local_id_array[i]];
	  owned_patch_count++;
	  index++;
  }
  rst_id->idx->agg_patch_array = (int*) realloc(rst_id->idx->agg_patch_array, owned_patch_count * sizeof(int));
  rst_id->idx->agg_owned_patch_count = owned_patch_count;

  rst_id->idx->agg_end = MPI_Wtime();

  /*********************** Write data out **********************/
  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  // write out files
  if (agg_cur_size > 0)
  {
	char *file_name;
	file_name = malloc(PATH_MAX * sizeof(*file_name));
	memset(file_name, 0, PATH_MAX * sizeof(*file_name));

	sprintf(file_name, "%s/time%09d/%d", directory_path, rst_id->idx->current_time_step, rank);
	int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

	uint64_t buffer_size =  agg_cur_size;
	uint64_t write_count = pwrite(fp, aggregate_buffer, buffer_size, 0);
	if (write_count != buffer_size)
	{
	  fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
	  return PIDX_err_io;
	}
	close(fp);
	free(file_name);
  }

  free(aggregate_buffer);
  free(directory_path);
  /***********************************************************/

  return PIDX_success;
}


PIDX_return_code PIDX_brick_res_precision_rst_buf_aggregated_read(PIDX_brick_res_precision_rst_id rst_id)
{
  int g = 0;
  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  char time_template[512];
  sprintf(time_template, "%%s/%s/%%d_%%d", rst_id->idx->filename_time_template);

  PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index];
  for (g = 0; g < var0->brick_res_precision_io_restructured_super_patch_count; ++g)
  {
    // loop through all groups
    char *file_name;
    file_name = malloc(PATH_MAX * sizeof(*file_name));
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));

    sprintf(file_name, time_template, directory_path, rst_id->idx->current_time_step, rst_id->idx_c->simulation_rank, g);
    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    int v_start = 0;
    int svi = rst_id->first_index;
    int evi = rst_id->last_index + 1;
    for (v_start = svi; v_start < evi; v_start = v_start + 1)
    {
      // copy the size and offset to output
      PIDX_variable var_start = rst_id->idx->variable[v_start];
      PIDX_patch out_patch = var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch;

      int bits = 0;
      PIDX_variable var = rst_id->idx->variable[v_start];
      bits = (var->bpv/8) * var->vps;

      uint64_t data_offset = 0, v1 = 0;
      for (v1 = 0; v1 < v_start; v1++)
        data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (rst_id->idx->variable[v1]->vps * (rst_id->idx->variable[v1]->bpv/8)));

      uint64_t buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
      uint64_t write_count = pread(fp, var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer, buffer_size, data_offset);
      if (write_count != buffer_size)
      {
        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
        return PIDX_err_io;
      }
    }
    close(fp);
    free(file_name);
  }

  free(directory_path);

  return PIDX_success;
}
