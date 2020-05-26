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
void wavelet_helper(unsigned char* buf, int step, int ng_step, int flag, int bits, uint64_t x, uint64_t y, uint64_t z, char* type_name)
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
void PIDX_wavelet_transform(unsigned char* buffer, uint64_t x, uint64_t y, uint64_t z, int bits, char* type_name)
{
	// Calculate the max wavelet level based on the min dimensional value
	uint64_t size[3] = {x, y, z};
	int min = size[0];
	for (int i = 1; i < 3; i++)
	{
		if (size[i] < min)
			min = size[i];
	}

	int max_wavelet_level = log2(min);

	// Get random wavelet level ( >= 1)
	time_t t;
	srand((unsigned) time(&t));
	int wavelet_level = rand()%max_wavelet_level;
	wavelet_level = (wavelet_level < 1) ? 1: wavelet_level;

	for (int level = 1; level <= wavelet_level; level++)
	{
		int step = pow(2, level);
		int ng_step = step/2;

        // Calculate x-dir
        wavelet_helper(buffer, step, ng_step, 0, bits, x, y, z, type_name);
        // Calculate y-dir
        wavelet_helper(buffer, step, ng_step, 1, bits, x, y, z, type_name);
        // Calculate z-dir
        wavelet_helper(buffer, step, ng_step, 2, bits, x, y, z, type_name);
	}
}


// ZFP compression
unsigned char* PIDX_compress_3D_float(unsigned char* buf, int dim_x, int dim_y, int dim_z, float param, int flag, char* type_name)
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
    unsigned char* output = (unsigned char*) malloc(max_compressed_bytes);
    bitstream* stream = stream_open(&output[0], max_compressed_bytes);
    zfp_stream_set_bit_stream(zfp, stream);
    size_t compressed_bytes = zfp_compress(zfp, field);
    if (compressed_bytes == 0)
        puts("ERROR: Something wrong happened during compression\n");
    output = (unsigned char*) realloc(output, compressed_bytes);
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    return output;
}


PIDX_return_code PIDX_brick_res_precision_rst_buf_aggregated_write(PIDX_brick_res_precision_rst_id rst_id)
{
  int g = 0;
  char *directory_path;
  directory_path = malloc(sizeof(*directory_path) * PATH_MAX);
  memset(directory_path, 0, sizeof(*directory_path) * PATH_MAX);
  strncpy(directory_path, rst_id->idx->filename, strlen(rst_id->idx->filename) - 4);

  PIDX_variable var0 = rst_id->idx->variable[rst_id->first_index];

//  printf("rank %d: number_of_bricks: %d\n", rst_id->idx_c->simulation_rank, var0->brick_res_precision_io_restructured_super_patch_count);

  for (g = 0; g < var0->brick_res_precision_io_restructured_super_patch_count; ++g)
  {
    // loop through all groups
    char *file_name;
    file_name = malloc(PATH_MAX * sizeof(*file_name));
    memset(file_name, 0, PATH_MAX * sizeof(*file_name));

    sprintf(file_name, "%s/time%09d/%d_%d_%d", directory_path, rst_id->idx->current_time_step, rst_id->idx_c->simulation_rank, g, var0->brick_res_precision_io_restructured_super_patch[g]->global_id);
    int fp = open(file_name, O_CREAT | O_WRONLY, 0664);

    int v_start = 0;
    int svi = rst_id->first_index;
    int evi = rst_id->last_index + 1;
    for (v_start = svi; v_start < evi; v_start = v_start + 1)
    {
      // copy the size and offset to output
      PIDX_variable var_start = rst_id->idx->variable[v_start];
      PIDX_patch out_patch = var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch;

//      printf("rank %d: %dx%dx%d\n", rst_id->idx_c->simulation_rank, out_patch->size[0], out_patch->size[1], out_patch->size[2]);


      int bits = 0;
      PIDX_variable var = rst_id->idx->variable[v_start];
      bits = (var->bpv/8) * var->vps;

      // If the patch size is less than the brick size (e.g., 32x24x32), this patch should be padding with 0 to be 32x32x32.
      unsigned char* buf = var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer;
      uint64_t patch_x = rst_id->restructured_grid->patch_size[0];
      uint64_t patch_y = rst_id->restructured_grid->patch_size[1];
      uint64_t patch_z = rst_id->restructured_grid->patch_size[2];

      uint64_t buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2];
      uint64_t size = patch_x * patch_y * patch_z;

      if (buffer_size < size)
      {
//    	  printf("%dx%dx%d\n", out_patch->size[0], out_patch->size[1], out_patch->size[2]);
    	  unsigned char* res_buf = calloc(size * bits, sizeof(unsigned char));
    	  int index1 = 0; int index2 = 0;
    	  for (int i = 0; i < out_patch->size[2]; i++)
    	  {
    		  for (int j = 0; j < out_patch->size[1]; j++)
    		  {
    			  memcpy(&res_buf[index2], &buf[index1], out_patch->size[0] * bits);
    			  index1 = i * out_patch->size[1] * out_patch->size[0] + j * out_patch->size[0];
    			  index2 = i * patch_y * patch_x + j * patch_x;
    		  }
    	  }
    	  // pass this restructured buffer to wavelet transform
    	  PIDX_wavelet_transform(res_buf, patch_x, patch_y, patch_z, bits, var->type_name);
    	  free(res_buf);
      }
      else
      {
    	  // pass original buffer to wavelet transform
    	  PIDX_wavelet_transform(buf, patch_x, patch_y, patch_z, bits, var->type_name);
      }

//      float a = 0;
//      if (var0->brick_res_precision_io_restructured_super_patch[g]->global_id == 0)
//      {
//    	  for (int i = 0; i < size; i++)
//    	  {
//    		  memcpy(&a, &buf[i*bits], bits);
//    		  printf("%f\n", a);
//    	  }
//      }

      uint64_t data_offset = 0, v1 = 0;
      for (v1 = 0; v1 < v_start; v1++)
        data_offset = data_offset + (out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * (rst_id->idx->variable[v1]->vps * (rst_id->idx->variable[v1]->bpv/8)));

//      uint64_t buffer_size =  out_patch->size[0] * out_patch->size[1] * out_patch->size[2] * bits;
//      uint64_t write_count = pwrite(fp, var_start->brick_res_precision_io_restructured_super_patch[g]->restructured_patch->buffer, buffer_size, data_offset);
//      if (write_count != buffer_size)
//      {
//        fprintf(stderr, "[%s] [%d] pwrite() failed.\n", __FILE__, __LINE__);
//        return PIDX_err_io;
//      }
    }
    close(fp);
    free(file_name);
  }

  // You can access max file size in
  // rst_id->idx->max_file_size

  free(directory_path);

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
