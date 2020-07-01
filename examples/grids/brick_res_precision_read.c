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

/*
  PIDX write example

  In this example we show how to write data using the PIDX library.

  We consider a global 3D regular grid domain that we will call
  global domain (g).
  This global domain represents the grid space where all the data are stored.

  In a parallel environment each core (e.g. MPI rank) owns a portion of the data
  that has to be written on the disk. We refer to this portion of the domain as
  local domain (l).

  In this example we well see how to execute parallel write with PIDX of a
  syntethic dataset.

  In the following picture is represented a sample domain decomposition
  of the global domain (l) in per-core local domains (l), sometimes referred
  as patches.
  In this example all the local domains have same dimesions for simplicity.
  PIDX supports different number and sizes of patches per core.

             *---------*--------*
           /         /         /| P7
          *---------*---------* |
         /         /         /| |
        *---------*---------* | *
        |         |         | |/|           --------->        IDX Data format
        |         |         | * |
        | P4      | P5      |/| | P3
        *---------*---------* | *
        |         |         | |/
        |         |         | *
        | P0      | P1      |/
        *---------*---------*

*/
#if !defined _MSC_VER
#include <unistd.h>
#endif
#include <stdarg.h>
#include <stdint.h>
#include <ctype.h>
#include <PIDX.h>

#if defined _MSC_VER
  #include "utils/PIDX_windows_utils.h"
#endif

#include "pidx_examples_utils.h"

char input_file_template[512];
char input_file[512];
int variable_index = 0;
int current_ts = 0;
int wavelet_level = 0;
int block_id = -1;
static int rst_box_size[PIDX_MAX_DIMENSIONS];
int bpv[MAX_VAR_COUNT];
char type_name[MAX_VAR_COUNT][512];
int vps[MAX_VAR_COUNT];
char var_name[MAX_VAR_COUNT][512];

static PIDX_point rst_box;
PIDX_variable* variable;


static char *usage = "Serial Usage: ./idx_read -t 0 -v 0 -f input_idx_file_name\n"
                     "Parallel Usage: mpirun -n 8 ./idx_read -t 0 -v 0 -f input_idx_file_name\n"
                     "  -f: IDX input filename\n"
                     "  -t: time step index to read\n"
                     "  -v: variable index to read\n"
		             "  -l: wavelet level to read\n"
		             "  -b: blocks id to read\n";

static void parse_args(int argc, char **argv);
static int parse_metadata();
static void create_pidx_point_and_access_new();
static void set_pidx_file(int ts);
static void set_pidx_variable(int var);


int main(int argc, char **argv)
{
	// Init MPI and MPI vars (e.g. rank and process_count)
	init_mpi(argc, argv);

	// corresponing variables
	parse_args(argc, argv);

	parse_metadata();

	printf("%d x %d x %d\n", global_box_size[0], global_box_size[1], global_box_size[2]);
	printf("%d x %d x %d\n", rst_box_size[0], rst_box_size[1], rst_box_size[2]);
	printf("%d %s %d %d\n", variable_count, type_name[0], bpv[0], vps[0]);

	create_pidx_point_and_access_new();
	variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
    memset(variable, 0, sizeof(*variable) * variable_count);

    set_pidx_file(current_ts);
    for (int var = 0; var < variable_count; var++)
      set_pidx_variable(var);

    // PIDX_close triggers the actual read on the disk
    // of the variables that we just set
//    PIDX_close(file);

    // Close access and free memory
    if (PIDX_close_access(p_access) != PIDX_success)
      terminate_with_error_msg("PIDX_close_access");

    free(variable);
    variable = 0;

	shutdown_mpi();

	return 0;
}

static void parse_args(int argc, char **argv)
{
  char flags[] = "f:t:v:l:b:";
  int one_opt = 0;

  while ((one_opt = getopt(argc, argv, flags)) != EOF)
  {
    /* postpone error checking for after while loop */
    switch (one_opt)
    {
		case('f'): // input file name
		  if (sprintf(input_file_template, "%s", optarg) < 0)
			terminate_with_error_msg("Invalid output file name template\n%s", usage);
		  sprintf(input_file, "%s%s", input_file_template, ".idx");
		  break;

		case('t'): // timesteps to read
		  if (sscanf(optarg, "%d", &current_ts) < 0)
			terminate_with_error_msg("Invalid variable file\n%s", usage);
		  break;

		case('v'): // variable to read
		  if (sscanf(optarg, "%d", &variable_index) < 0)
			terminate_with_error_msg("Invalid variable file\n%s", usage);
		  break;

		case('l'): // wavelet level to read
		  if (sscanf(optarg, "%d", &wavelet_level) < 0)
			terminate_with_error_msg("Invalid variable file\n%s", usage);
		  break;

		case('b'): // block id to read
		  if (sscanf(optarg, "%d", &block_id) < 0)
			terminate_with_error_msg("Invalid brick size\n%s", usage);
		  break;

		default:
		  terminate_with_error_msg("Wrong arguments\n%s", usage);
    }
  }
}

static int parse_metadata()
{
  FILE *fp = fopen(input_file, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Error Opening %s\n", input_file);
    return PIDX_err_file;
  }
  int i = 0;
  int variable_counter = 0, count = 0, len = 0;
  char *pch1;
  char line [ 512 ];
  while (fgets(line, sizeof (line), fp) != NULL)
  {
	line[strcspn(line, "\r\n")] = 0;
	if (strcmp(line, "(box)") == 0)
	{
	  if ( fgets(line, sizeof line, fp) == NULL)
	    return PIDX_err_file;
	  line[strcspn(line, "\r\n")] = 0;
	  pch1 = strtok(line, " ");
	  while (pch1 != NULL)
	  {
		if (strcmp(pch1, "0") != 0)
		{
		  global_box_size[i] = atoi(pch1) + 1;
		  i++;
		}
		pch1 = strtok(NULL, " ");
	  }
	}
	if (strcmp(line, "(restructure box size)") == 0)
	{
	  if ( fgets(line, sizeof line, fp) == NULL)
	    return PIDX_err_file;
	  line[strcspn(line, "\r\n")] = 0;
	  pch1 = strtok(line, " ");
	  i = 0;
	  while (pch1 != NULL)
	  {
		rst_box_size[i] = atoi(pch1);
		i++;
		pch1 = strtok(NULL, " ");
	  }
	}
	if (strcmp(line, "(fields)") == 0)
	{
	  if ( fgets(line, sizeof line, fp) == NULL)
		return PIDX_err_file;
	  line[strcspn(line, "\r\n")] = 0;
	  count = 0;
	  variable_counter = 0;
	  while (line[X] != '(')
	  {
		pch1 = strtok(line, " +");
		while (pch1 != NULL)
		{
		  if (count == 0)
		  {
			char* temp_name = strdup(pch1);
			strcpy(var_name[variable_counter], temp_name);
			free(temp_name);
		  }

		  if (count == 1)
		  {
			len = strlen(pch1) - 1;
			if (pch1[len] == '\n')
			  pch1[len] = 0;

			strcpy(type_name[variable_counter], pch1);
			int ret;
			int bits_per_sample = 0;
			int sample_count = 0;
			ret = PIDX_values_per_datatype(type_name[variable_counter], &sample_count, &bits_per_sample);
			if (ret != PIDX_success)  return PIDX_err_file;

			bpv[variable_counter] = bits_per_sample;
			vps[variable_counter] = sample_count;
		  }
		  count++;
		  pch1 = strtok(NULL, " +");
		}
		count = 0;
		if ( fgets(line, sizeof line, fp) == NULL)
		  return PIDX_err_file;
		line[strcspn(line, "\r\n")] = 0;
		variable_counter++;
	  }
	  variable_count = variable_counter;
	}
  }
  fclose(fp);
  return PIDX_success;
}

static void create_pidx_point_and_access_new()
{
  PIDX_set_point(global_size, global_box_size[X], global_box_size[Y], global_box_size[Z]);
  PIDX_set_point(rst_box, rst_box_size[X], rst_box_size[Y], rst_box_size[Z]);

  //  Creating access
  PIDX_create_access(&p_access);
  PIDX_set_mpi_access(p_access, MPI_COMM_WORLD);

  PIDX_create_metadata_cache(&cache);

  return;
}

static void set_pidx_file(int ts)
{
  PIDX_return_code ret;

  // Create IDX file
  ret = PIDX_file_create(input_file, PIDX_MODE_CREATE, p_access, global_size, &file);
  if (ret != PIDX_success)
    terminate_with_error_msg("PIDX_file_create\n");

  // Set the current timestep
  PIDX_set_current_time_step(file, ts);
  // Set the number of variables
  PIDX_set_variable_count(file, variable_count);

  // Select I/O mode (PIDX_IDX_IO for the multires, PIDX_RAW_IO for non-multires)
  PIDX_set_io_mode(file, PIDX_BRICK_RES_PRECISION_IO);

  PIDX_set_restructuring_box(file, rst_box);

  // ADD BY KE
  PIDX_set_required_wavelet_level(file, wavelet_level);

  PIDX_set_required_block_id(file, block_id);

  return;
}


static void set_pidx_variable(int var)
{
  PIDX_return_code ret = 0;

  // Set variable name, number of bits, typename
  ret = PIDX_variable_create(var_name[var], bpv[var] * vps[var], type_name[var], &variable[var]);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_create");

//  // Set the variable offset and size of the local domain,
//  // where the data is in memory (data) and what is its layout in memory (row major)
//  ret = PIDX_variable_write_data_layout(variable[var], local_offset, local_size, data[var], PIDX_row_major);
//  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_variable_write_data_layout");

//  // Tell PIDX that we want to write this variable
//  ret = PIDX_append_and_write_variable(file, variable[var]);
//  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_append_and_write_variable");

  return;
}
