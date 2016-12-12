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

#include "PIDX_file_handler.h"


PIDX_return_code PIDX_debug_disable_restructuring(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_rst = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_chunking(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_chunk = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_compression(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_compress = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_hz(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_hz = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_agg(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_agg = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_disable_io(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_do_io = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_rst(PIDX_file file, int debug_rst)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_rst = debug_rst;

  return PIDX_success;
}



PIDX_return_code PIDX_debug_hz(PIDX_file file, int debug_hz)
{
  if(!file)
    return PIDX_err_file;

  file->idx_dbg->debug_hz = debug_hz;

  return PIDX_success;
}


PIDX_return_code PIDX_disable_rst(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx->enable_rst = 0;

  return PIDX_success;
}



PIDX_return_code PIDX_disable_agg(PIDX_file file)
{
  if(!file)
    return PIDX_err_file;

  file->idx->enable_agg = 0;

  return PIDX_success;
}


PIDX_return_code PIDX_dump_rst_info(PIDX_file file, int dump_rst_info)
{
  if(!file)
    return PIDX_err_file;

  if (dump_rst_info == PIDX_RST_DUMP_INFO)
  {
    file->idx_dbg->dump_rst_info = 1;
  }
  else if (dump_rst_info == PIDX_SIMULATE_RST_AND_DUMP_INFO)
  {
    file->idx_dbg->simulate_rst = 1;
    file->idx_dbg->dump_rst_info = 1;
  }
  else if (dump_rst_info == PIDX_NO_IO_AND_SIMULATE_RST_AND_DUMP_INFO)
  {
    file->idx_dbg->simulate_rst = 1;
    file->idx_dbg->simulate_rst_io = 1;
    file->idx_dbg->dump_rst_info = 1;
  }

  char filename_skeleton[512];
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_dbg->rst_dump_dir_name, "%s_%d_rst_dump", filename_skeleton, file->idx->current_time_step);

  return PIDX_success;
}


PIDX_return_code PIDX_dump_agg_info(PIDX_file file, int dump_agg_info)
{
  if(!file)
    return PIDX_err_file;

  char filename_skeleton[512];
  file->idx_dbg->dump_agg_info = dump_agg_info;
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_dbg->agg_dump_dir_name, "%s_agg_dump", filename_skeleton);

  return PIDX_success;
}



PIDX_return_code PIDX_dump_io_info(PIDX_file file, int dump_io_info)
{
  if(!file)
    return PIDX_err_file;

  char filename_skeleton[512];
  file->idx_dbg->dump_io_info = dump_io_info;
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_dbg->io_dump_dir_name, "%s_io_dump", filename_skeleton);

  return PIDX_success;
}



PIDX_return_code PIDX_dump_process_state(PIDX_file file, int process_state)
{
  if(!file)
    return PIDX_err_file;

  char filename_skeleton[512];
  file->idx_dbg->dump_process_state = process_state;
  strncpy(filename_skeleton, file->idx->filename, strlen(file->idx->filename) - 4);
  filename_skeleton[strlen(file->idx->filename) - 4] = '\0';
  sprintf(file->idx_dbg->process_state_dump_dir_name, "%s_%d_process_state_dump", filename_skeleton, file->idx->current_time_step);

  return PIDX_success;
}
