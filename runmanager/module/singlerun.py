#!/usr/bin/env python3

import configparser
import copy
import logging
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree

import bjob
import bsub

logger = logging.getLogger('__main__').getChild(__name__)
SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAX_NPROC = 18 # not 20

#______________________________________________________________________________
class SingleRun(object):
  ''' Manager class of single run '''

  #____________________________________________________________________________
  def __init__(self, tag, runinfo):
    self.__status = 'INIT'
    # true: success, false: failure, 0: running, 1: killed
    self.__bjob_status  = None
    self.__merge_status = 'INIT'
    self.__tag = tag
    self.__key = runinfo['key']
    self.__bin_path = runinfo['bin']
    self.__is_dstgenfit = ('Dst' in self.__bin_path or 'Genfit' in self.__bin_path)
    self.__is_dst = ( self.__is_dstgenfit and
                      'dstin' in runinfo and
                      len(runinfo['dstin']) > 0)
    if self.__is_dst:
      self.__dstin_path = runinfo['dstin']
    else:
      self.__dstin_path = None
    self.__conf_path = runinfo['conf']
    self.__root_path = runinfo['root']
    self.__macro_path = None
    self.__log_path = None
    self.__merge_log_path = None
    self.__nproc = runinfo['nproc']
    if self.__nproc > MAX_NPROC:
      logger.warning(f'{self.__key} nproc must be {MAX_NPROC} or less'+
                     f'=> change {MAX_NPROC}')
      self.__nproc = MAX_NPROC
    self.__buff_path = runinfo['buff']
    self.__nevents = runinfo['nevents']
    self.__div_unit = runinfo['unit']
    self.__queue = runinfo['queue']
    self.__qmerge = runinfo['qmerge']
    self.__option = ''
    self.__start_time = time.time()
    self.__elapsed_time = 0
    self.__basename = os.path.splitext(os.path.basename(self.__root_path))[0]
    tmp_dir = os.path.join(SCRIPT_DIR, 'tmp')
    if not os.path.exists(tmp_dir):
      os.mkdir(tmp_dir)
    self.__dummy_dir = tempfile.TemporaryDirectory(dir=tmp_dir)
    self.__elem_list = list()
    self.__root_list = list()
    self.__log_list = list()
    self.__bjob_list = list()
    self.__jobid_list = list()
    self.__jobstat_list = list()
    self.__merging_process = None
    self.__merging_job = None
    self.__make_log()
    # self.__make_prefetch()
    self.__dump_init_info()
    self.__make_element()
    self.__ncomplete = 0

  #____________________________________________________________________________
  def clear(self):
    ''' Clear divided root files and keep logs. '''
    if not self.__qmerge is None:
      self.clear_divided_root_file()

  #____________________________________________________________________________
  def clear_all(self):
    ''' Clear divided root files and logs. '''
    self.clear()
    self.clear_all_log()

  #____________________________________________________________________________
  def clear_all_log(self):
    ''' Clear all logs. '''
    self.clear_log_list()
    self.clear_main_log()
    self.clear_merging_log()

  #____________________________________________________________________________
  def clear_divided_root_file(self):
    ''' Clear divided root files. '''
    for item in self.__root_list:
      if os.path.exists(item):
        os.remove(item)

  #____________________________________________________________________________
  def clear_log_list(self):
    ''' Clear logs for divided runs. '''
    for item in self.__log_list:
      if os.path.exists(item):
        os.remove(item)

  #____________________________________________________________________________
  def clear_main_log(self):
    ''' Clear main log. '''
    if self.__log_path is not None:
      if os.path.exists(self.__log_path):
        os.remove(self.__log_path)

  #____________________________________________________________________________
  def clear_merging_log(self):
    ''' Clear merging log. '''
    if self.__merge_log_path is not None:
      if os.path.exists(self.__merge_log_path):
        os.remove(self.__merge_log_path)

  #____________________________________________________________________________
  def execute(self):
    ''' Execute batch job. '''
    if self.__bjob_status is not None:
      return
    self.__make_macro()
    self.__make_root_list()
    self.__make_log_list()
    if self.__bjob_status is False:
      return
    for i, elem in enumerate(self.__elem_list):
      if os.path.exists(self.__root_list[i]):
        os.remove(self.__root_list[i])
      job = bsub.BSub(self, elem,
                      self.__conf_path,
                      self.__root_list[i],
                      self.__log_list[i])
      if self.__is_dst is True:
        job.set_dstin_list(self.__dstin_path)
      job.execute()
      self.__bjob_list.append(job)
      self.__jobid_list.append(None)
      self.__jobstat_list.append(job.get_status())
    for i, item in enumerate(self.__bjob_list):
      self.__dump_log(f'pid[bjob({i})]', item.get_process_id())
    self.__dump_log(None, '_'*80)
    self.__bjob_status = 0
    # self.__status = 0

  #____________________________________________________________________________
  def finalize(self, keep_log=False):
    ''' Finalize. '''
    self.update_status()
    self.__dump_log('end', time.ctime(time.time()))
    if self.__status == 'DONE':
      self.clear_all()
    elif self.__status == 'TERMINATED':
      if keep_log:
        self.clear()
      else:
        self.clear_all()
    else:
      self.clear()

  #____________________________________________________________________________
  def get_bin_path(self):
    ''' Get bin path. '''
    return self.__bin_path

  #____________________________________________________________________________
  def get_elapsed_time(self):
    ''' Get elapsed time. '''
    # self.update_status()
    if (self.__status != 'DONE' and
        self.__status != 'FAILED' and
        self.__status != 'TERMINATED'):
      self.__elapsed_time = time.time() - self.__start_time
    return self.__elapsed_time

  #____________________________________________________________________________
  def get_info(self):
    ''' Get run info. '''
    # self.update_status()
    info = dict()
    info['queue'] = self.__queue
    info['qmerge'] = self.__qmerge
    info['nproc'] = self.__nproc
    info['unit'] = self.__div_unit
    info['nev'] = self.__nevents
    info['bin'] = self.__bin_path
    info['conf'] = self.__conf_path
    info['root'] = self.__root_path
    info['time'] = self.get_elapsed_time()
    info['stat'] = self.get_status()
    info['nseg'] = self.get_nsegs()
    info['prog'] = self.get_progress()
    return {format(self.__key): info}

  #____________________________________________________________________________
  def get_key(self):
    ''' Get key. '''
    return self.__key

  #____________________________________________________________________________
  def get_macro_path(self):
    ''' Get macro path. '''
    return self.__macro_path

  #____________________________________________________________________________
  def get_nsegs(self):
    ''' Get number of elements. '''
    return len(self.__elem_list)

  #____________________________________________________________________________
  def get_option(self):
    ''' Get option. '''
    return self.__option

  #____________________________________________________________________________
  def get_progress(self):
    ''' Get progress. '''
    return self.__ncomplete

  #____________________________________________________________________________
  # def get_prefetch_path(self):
  #   ''' Get prefetch path. '''
  #   return self.__prefetch_path

  #____________________________________________________________________________
  def get_queue_merge(self):
    ''' Get queue for merging job. '''
    return self.__qmerge

  #____________________________________________________________________________
  def get_queue(self):
    ''' Get queue. '''
    return self.__queue

  #____________________________________________________________________________
  def get_start_time(self):
    ''' Get starting time. '''
    return time.ctime(self.__start_time)

  #____________________________________________________________________________
  def get_status(self):
    ''' Get status. '''
    # if (self.__status == 'TERMINATED' or
    #     self.__status == 'DONE' or
    #     self.__status == 'FAILED'):
    #   return self.__status
    # self.update_status()
    return self.__status

  #____________________________________________________________________________
  def kill_all(self):
    ''' Kill all processes. '''
    self.update_status()
    if self.__status == 'DONE':
      return
    self.kill_bjob()
    self.kill_merge()
    self.update_status()

  #____________________________________________________________________________
  def kill_bjob(self):
    ''' Kill bjob. '''
    if self.__bjob_status is True:
      return
    for job in self.__bjob_list:
      stat = job.get_status()
      if stat == 0 or stat == 1:
        job.kill()
        self.__bjob_status = 1
        jpid = job.get_process_id() if stat == 0 else job.get_job_id()
        buff = f'Process was killed [jid/pid: {jpid}]'
        self.__dump_log('kill_bjob', buff)
    if self.__bjob_status == 1:
      self.__dump_log(None, '_'*80)

  #____________________________________________________________________________
  def kill_merge(self):
    ''' Kill merging job. '''
    self.__update_merging_status()
    if self.__merge_status == 'DONE' or self.__merge_status == 'FAILED':
      return
    if self.__merge_status == 'RUNNING':
      logger.error('Killing merging job')
      self.__merging_job.kill()
      self.__merge_status = 'TERMINATED'
      buff = f'merging job was killed [jid: {self.__merging_job.get_job_id()}]'
      self.__dump_log('kill_merge', buff)
      self.__dump_log(None, '_'*80)

  #____________________________________________________________________________
  def merge_root_files(self):
    ''' Merge root files. '''
    if self.__merge_status != 'INIT':
      return
    if self.__qmerge is None:
      self.__merge_status = 'DONE'
      return
    if len(self.__elem_list) == 1:
      os.rename(self.__root_list[0], self.__root_path)
      self.__merge_status = 'DONE'
      return
    # size = 0
    # for item in self.__root_list:
    #   size += os.path.getsize(item)
    # if size > 3.8*(10**9):
    #   qOpt = '-q sx'
    # else:
    #   qOpt = '-q s'
    # qOpt = '-q sx' if size > 3.8*(10**9) else '-q s'
    qOpt = f'-q {self.__qmerge}'
    pOpt = f'-j {self.__nproc}' if self.__nproc > 1 else ''
    bOpt = f'-d {self.__buff_path}' if self.__buff_path is not None else ''
    cmd = shlex.split(f'bsub {qOpt} -o {self.__merge_log_path} '+
                      f'hadd -ff {pOpt} {bOpt}')
    cmd.append(self.__root_path)
    cmd.extend(self.__root_list)
    self.__merging_process = subprocess.run(cmd,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            check=True)
    time.sleep(1)
    buff = self.__merging_process.stdout
    job_id = bjob.BJob.read_job_id(buff.decode())
    if job_id is not None:
      btop_proc = subprocess.run(['btop', str(job_id)],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
      self.__merging_job = bjob.BJob(job_id)
      self.__dump_log('jid[merge]', self.__merging_job.get_job_id())
      self.__dump_log(None, '_'*80)
      self.__merge_status = 'RUNNING'
    else:
      logger.warning(f'failed to read jobid ({cmd})')
      self.__merge_status = 'FAILED'
  #____________________________________________________________________________
  def set_option(self, option):
    ''' Set option. '''
    self.__option = option

  #____________________________________________________________________________
  def update_status(self):
    ''' Update all status. '''
    self.__update_job_status()
    self.__update_merging_status()
    logger.debug(f'status is {self.__bjob_status}')
    if self.__bjob_status is None: # initial
      pass
    elif self.__bjob_status is 0: # running
      self.__status = 'BJOB-RUNNING'
    elif self.__bjob_status is True: # complete
      self.__status = 'BJOB-DONE'
    elif self.__bjob_status is False: # failure
      self.__status = 'FAILED'
      return
    elif self.__bjob_status == 1: # killed
      self.__status = 'TERMINATED'
      return
    else: # unknown
      self.__status = 'UNKNOWN'
    if self.__merge_status == 'INIT':
      pass
    elif self.__merge_status == 'RUNNING':
      self.__status = 'MERGING'
    elif self.__merge_status == 'DONE':
      self.__status = 'DONE'
      return
    elif self.__merge_status == 'FAILED':
      self.__status = 'FAILED'
      return
    elif self.__merge_status == 'TERMINATED':
      self.__status = 'TERMINATED'
      return
    else: # unknown
      self.__status = 'UNKNOWN'

  #____________________________________________________________________________
  def __dump_init_info(self):
    ''' Dump initial run info. '''
    self.__dump_log('start', time.ctime(self.__start_time))
    self.__dump_log(None, '_'*80)
    self.__dump_log('key', self.__key)
    self.__dump_log(None, '_'*80)
    self.__dump_log('bin', self.__bin_path)
    self.__dump_log('conf', self.__conf_path)
    self.__dump_log('unit', self.__div_unit)
    self.__dump_log('nevent', self.__nevents)
    self.__dump_log('out', self.__root_path)
    self.__dump_log('nproc', self.__nproc)
    self.__dump_log('buff_path', self.__buff_path)
    self.__dump_log('queue', self.__queue)
    self.__dump_log('qmerge', self.__qmerge)
    self.__dump_log('dirdummy', self.__dummy_dir.name)
    self.__dump_log(None, '_'*80)

  #____________________________________________________________________________
  def __dump_log(self, key=None, msg=None):
    ''' Dump log. '''
    buff = str()
    if key is None:
      if msg is not None:
        buff = str(msg)
    else:
      if msg is None:
        buff = str(key) + ':'
      elif isinstance(msg, list) or isinstance(msg, tuple):
        for i, item in enumerate(msg):
          buff += f'{key}({i}):'.ljust(16) + str(item)
          if i < len(msg) - 1:
            buff += '\n'
      else:
        buff = f'{key}:'.ljust(16) + str(msg)
    with open(self.__log_path, 'a') as flog:
      flog.write(buff + '\n')

  #____________________________________________________________________________
  def __make_element(self):
    ''' Make element list. '''
    if self.__div_unit <= 0:
      nsegs = 1
    else:
      nsegs = (1 if self.__nevents is None else
              (self.__nevents -1) // self.__div_unit + 1)
    for i in range(nsegs):
      self.__elem_list.append(f'{self.__basename}_{i:05d}')
    self.__dump_log('elem', self.__elem_list)
    self.__dump_log(None, '_'*80)
    return nsegs

  #____________________________________________________________________________
  def __make_log(self):
    ''' Make log. '''
    log_dir = os.path.join(SCRIPT_DIR, 'log')
    if not os.path.exists(log_dir):
      os.mkdir(log_dir)
    self.__log_path = os.path.join(log_dir, self.__basename+'.log')
    if os.path.exists(self.__log_path):
      os.remove(self.__log_path)

  #____________________________________________________________________________
  def __make_log_list(self):
    ''' Make log list. '''
    log_dir = os.path.join(SCRIPT_DIR, 'log')
    if not os.path.exists(log_dir):
      os.mkdir(log_dir)
    for item in self.__elem_list:
      path = os.path.join(log_dir, item+'.log')
      if os.path.exists(path):
        os.remove(path)
      self.__log_list.append(path)
    if 1 < len(self.__elem_list):
      self.__merge_log_path = os.path.join(log_dir,
                                          self.__basename+'_merge.log')
      if os.path.exists(self.__merge_log_path):
        os.remove(self.__merge_log_path)
    self.__dump_log('log', self.__log_list)
    self.__dump_log('mergelog', self.__merge_log_path)
    self.__dump_log(None, '_'*80)

  #____________________________________________________________________________
  def __make_macro(self):
    macro_path = os.path.join(self.__dummy_dir.name, 'tmp.mac')
    with open(macro_path, 'w') as f:
      f.write(f'/run/beamOn {min(self.__nevents, self.__div_unit)}')
    self.__macro_path = macro_path

  #____________________________________________________________________________
  def __make_path(self, path):
    if not os.path.exists(path):
      ret = os.path.join(os.path.dirname(self.__conf_path), path)
    ret = os.path.abspath(path)
    if os.path.exists(ret):
      return ret
    else:
      return path

  #____________________________________________________________________________
  def __make_root_list(self):
    ''' Make root list. '''
    for item in self.__elem_list:
      path = os.path.join(os.path.dirname(self.__root_path), item+'.root')
      self.__root_list.append(path)
    self.__dump_log('out', self.__root_list)
    self.__dump_log(None, '_'*80)

  #____________________________________________________________________________
  def __update_job_status(self):
    ''' Update job status. '''
    if (self.__bjob_status is not None and self.__bjob_status != 0):
      return
    n_complete = 0
    for i, job in enumerate(self.__bjob_list):
      stat = job.get_status()
      if stat is True:
        n_complete += 1
        if self.__jobstat_list[i] is 1:
          self.__dump_log(f'time[bjob({i})]',
                         self.decode_time(job.get_run_time()))
          self.__dump_log(None, '_'*80)
      elif stat is False and self.__jobstat_list[i] is 1:
        self.__bjob_status = False
        self.__dump_log('error', f'error at {job.get_tag()}')
        self.__dump_log(None, '_'*80)
      elif stat == 1 and self.__jobstat_list[i] != stat:
        jid = job.get_job_id()
        self.__jobid_list[i] = jid
        self.__dump_log(f'jid[bjob({i})]', jid)
        self.__dump_log(None, '_'*80)
      self.__jobstat_list[i] = stat
    if len(self.__elem_list) == n_complete:
      self.__bjob_status = True
    self.__ncomplete = n_complete

  #____________________________________________________________________________
  def __update_merging_status(self):
    ''' Update merging status. '''
    if self.__merge_status != 'RUNNING':
      return
    stat = self.__merging_job.get_status()
    if stat == 'INIT' or stat == 'PEND' or stat == 'RUN':
      return
    elif stat == 'DONE':
      self.__merge_status = 'DONE'
    elif stat == 'EXIT':
      self.__merge_status = 'FAILED'
      self.__dump_log('error', 'merging error')
      self.__dump_log(None, '_'*80)
    else:
      self.__merge_status = 'UNKNOWN'

  #____________________________________________________________________________
  @staticmethod
  def decode_status(data):
    ''' Decode status. '''
    buff = ''
    stat = data['stat']
    nseg = data['nseg']
    prog = data['prog']
    queue = data['queue']
    qmerge = data['qmerge'] if 'qmerge' in data else '?'
    if stat == 'INIT':
      buff = 'initialized'
    elif stat == 'PURGED':
      buff = 'purged'
    elif stat == 'BJOB-RUNNING':
      buff = f'running:{queue}({prog}/{nseg})'
    elif stat == 'BJOB-DONE':
      buff = f'running:{queue}({prog}/{nseg})'
    elif stat == 'MERGING':
      buff = f'merging:{qmerge}({nseg})'
    elif stat == 'TERMINATED':
      buff = 'terminated'
    elif stat == 'DONE':
      buff = f'done({nseg})'
    elif stat == 'FAILED':
      buff = 'failed'
    else:
      buff = 'unknown'
    return buff

  #____________________________________________________________________________
  @staticmethod
  def decode_time(data):
    ''' Decode time. '''
    second = 0
    if isinstance(data, dict):
      second = int(data['time'])
    elif isinstance(data, float):
      second = int(data)
    elif isinstance(data, int):
      second = int(data)
    hour    = second // 3600
    second -= hour    * 3600
    minute  = second // 60
    second -= minute  * 60
    return f'{hour}:{minute:02d}:{second:02d}'
