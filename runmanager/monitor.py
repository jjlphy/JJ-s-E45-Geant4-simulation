#!/usr/bin/env python3

import argparse
import fcntl
import json
import logging
import logging.config
import os
import sys
import time
import yaml

top_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(top_dir, 'module'))

import runmanager
import singlerun
import utility
from utility import pycolor as cl

logger = logging.getLogger(__name__)
DISPLAY_PERIOD = 5  # second
parsed = None

#______________________________________________________________________________
def display(filename):
  buff = str()
  with open(filename, 'r') as f:
    try:
      fcntl.flock(f.fileno(), fcntl.LOCK_SH)
    except IOError as e:
      logger.error(e)
      return False
    else:
      buff = f.read()
    finally:
      fcntl.flock(f.fileno(), fcntl.LOCK_UN)
  try:
    info = json.loads(buff)
  except ValueError as e:
    logger.error(e)
    return False
  os.system('clear')
  buff = (cl.reverce + cl.bold
          + 'KEY'.ljust(8) + '  '
          + 'STAT'.ljust(18) + '  '
          + 'BIN'.ljust(16) + '  '
          + 'CONF'.ljust(16) + '  '
          + '#EVENT(UNIT)'.ljust(24) + '  '
          + 'ROOT'.ljust(16) + '  '
          + 'TIME'.ljust(8)
          + cl.end)
  print(buff)
  n_unfinished = 0
  for key, item in sorted(info.items(),
                          key=lambda x:int(x[0]) if isinstance(x[0], int) else x[0]):
    status = singlerun.SingleRun.decode_status(item)
    if parsed.skip and 'done' in status:
      continue
    n_unfinished += 1
    ptime = singlerun.SingleRun.decode_time(item)
    nev = None
    if 'nev' in item:
      nev = item['nev']
    if 'unit' in item:
      unit = item['unit']
    inbuf = f'{str(nev)} ({str(unit)})'
    buff = (cl.bold + key[:8].ljust(8) + cl.end + '  '
            + f'{cl.reverce}{cl.bold}{cl.red}{status}{cl.end}'.ljust(16 + 20)
            + '  '
            + os.path.basename(item['bin'])[-16:].ljust(16) + '  '
            + os.path.basename(item['conf'])[-16:].ljust(16) + '  '
            + inbuf[-24:].ljust(24)
            + '  '
            + os.path.basename(item['root'])[-16:].ljust(16) + '  '
            + ptime.rjust(8))
    print(buff)
  # if n_unfinished == 0:
  #   print('finished')
  #   return False
  # else:
  return True
  # buff = cl.reverce + cl.bold + "Press 'Ctrl-C' to exit" + cl.end
  # print(buff)

#______________________________________________________________________________
def main(path):
  ptime = time.time()
  try:
    while display(path):
      dtime = DISPLAY_PERIOD - (time.time() - ptime)
      if dtime > 0:
        time.sleep(dtime)
      ptime = time.time()
  except KeyboardInterrupt as e:
    logger.info(e)
  except FileNotFoundError as e:
    logger.error(e)

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-s', '--skip', help='skip runs with status of "done"',
                      action='store_true')
  parser.add_argument('stat_json', help='stat json file path')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(top_dir, 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  main(parsed.stat_json)
