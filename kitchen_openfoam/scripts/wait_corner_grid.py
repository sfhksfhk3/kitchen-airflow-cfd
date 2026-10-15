"""Prepare the independent corner mesh sensitivity once its saved state exists."""
from pathlib import Path
import time
from foam_io import times
from prepare_sensitivity import prepare
from run_batch import run_case
ROOT=Path(__file__).resolve().parents[1]
while not any(90<=t<=110 for t in times(ROOT/'cases/corner_both_closed')):time.sleep(5)
p=prepare('corner','both_closed','grid');print(run_case(p),flush=True)
