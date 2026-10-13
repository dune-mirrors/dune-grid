* Activate a venv that contains Dune

* Make sure sphinx is installed in that venv (by calling `which sphinx`).

* If it is not: `pip install sphinx`.

* For some strange reason, on my computer `sphinx-build` from the venv still uses
  the system `python3`, even though `python` is the correct one
  from the venv.

  If that happens for you, too, instead of calling `sphinx-build` directly,
  do it via

    python -m sphinx -b html . ../build/html/

  (to be called from the directory that contains `conf.py`).
