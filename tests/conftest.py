"""Shared pytest setup: stop PETSc from reporting pytest's own flags as unused options."""

import os

os.environ.setdefault("FIREDRAKE_DISABLE_OPTIONS_LEFT", "1")
