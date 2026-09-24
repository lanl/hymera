#!/usr/bin/env bash
# Provision the Python environment the regression harness needs.
#
# The harness reads the solver's binary and ASCII output and plots it, so it
# needs numpy and matplotlib. Neither is linked into the solver -- this is
# tooling only, so a virtualenv is sufficient and avoids touching the Spack
# environment that builds the code.
#
# Idempotent: re-running is cheap and safe. Source-able for the VENV path.
#
#   ./tests/regression/setup_env.sh          # create/verify
#   . ./tests/regression/setup_env.sh        # ...and export HARNESS_PYTHON
set -euo pipefail

VENV="${HYMERA_HARNESS_VENV:-$HOME/.venv-hymera}"
PYTHON="$VENV/bin/python"

if [ ! -x "$PYTHON" ]; then
  echo "creating harness venv at $VENV"
  python3 -m venv "$VENV"
fi

if ! "$PYTHON" -c 'import numpy, matplotlib' 2>/dev/null; then
  echo "installing numpy and matplotlib"
  "$VENV/bin/pip" install --quiet --upgrade pip
  "$VENV/bin/pip" install --quiet numpy matplotlib
fi

"$PYTHON" - <<'EOF'
import matplotlib
matplotlib.use("Agg")
import numpy
print(f"harness python ready: numpy {numpy.__version__}, "
      f"matplotlib {matplotlib.__version__}")
EOF

export HARNESS_PYTHON="$PYTHON"
echo "HARNESS_PYTHON=$PYTHON"
