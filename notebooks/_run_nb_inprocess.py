"""Execute a notebook in-process (no ZMQ kernel; the sandbox forbids socket bind).
Captures stdout + rich display outputs into the .ipynb and fails loudly on any cell error."""
import sys, os, nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

path = sys.argv[1]
os.chdir(os.path.dirname(os.path.abspath(path)) or ".")
nb = nbformat.read(os.path.basename(path), as_version=4)
shell = InteractiveShell.instance()
shell.run_cell("import matplotlib; matplotlib.use('Agg')", store_history=False)

n = 0
for cell in nb.cells:
    if cell.cell_type != "code":
        continue
    n += 1
    with capture_output() as cap:
        res = shell.run_cell(cell.source, store_history=False)
    outs = []
    if cap.stdout:
        outs.append(nbformat.v4.new_output("stream", name="stdout", text=cap.stdout))
    if cap.stderr:
        outs.append(nbformat.v4.new_output("stream", name="stderr", text=cap.stderr))
    for ro in cap.outputs:
        outs.append(nbformat.v4.new_output("display_data", data=dict(ro.data),
                                           metadata=dict(ro.metadata or {})))
    cell.outputs = outs
    cell.execution_count = n
    if res.error_in_exec is not None:
        e = res.error_in_exec
        import traceback as tb
        outs.append(nbformat.v4.new_output("error", ename=type(e).__name__,
                    evalue=str(e), traceback=tb.format_exception(type(e), e, e.__traceback__)))
        nbformat.write(nb, os.path.basename(path))
        print(f"CELL {n} ERROR: {type(e).__name__}: {e}")
        sys.exit(1)

nbformat.write(nb, os.path.basename(path))
print(f"EXECUTED {n} code cells OK -> {os.path.basename(path)}")
