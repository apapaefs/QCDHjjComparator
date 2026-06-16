# QCDHjjComparator

Small diagnostic scripts for comparing matrix-element-squared values for
loop-induced `g g -> H j j` subprocesses between OpenLoops and MadGraph5_aMC.

The scripts assume prebuilt MadGraph standalone subprocess directories. They do
not generate or build MadGraph processes automatically.

## Contents

- `compare_openloops_mg_hjj.py`: generates shared `g g -> H + two partons`
  phase-space points and compares OpenLoops against MadGraph.
- `compare_openloops_mg_hjj_user_momenta.py`: compares one user-supplied
  momentum point.
- `hgg_mg_eval.f`: tiny Fortran driver for MadGraph standalone subprocesses.
- `mg_eval.mk`: makefile fragment that links `hgg_mg_eval.f` against the
  generated MadLoop objects and model libraries.
- `compare_openloops_mg_hjj.example.json`: example provider configuration.
- `user_momenta.example.dat`: example five-line momentum input.
- `cards/`: example MadGraph process cards.

## Generate MadGraph Processes

Use explicit subprocesses rather than `j j` if you want point-by-point
subprocess comparisons.

For `g g -> H g g`:

```bash
mg5_aMC cards/qcdhgg_proc_card.mg5
```

For `g g -> H u~ u`:

```bash
mg5_aMC cards/qcdhuux_proc_card.mg5
```

Edit the `output` paths in those cards before running them.

The important MadGraph syntax is:

```mg5
import model loop_sm-no_b_mass
generate g g > h g g [noborn=QCD] @1
```

and, for the quark channel:

```mg5
import model loop_sm-no_b_mass
generate g g > h u~ u [noborn=QCD] @1
```

## Build The MadGraph Evaluator

After MadGraph has generated a process, enter the subprocess directory, for
example:

```bash
cd /path/to/MG5_runs/QCDHgg/SubProcesses/PV1_0_1_gg_hgg
cp /path/to/QCDHjjComparator/hgg_mg_eval.f .
make -f Makefile -f /path/to/QCDHjjComparator/mg_eval.mk hgg_mg_eval
```

Do not use plain `make hgg_mg_eval`: the generated MadGraph Makefile does not
define that target, so Make falls back to compiling `hgg_mg_eval.f` by itself.
That produces undefined references to `ML5_...`, `SETPARA`, and
`UPDATE_AS_PARAM2`. The `mg_eval.mk` fragment adds the missing target and links
against the same generated MadLoop objects used by MadGraph's `check`
executable. It will also ask the generated `Source` makefile to build
`libdhelas.a` and `libmodel.a` if they are missing.

The evaluator reads `hgg_mg_points.dat` and writes `hgg_mg_eval.out`. It expects
the momentum order:

```text
incoming g, incoming g, H, parton, parton
```

The `hgg_mg_eval` name is only historical. The executable evaluates whatever
MadGraph subprocess directory it is linked in. If you copy the same driver into
`PV1_0_1_gg_huxu` and build it there, it evaluates `g g -> H u~ u`; if you copy
it into `PV1_0_1_gg_hgg`, it evaluates `g g -> H g g`.

## Configure Providers

Do not run `compare_openloops_mg_hjj.example.json` directly. It contains
placeholder `mg_dir` paths. Copy it first and replace the `mg_dir` entries with
the generated subprocess directories:

```bash
cp compare_openloops_mg_hjj.example.json compare_openloops_mg_hjj.json
```

For example:

```json
"mg_dir": "/path/to/MG5_runs/QCDHgg/SubProcesses/PV1_0_1_gg_hgg"
```

Also check `openloops_root`; it should point to an OpenLoops installation with
the relevant `pphjj2` process library installed.

## Compare Generated Points

```bash
python3 compare_openloops_mg_hjj.py \
  --config compare_openloops_mg_hjj.json \
  --subprocess hgg \
  --n 10 \
  --seed 260616 \
  --output hgg_ol_mg.tsv
```

The output columns include raw provider values, scaled values, the ratio
`ol_over_mg`, and the provider accuracy/return-code diagnostics.

## Compare A User-Supplied Point

Create a five-line file with columns:

```text
E px py pz
```

in the provider order `g, g, H, parton, parton`. See
`user_momenta.example.dat`.

Then run:

```bash
python3 compare_openloops_mg_hjj_user_momenta.py \
  --config compare_openloops_mg_hjj.json \
  --subprocess hgg \
  --momenta user_momenta.example.dat \
  --output user_point.tsv
```

## Tests

The unit tests do not require a real MadGraph build:

```bash
python3 -m unittest discover -s tests
python3 -m py_compile compare_openloops_mg_hjj.py compare_openloops_mg_hjj_user_momenta.py
```
