# Two-Sensor Synthetic Fixture

This fixture follows a common test-data layout:
- `tests/fixtures/two_sensor_case/input/`: deterministic synthetic inputs
- `tests/fixtures/two_sensor_case/expected/`: checked expected outputs used by regression tests

## Contents
- `generate_two_sensor_case.py`: deterministic generator for synthetic LAZ + sensor table
- `case_config.json`: shared case metadata (plot bounds, voxel size, filename parsing indices)
- `input/laz/*.laz`: one LAZ per static sensor
- `input/sensor_positions.txt`: scan-position table for `OccPy.define_sensor_pos(...)`
- `expected/*.npy`: expected `Nhit/Nmiss/Nocc/Classification` arrays (to be generated once and committed)

## Generate Input
Run from repository root:

```bash
python tests/fixtures/two_sensor_case/generate_two_sensor_case.py
```

## Produce Expected Output Once
Use your normal OccPy workflow on this fixture input and write:
- `expected/Nhit.npy`
- `expected/Nmiss.npy`
- `expected/Nocc.npy`
- `expected/Classification.npy`

The regression test will compare current outputs against these files.
