To run Python script with uv:
```bash
uv run python slift.py`
```
To lauch Jupyter Lab:
```bash
uv run --with jupyter jupyter lab`
```
To install gurobipy into the current uv environment without adding it to pyproject.toml:
```bash
uv pip install gurobipy
```
or from inside your uv environment/project, install it with:
```bash
uv add gurobipy
```
