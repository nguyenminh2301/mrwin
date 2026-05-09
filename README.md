# mrwin

`mrwin` is a work-in-progress code repository for causal win-statistic methods.

The repository is being organized around an R package as the primary interface. The current Python implementation is retained only as a prototype and test harness while the R package code is developed.

## Current layout

```text
mrwin/
├── DESCRIPTION          R package metadata scaffold
├── NAMESPACE            R namespace scaffold
├── R/                   Future R package code
├── python/              Python prototype code used for validation
├── tests/
│   ├── python/          Python prototype tests
│   └── testthat/        Future R package tests
├── requirements.txt     Python prototype dependencies
├── pyproject.toml       Python prototype packaging/test config
└── run_all.sh           Optional Python prototype runner
```

## Python prototype

The Python code is not the intended long-term package interface. It is kept to preserve the current computational checks while the R implementation is written.

```bash
python -m pip install -r requirements.txt
python -m pytest
```

## R package scaffold

R package code should be added under `R/`. Tests for the R package should go under `tests/testthat/`.

The scaffold is intentionally minimal so manuscript drafts, replication outputs, and generated report files do not become part of the code repository.
