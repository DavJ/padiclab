# Preparing the first PyPI release

The package can be built and validated locally without a PyPI account:

```bash
python -m pip install '.[dev]' twine
python -m pytest
python -m build
python -m twine check dist/*
```

The first release uses the existing distribution name **`padic-tools`**,
import name **`padic_tools`**, and version **`0.1.0`**. The repository name
`padiclab` does not change these Python package names.

Before uploading, select a license. Add `LICENSE` and the corresponding
`project.license` metadata. No license
has been assigned on the owner's behalf. MIT and BSD-3-Clause are common
permissive options; the actual selection belongs to the owner.

Configure a pending Trusted Publisher on TestPyPI with:

| Setting | Value |
|---|---|
| Project | `padic-tools` |
| Owner | `DavJ` |
| Repository | `padiclab` |
| Workflow | `publish.yml` |
| Environment | `testpypi` |

Create the matching GitHub environment, using the desired release-review
settings. Tag a reviewed commit as `v0.1.0` after tests and mathematical
review, and manually dispatch the workflow **from that tag** to TestPyPI.
Verify the uploaded wheel in a fresh environment:

```bash
python -m venv /tmp/padic-testpypi-check
/tmp/padic-testpypi-check/bin/python -m pip install \
  --index-url https://test.pypi.org/simple/ --no-deps padic-tools==0.1.0
/tmp/padic-testpypi-check/bin/python -c \
  'from padic_tools import Qp; assert (Qp(5)("1/3")*3).is_congruent(1)'
```

For the real release, configure a separate publisher on PyPI with environment
`pypi`, then manually dispatch the same reviewed tag to PyPI. The default
target is TestPyPI. Ordinary pushes and pull requests do not publish.

The production publisher must use these exact values:

| Setting | Value |
|---|---|
| PyPI project | `padic-tools` |
| GitHub owner | `DavJ` |
| GitHub repository | `padiclab` |
| Workflow filename | `publish.yml` |
| Environment | `pypi` |

For a new project, add the pending publisher under the PyPI account's
**Publishing** settings. For an existing project, use the project's
**Publishing** settings. The workflow exchanges its GitHub Actions identity
for a short-lived publishing token; no long-lived PyPI token is needed.

After a successful production upload, verify the release from PyPI in a
separate environment:

```bash
python -m venv /tmp/padic-pypi-check
/tmp/padic-pypi-check/bin/python -m pip install \
  --index-url https://pypi.org/simple/ --no-deps padic-tools==0.1.0
/tmp/padic-pypi-check/bin/python -c \
  'from padic_tools import Qp; assert (Qp(5)("1/3")*3).is_congruent(1)'
```

Compare the wheel and source archive SHA-256 hashes shown by PyPI with the
files in the `release-distributions` artifact from the publishing run.

The manual release gate checks the matching version tag and the presence of
the owner-selected license. This is a project policy, not an assertion that
PyPI technically mandates an open-source license. No upload has been made.

See [PyPI publishing documentation](https://docs.pypi.org/publishing/publishing-to-pypi/)
and [Trusted Publishing setup](https://docs.pypi.org/trusted-publishers/using-a-publisher/).
