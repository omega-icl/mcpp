# Contributing to MC++

Thank you for considering a contribution.  MC++ is developed in the open at
<https://github.com/omega-icl/mcpp>, and issues, questions and pull requests are
all welcome.

## Reporting a problem

Open an issue at <https://github.com/omega-icl/mcpp/issues>.  A report is much
easier to act on when it says:

* the MC++ or `pymcpp` version (`pymcpp.__version__`, or the tag you built);
* the platform, compiler and Python version;
* whether you installed a PyPI wheel or built from source, and with which CMake
  options and interval backend;
* a small factorable function or script that reproduces the behaviour, with what
  you expected and what you got.

Bounds and relaxations have a particular failure mode worth naming: an enclosure
that is valid but wide is not a bug, whereas an enclosure that *excludes* a
value the function attains is.  Please say which one you are reporting.

## Asking a question

Usage questions belong in an issue too — there is no separate forum, and an
answered issue helps the next person.  The notebooks in `notebook/` and the
documentation at <https://pymcpp.readthedocs.io> cover each arithmetic and are
the fastest route to an answer.  For matters that do not belong in public, write
to the project lead (see the repository's Contacts section).

## Contributing code

1. Open an issue first for anything beyond a small fix, so the approach can be
   agreed before you spend time on it.
2. Fork the repository and branch from the default branch — whichever branch
   GitHub shows first, and the one `git clone` checks out for you.
3. Follow the surrounding style.  `.clang-format` and `.pre-commit-config.yaml`
   are in the repository: `pre-commit install` applies them as you commit.
4. Make the release checks pass before opening the pull request:

   ```
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --target check
   ```

   This gates the version written in three places, the recorded ABI
   fingerprint, the Python checks in `src/pymcpp/`, and a warning-free Doxygen
   build.  Install Doxygen to get the last of these; without it the step is
   skipped.
5. New behaviour needs a check that would fail without it.  Add it to
   `src/pymcpp/check_*.py` for anything reachable from Python, or to the
   relevant program under `test/` for C++-only code.
6. Add yourself to `AUTHORS`, in the group that fits your contribution, as part
   of the pull request.

## A new arithmetic

An arithmetic is a type that can be propagated through a DAG.  To add one,
specialise `mc::Op<T>` for it (`src/mc/mcop.hpp` states what the specialisation
must provide) and bind it in `src/pymcpp/` if it should be reachable from
Python.  Existing arithmetics — `interval.hpp`, `mccormick.hpp`, `scmodel.hpp`,
`supmodel.hpp` — are the models to follow.  Note that `pymcpp` exchanges types
with dependent modules such as CRONOS, so a change to a shared class layout
breaks their ABI; `INSTALL.md` ("Releasing") states the rule, and the
fingerprint in `tools/abi/` enforces it.

## Licence

Contributions are accepted under the Eclipse Public License 2.0 with the GNU
General Public License v2 or later as a Secondary License, the terms MC++ is
distributed under; see `LICENSE`.  By opening a pull request you confirm that
you may contribute the code under those terms.
