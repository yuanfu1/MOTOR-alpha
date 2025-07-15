(test_harness)=

# PETSc Testing System

The PETSc test system consists of

- Formatted comments at the bottom of the tutorials and test source files that describes the tests to be run.
- The *test generator* (`config/gmakegentest.py`) that parses the tutorial and test source files and generates the makefiles and shell scripts. This is run
  automatically by the make system and rarely is run directly.
- The *PETSc test harness* that consists of makefile and shell scripts that runs the executables with several logging and reporting features.

Details on using the harness may be found in the {ref}`user's manual <sec_runningtests>`. The testing system is used by {any}`pipelines`.

## PETSc Test Description Language

PETSc tests and tutorials contain at the bottom of the their source files a simple language to
describe tests and subtests required to run executables associated with
compilation of that file. The general skeleton of the file is

```
static const char help[] = "A simple MOAB example\n";

...
<source code>
...

/*TEST
   build:
     requires: moab
   testset:
     suffix: 1
     requires: !complex
   testset:
     suffix: 2
     args: -debug -fields v1,v2,v3
     test:
     test:
       args: -foo bar
TEST*/
```

For our language, a *test* is associated with the following

- A single shell script

- A single makefile

- An output file that represents the *expected results*. It is also possible -- though unusual -- to have multiple output files for a single test

- Two or more command tests, usually:

  - one or more `mpiexec` tests that run the executable
  - one or more `diff` tests to compare output with the expected result

Our language also supports a *testset* that specifies either a new test
entirely or multiple executable/diff tests within a single test. At the
core, the executable/diff test combination will look something like
this:

```sh
mpiexec -n 1 ../ex1 1> ex1.tmp 2> ex1.err
diff ex1.tmp output/ex1.out 1> diff-ex1.tmp 2> diff-ex1.err
```

In practice, we want to do various logging and counting by the test
harness; as are explained further below. The input language supports
simple yet flexible test control.

(test_harness_data)=

### Datafiles needed for some tests

Some tests require matrices or meshes that are too large for the primary PETSc Git repository.
The repository [datafiles](https://gitlab.com/petsc/datafiles) contains all the test files needed for the test suite.
To run these tests one must first clone the datafiles repository and then set the environmental variable `DATAFILESPATH`.
For these tests `requires: datafilespath` should be specified.

### Runtime Language Options

At the end of each test file, a marked comment block is
inserted to describe the test(s) to be run. The elements of the test are
done with a set of supported key words that sets up the test.

The goals of the language are to be

- as minimal as possible with the simplest test requiring only one keyword,
- independent of the filename such that a file can be renamed without rewriting the tests, and
- intuitive.

In order to enable the second goal, the *basestring* of the filename is
defined as the filename without the extension; for example, if the
filename is `ex1.c`, then `basestring=ex1`.

With this background, these keywords are as follows.

- **testset** or **test**: (*Required*)

  - At the top level either a single test or a test set must be
    specified. All other keywords are sub-entries of this keyword.

- **suffix**: (*Optional*; *Default:* `suffix=""`)

  - The test name is given by `testname = basestring` if the suffix
    is set to an empty string, and by
    `testname = basestring + "_" + suffix` otherwise.
  - This can be specified only for top level test nodes.

- **output_file**: (*Optional*; *Default:*
  `output_file = "output/" + testname + ".out"`)

  - The output of the test is to be compared with an *expected result*
    whose name is given by `output_file`.
  - This file is described relative to the source directory of the
    source file and should be in the output subdirectory (for example,
    `output/ex1.out`)

- **nsize**: (*Optional*; *Default:* `nsize=1`)

  - This integer is passed to mpiexec; i.e., `mpiexec -n nsize`

- **args**: (*Optional*; *Default:* `""`)

  - These arguments are passed to the executable.

- **diff_args**: (*Optional*; *Default:* `""`)

  - These arguments are passed to the `lib/petsc/bin/petscdiff` script that
    is used in the diff part of the test. For example, `-j` enables testing
    the floating point numbers.

- **TODO**: (*Optional*; *Default:* `False`)

  - Setting this Boolean to True will tell the test to appear in the
    test harness but report only TODO per the TAP standard. Optionally
    provide a string indicating why it is todo.
  - A runscript will be generated and can easily be modified by hand
    to run.

- **filter**: (*Optional*; *Default:* `""`)

  - Sometimes only a subset of the output is meant to be tested
    against the expected result. If this keyword is used, it filters
    the executable output to
    compare with `output_file`.
  - The value of this is the command to be run, for example,
    `grep foo` or `sort -nr`.
  - **NOTE: this method of testing error output is NOT recommended. See section on**
    {ref}`testing errors <sec_testing_error_testing>` **instead.** If the filter begins
    with `Error:`, then the test is assumed to be testing the `stderr` output, and the
    error code and output are set up to be tested.

- **filter_output**: (*Optional*; *Default:* `""`)

  - Sometimes filtering the output file is useful for standardizing
    tests. For example, in order to handle the issues related to
    parallel output, both the output from the test example and the
    output file need to be sorted (since sort does not produce the
    same output on all machines). This works the same as filter to
    implement this feature

- **localrunfiles**: (*Optional*; *Default:* `""`)

  - Some tests
    require runtime files that are maintained in the source tree.
    Files in this (space-delimited) list will be copied over to the
    testing directory so they will be found by the executable. If you
    list a directory instead of files, it will copy the entire
    directory (this is limited currently to a single directory)
  - The copying is done by the test generator and not by creating
    makefile dependencies.

- **temporaries**: (*Optional*; *Default:* `""`)

  - Some tests produce temporary files that are read by the filter
    to compare to expected results.
    Files in this (space-delimited) list will cleared before
    the test is run to ensure that stale temporary files are not read.

- **requires**: (*Optional*; *Default:* `""`)

  - This is a space-delimited list of run requirements (not build
    requirements; see Build requirements below).
  - In general, the language supports `and` and `not` constructs
    using `! => not` and `, => and`.
  - MPIUNI should work for all -n 1 examples so this need not be in
    the requirements list.
  - Some tests require matrices or meshes contained in the
    directory given by the environmental variable `DATAFILESPATH`.
    For these tests `requires: datafilespath` is
    specified. See {any}`test harness data<test_harness_data>`
  - Packages are indicated with lower-case specification, for example,
    `requires: superlu_dist`.
  - Any defined variable in petscconf.h can be specified with the
    `defined(...)` syntax, for example, `defined(PETSC_USE_INFO)`.
  - Any definition of the form `PETSC_HAVE_FOO` can just use
    `requires: foo` similar to how third-party packages are handled.

- **timeoutfactor**: (*Optional*; *Default:* `"1"`)

  - This parameter allows you to extend the default timeout for an
    individual test such that the new timeout time is
    `timeout=(default timeout) x (timeoutfactor)`.
  - Tests are limited to a set time that is found at the top of
    `"config/petsc_harness.sh"` and can be overwritten by passing in
    the `TIMEOUT` argument to `gmakefile`

- **env**: (*Optional*; *Default:* `env=""`)

  - Allows you to set environment variables for the test. Values are copied verbatim to
    the runscript and defined and exported prior to all other variables.

  - Variables defined within `env:` blocks are expanded and processed by the shell that
    runs the runscript. No prior preprocessing (other than splitting the lines into
    separate declarations) is done. This means that any escaping of special characters
    must be done in the text of the `TEST` block.

  - Defining the `env:` keyword more than once is allowed. Subsequent declarations are
    then appended to prior list of declarations . Multiple environment variables may also
    be defined in the same `env:` block, i.e. given a test `ex1.c` with the following
    spec:

    ```yaml
    test:
      env: FOO=1 BAR=1

    # equivalently
    test:
      env: FOO=1
      env: BAR=1
    ```

    results in

    ```console
    $ export FOO=1; export BAR=1; ./ex1
    ```

  - Variables defined in an `env:` block are evaluated by the runscript in the order in
    which they are defined in the `TEST` block. Thus it is possible for later variables
    to refer to previously defined ones:

    ```yaml
    test:
      env: FOO='hello' BAR=${FOO}
    ```

    results in

    ```console
    $ export FOO='hello'; export BAR=${FOO}; ./ex1
    # expanded by shell to
    $ export FOO='hello'; export BAR='hello'; ./ex1
    ```

    Note this also implies that

    ```yaml
    test:
      env: FOO=1 FOO=0
    ```

    results in

    ```console
    $ export FOO=1; export FOO=0; ./ex1
    ```

### Additional Specifications

In addition to the above keywords, other language features are
supported.

- **for loops**: Specifying `{{list of values}}` will generate a loop over
  an enclosed space-delimited list of values.
  It is supported within `nsize` and `args`. For example,

  ```
  nsize: {{1 2 4}}
  args: -matload_block_size {{2 3}shared output}
  ```

  Here the output for each `-matload_block_size` value is assumed to be
  the same so that only one output file is needed.

  If the loop causes different output for each loop iteration, then `separate output` needs to be used:

  ```
  args: -matload_block_size {{2 3}separate output}
  ```

  In this case, each loop value generates a separate script,
  and uses a separate output file for comparison.

  Note that `{{...}}` is equivalent to `{{...}shared output}`.

(sec_testing_error_testing)=

### Testing Errors And Exceptional Code

It is possible (and encouraged!) to test error conditions within the test harness. Since
error messages produced by `SETERRQ()` and friends are not portable between systems,
additional arguments must be passed to tests to modify error handling, specifically:

```yaml
args: -petsc_ci_portable_error_output -error_output_stdout
```

These arguments have the following effect:

- `-petsc_ci_portable_error_output`: Strips system or configuration-specific information
  from error messages. Specifically this:

  - Removes all path components except the file name from the traceback
  - Removes line and column numbers from the traceback
  - Removes PETSc version information
  - Removes `configure` options used
  - Removes system name
  - Removes hostname
  - Removes date

  With this option error messages will be identical across systems, runs, and PETSc
  configurations (barring of course configurations in which the error is not raised).

  Furthermore, this option also changes the default behavior of the error handler to
  **gracefully** exit where possible. For single-ranked runs this means returning with
  exit-code `0` and calling `MPI_Finalize()` instead of `MPI_Abort()`. Multi-rank
  tests will call `MPI_Abort()` on errors raised on `PETSC_COMM_SELF`, but will call
  `MPI_Finalize()` otherwise.

- `-error_output_stdout`: Forces `SETERRQ()` and friends to dump error messages to
  `stdout` instead of `stderr`. While using `stderr` (alongside the `Error:`
  sub-directive under `filter:`) also works it appears to be unstable under heavy
  load, especially in CI.

Using both options in tandem allows one to use the normal `output:` mechanism to compare
expected and actual error outputs.

When writing ASCII output that may be not portable, so one wants `-petsc_ci_portable_error_output` to
cause the output to be skipped, enclose the output with code such as

```
if (!PetscCIEnabledPortableErrorOutput)
```

to prevent it from being output when the CI test harness is running.

### Test Block Examples

The following is the simplest test block:

```yaml
/*TEST
  test:
TEST*/
```

If this block is in `src/a/b/examples/tutorials/ex1.c`, then it will
create `a_b_tutorials-ex1` test that requires only one
process, with no arguments, and diff the resultant output with
`src/a/b/examples/tutorials/output/ex1.out`.

For Fortran, the equivalent is

```fortran
!/*TEST
!  test:
!TEST*/
```

A more complete example, showing just the lines between `/*TEST` and `TEST*/`:

```yaml
test:
test:
  suffix: 1
  nsize: 2
  args: -t 2 -pc_type jacobi -ksp_monitor_short -ksp_type gmres
  args: -ksp_gmres_cgs_refinement_type refine_always -s2_ksp_type bcgs
  args: -s2_pc_type jacobi -s2_ksp_monitor_short
  requires: x
```

This creates two tests. Assuming that this is
`src/a/b/examples/tutorials/ex1.c`, the tests would be
`a_b_tutorials-ex1` and `a_b_tutorials-ex1_1`.

Following is an example of how to test a permutation of arguments
against the same output file:

```yaml
testset:
  suffix: 19
  requires: datafilespath
  args: -f0 ${DATAFILESPATH}/matrices/poisson1
  args: -ksp_type cg -pc_type icc -pc_factor_levels 2
  test:
  test:
    args: -mat_type seqsbaij
```

Assuming that this is `ex10.c`, there would be two mpiexec/diff
invocations in `runex10_19.sh`.

Here is a similar example, but the permutation of arguments creates
different output:

```yaml
testset:
  requires: datafilespath
  args: -f0 ${DATAFILESPATH}/matrices/medium
  args: -ksp_type bicg
  test:
    suffix: 4
    args: -pc_type lu
  test:
    suffix: 5
```

Assuming that this is `ex10.c`, two shell scripts will be created:
`runex10_4.sh` and `runex10_5.sh`.

An example using a for loop is:

```yaml
testset:
  suffix: 1
  args: -f ${DATAFILESPATH}/matrices/small -mat_type aij
  requires: datafilespath
testset:
  suffix: 2
  output_file: output/ex138_1.out
  args: -f ${DATAFILESPATH}/matrices/small
  args: -mat_type baij -matload_block_size {{2 3}shared output}
  requires: datafilespath
```

In this example, `ex138_2` will invoke `runex138_2.sh` twice with
two different arguments, but both are diffed with the same file.

Following is an example showing the hierarchical nature of the test
specification.

```yaml
testset:
  suffix:2
  output_file: output/ex138_1.out
  args: -f ${DATAFILESPATH}/matrices/small -mat_type baij
  test:
    args: -matload_block_size 2
  test:
    args: -matload_block_size 3
```

This is functionally equivalent to the for loop shown above.

Here is a more complex example using for loops:

```yaml
testset:
  suffix: 19
  requires: datafilespath
  args: -f0 ${DATAFILESPATH}/matrices/poisson1
  args: -ksp_type cg -pc_type icc
  args: -pc_factor_levels {{0 2 4}separate output}
  test:
  test:
    args: -mat_type seqsbaij
```

If this is in `ex10.c`, then the shell scripts generated would be

- `runex10_19_pc_factor_levels-0.sh`
- `runex10_19_pc_factor_levels-2.sh`
- `runex10_19_pc_factor_levels-4.sh`

Each shell script would invoke twice.

### Build Language Options

You can specify issues related to the compilation of the source file
with the `build:` block. The language is as follows.

- **requires:** (*Optional*; *Default:* `""`)

  - Same as the runtime requirements (for example, can include
    `requires: fftw`) but also requirements related to types:

    1. Precision types: `single`, `double`, `quad`, `int32`
    2. Scalar types: `complex` (and `!complex`)

  - In addition, `TODO` is available to allow you to skip the build
    of this file but still maintain it in the source tree.

- **depends:** (*Optional*; *Default:* `""`)

  - List any dependencies required to compile the file

A typical example for compiling for only real numbers is

```
/*TEST
  build:
    requires: !complex
  test:
TEST*/
```

## Running the tests

The make rules for running tests are contained in `gmakefile.test` in the PETSc root directory. They can usually be accessed by
simply using commands such as

```console
$ make test
```

or, for a list of test options,

```console
$ make help-test
```

### Determining the failed jobs of a given run

The running of the test harness will show which tests fail, but you may not have
logged the output or run without showing the full error. The best way of
examining the errors is with this command:

```console
$ $EDITOR $PETSC_DIR/$PETSC_ARCH/tests/test*err.log
```

This method can also be used for the PETSc continuous integration (CI) pipeline jobs. For failed jobs you can download the
log files from the `artifacts download` tab on the right side:

:::{figure} /images/developers/test-artifacts.png
:alt: Test Artifacts at Gitlab

Test artifacts can be downloaded from GitLab.
:::

To see the list of all tests that failed from the last run, you can also run this command:

```console
$ make print-test test-fail=1
```

To print it out in a column format:

```console
$ make print-test test-fail=1 | tr ' ' '\n' | sort
```

Once you know which tests failed, the question is how to debug them.

### Introduction to debugging workflows

Here, two different workflows on developing with the test harness are presented,
and then the language for adding a new test is described. Before describing the
workflow, we first discuss the output of the test harness and how it maps onto
makefile targets and shell scripts.

Consider this line from running the PETSc test system:

```
TEST arch-ci-linux-uni-pkgs/tests/counts/vec_is_sf_tests-ex1_basic_1.counts
```

The string `vec_is_sf_tests-ex1_basic_1` gives the following information:

- The file generating the tests is found in `$PETSC_DIR/src/vec/is/sf/tests/ex1.c`
- The makefile target for the *test* is `vec_is_sf_tests-ex1_basic_1`
- The makefile target for the *executable* is `$PETSC_ARCH/tests/vec/is/sf/tests/ex1`
- The shell script running the test is located at: `$PETSC_DIR/$PETSC_ARCH/tests/vec/is/sf/tests/runex1_basic_1.sh`

Let's say that you want to debug a single test as part of development. There
are two basic methods of doing this: 1) use shell script directly in test
directory, or 2) use the gmakefile.test from the top level directory. We present both
workflows.

### Debugging a test using shell the generated scripts

First, look at the working directory and the options for the
scripts:

```console
$ cd $PETSC_ARCH/tests/vec/is/sf/tests
$ ./runex1_basic_1.sh -h
Usage: ./runex1_basic_1.sh [options]

OPTIONS
  -a <args> ......... Override default arguments
  -c ................ Cleanup (remove generated files)
  -C ................ Compile
  -d ................ Launch in debugger
  -e <args> ......... Add extra arguments to default
  -f ................ force attempt to run test that would otherwise be skipped
  -h ................ help: print this message
  -n <integer> ...... Override the number of processors to use
  -j ................ Pass -j to petscdiff (just use diff)
  -J <arg> .......... Pass -J to petscdiff (just use diff with arg)
  -m ................ Update results using petscdiff
  -M ................ Update alt files using petscdiff
  -o <arg> .......... Output format: 'interactive', 'err_only'
  -p ................ Print command: Print first command and exit
  -t ................ Override the default timeout (default=60 sec)
  -U ................ run cUda-memcheck
  -V ................ run Valgrind
  -v ................ Verbose: Print commands
```

We will be using the `-C`, `-V`, and `-p` flags.

A basic workflow is something similar to:

```console
$ <edit>
$ runex1_basic_1.sh -C
$ <edit>
$ ...
$ runex1_basic_1.sh -m # If need to update results
$ ...
$ runex1_basic_1.sh -V # Make sure valgrind clean
$ cd $PETSC_DIR
$ git commit -a
```

For loops it sometimes can become onerous to run the whole test.
In this case, you can use the `-p` flag to print just the first
command. It will print a command suitable for running from
`$PETSC_DIR`, but it is easy to modify for execution in the test
directory:

```console
$ runex1_basic_1.sh -p
```

### Debugging a PETSc test using the gmakefile.test

First recall how to find help for the options:

```console
$ make help-test
Test usage:
   /usr/bin/gmake --no-print-directory test <options>

Options:
  NO_RM=1           Do not remove the executables after running
  REPLACE=1         Replace the output in PETSC_DIR source tree (-m to test scripts)
  OUTPUT=1          Show only the errors on stdout
  ALT=1             Replace 'alt' output in PETSC_DIR source tree (-M to test scripts)
  DIFF_NUMBERS=1    Diff the numbers in the output (-j to test scripts and petscdiff)
  CUDAMEMCHECK=1    Execute the tests using cuda-memcheck (-U to test scripts)
                    Use PETSC_CUDAMEMCHECK_COMMAND to change the executable to run and
                    PETSC_CUDAMEMCHECK_ARGS to change the arguments (note: both
                    cuda-memcheck and compute-sanitizer are supported)
  VALGRIND=1        Execute the tests using valgrind (-V to test scripts)
  DEBUG=1           Launch tests in the debugger (-d to the scripts)
  NP=<num proc>     Set a number of processors to pass to scripts.
  FORCE=1           Force SKIP or TODO tests to run
  PRINTONLY=1       Print the command, but do not run.  For loops print first command
  TIMEOUT=<time>    Test timeout limit in seconds (default in config/petsc_harness.sh)
  TESTDIR='tests'   Subdirectory where tests are run ($PETSC_DIR/$PETSC_ARCH
                    or /
                    or /share/petsc/examples/)
  TESTBASE='tests'   Subdirectory where tests are run ($PETSC_DIR/$PETSC_ARCH)
  OPTIONS='<args>'  Override options to scripts (-a to test scripts)
  EXTRA_OPTIONS='<args>'  Add options to scripts (-e to test scripts)

Special options for macOS:
  MACOS_FIREWALL=1  Add each built test to the macOS firewall list to prevent popups. Configure --with-macos-firewall-rules to make this default

Tests can be generated by searching with multiple methods
  For general searching (using config/query_tests.py):
    /usr/bin/gmake --no-print-directory test search='sys*ex2*'
   or the shortcut using s
    /usr/bin/gmake --no-print-directory test s='sys*ex2*'
  You can also use the full path to a file directory
    /usr/bin/gmake --no-print-directory test s='src/sys/tests/'
   or a file
    /usr/bin/gmake --no-print-directory test s='src/sys/tests/ex1.c'

  To search for fields from the original test definitions:
    /usr/bin/gmake --no-print-directory test query='requires' queryval='*MPI_PROCESS_SHARED_MEMORY*'
   or the shortcut using q and qv
    /usr/bin/gmake --no-print-directory test q='requires' qv='*MPI_PROCESS_SHARED_MEMORY*'
  To filter results from other searches, use searchin
    /usr/bin/gmake --no-print-directory test s='src/sys/tests/' searchin='*options*'

  To re-run the last tests which failed:
    /usr/bin/gmake --no-print-directory test test-fail='1'

  To see which targets match a given pattern (useful for doing a specific target):
    /usr/bin/gmake --no-print-directory print-test search=sys*

  To build an executable, give full path to location:
    /usr/bin/gmake --no-print-directory ${PETSC_ARCH}/tests/sys/tests/ex1
  or make the test with NO_RM=1
```

To compile the test and run it:

```console
$ make test search=vec_is_sf_tests-ex1_basic_1
```

This can consist of your basic workflow. However,
for the normal compile and edit, running the entire harness with search can be
cumbersome. So first get the command:

```console
$ make vec_is_sf_tests-ex1_basic_1 PRINTONLY=1
<copy command>
<edit>
$ make $PETSC_ARCH/tests/vec/is/sf/tests/ex1
$ /scratch/kruger/contrib/petsc-mpich-cxx/bin/mpiexec -n 1 arch-mpich-cxx-py3/tests/vec/is/sf/tests/ex1
...
$ cd $PETSC_DIR
$ git commit -a
```

### Advanced searching

For forming a search, it is recommended to always use `print-test` instead of
`test` to make sure it is returning the values that you want.

The three basic and recommended arguments are:

- `search` (or `s`)

  - Searches based on name of test target (see above)

  - Use the familiar glob syntax (like the Unix `ls` command). Example:

    ```console
    $ make print-test search='vec_is*ex1*basic*1'
    ```

    Equivalently:

    ```console
    $ make print-test s='vec_is*ex1*basic*1'
    ```

  - It also takes full paths. Examples:

    ```console
    $ make print-test s='src/vec/is/tests/ex1.c'
    ```

    ```console
    $ make print-test s='src/dm/impls/plex/tests/'
    ```

    ```console
    $ make print-test s='src/dm/impls/plex/tests/ex1.c'
    ```

- `query` and `queryval` (or `q` and `qv`)

  - `query` corresponds to test harness keyword, `queryval` to the value. Example:

    ```console
    $ make print-test query='suffix' queryval='basic_1'
    ```

  - Invokes `config/query_tests.py` to query the tests (see
    `config/query_tests.py --help` for more information).

  - See below for how to use as it has many features

- `searchin` (or `i`)

  - Filters results of above searches. Example:

    ```console
    $ make print-test s='src/dm/impls/plex/tests/ex1.c' i='*refine_overlap_2d*'
    ```

Searching using GNU make's native regexp functionality is kept for people who like it, but most developers will likely prefer the above methods:

- `gmakesearch`

  - Use GNU make's own filter capability.

  - Fast, but requires knowing GNU make regex syntax which uses `%` instead of `*`

  - Also very limited (cannot use two `%`'s for example)

  - Example:

    ```console
    $ make test gmakesearch='vec_is%ex1_basic_1'
    ```

- `gmakesearchin`

  - Use GNU make's own filter capability to search in previous results. Example:

    ```console
    $ make test gmakesearch='vec_is%1' gmakesearchin='basic'
    ```

### Query-based searching

Note the use of glob style matching is also accepted in the value field:

```console
$ make print-test query='suffix' queryval='basic_1'
```

```console
$ make print-test query='requires' queryval='cuda'
```

```console
$ make print-test query='requires' queryval='defined(PETSC_HAVE_MPI_GPU_AWARE)'
```

```console
$ make print-test query='requires' queryval='*GPU_AWARE*'
```

Using the `name` field is equivalent to the search above:

- Example:

  ```console
  $ make print-test query='name' queryval='vec_is*ex1*basic*1'
  ```

- This can be combined with union/intersect queries as discussed below

Arguments are tricky to search for. Consider

```none
args: -ksp_monitor_short -pc_type ml -ksp_max_it 3
```

Search terms are

```none
ksp_monitor, pc_type ml, ksp_max_it
```

Certain items are ignored:

- Numbers (see `ksp_max_it` above), but floats are ignored as well.
- Loops: `args: -pc_fieldsplit_diag_use_amat {{0 1}}` gives `pc_fieldsplit_diag_use_amat` as the search term
- Input files: `-f *`

Examples of argument searching:

```console
$ make print-test query='args' queryval='ksp_monitor'
```

```console
$ make print-test query='args' queryval='*monitor*'
```

```console
$ make print-test query='args' queryval='pc_type ml'
```

Multiple simultaneous queries can be performed with union (`,`), and intersection
(`|`) operators in the `query` field. One may also use their alternate spellings
(`%AND%` and `%OR%` respectively). The alternate spellings are useful in cases where
one cannot avoid (possibly multiple) shell expansions that might otherwise interpret the
`|` operator as a shell pipe. Examples:

- All examples using `cuda` and all examples using `hip`:

  ```console
  $ make print-test query='requires,requires' queryval='cuda,hip'
  # equivalently
  $ make print-test query='requires%AND%requires' queryval='cuda%AND%hip'
  ```

- Examples that require both triangle and ctetgen (intersection of tests)

  ```console
  $ make print-test query='requires|requires' queryval='ctetgen,triangle'
  # equivalently
  $ make print-test query='requires%OR%requires' queryval='ctetgen%AND%triangle'
  ```

- Tests that require either `ctetgen` or `triangle`

  ```console
  $ make print-test query='requires,requires' queryval='ctetgen,triangle'
  # equivalently
  $ make print-test query='requires%AND%requires' queryval='ctetgen%AND%triangle'
  ```

- Find `cuda` examples in the `dm` package.

  ```console
  $ make print-test query='requires|name' queryval='cuda,dm*'
  # equivalently
  $ make print-test query='requires%OR%name' queryval='cuda%AND%dm*'
  ```

Here is a way of getting a feel for how the union and intersect operators work:

```console
$ make print-test query='requires' queryval='ctetgen' | tr ' ' '\n' | wc -l
170
$ make print-test query='requires' queryval='triangle' | tr ' ' '\n' | wc -l
330
$ make print-test query='requires,requires' queryval='ctetgen,triangle' | tr ' ' '\n' | wc -l
478
$ make print-test query='requires|requires' queryval='ctetgen,triangle' | tr ' ' '\n' | wc -l
22
```

The total number of tests for running only ctetgen or triangle is 500. They have 22 tests in common, and 478 that
run independently of each other.

The union and intersection have fixed grouping. So this string argument

```none
query='requires,requires|args' queryval='cuda,hip,*log*'
# equivalently
query='requires%AND%requires%OR%args' queryval='cuda%AND%hip%AND%*log*'
```

will can be read as

```none
requires:cuda && (requires:hip || args:*log*)
```

which is probably not what is intended.

`query/queryval` also support negation (`!`, alternate `%NEG%`), but is limited.
The negation only applies to tests that have a related field in it. So for example, the
arguments of

```console
query=requires queryval='!cuda'
# equivalently
query=requires queryval='%NEG%cuda'
```

will only match if they explicitly have:

```
requires: !cuda
```

It does not match all cases that do not require cuda.

### Debugging for loops

One of the more difficult issues is how to debug for loops when a subset of the
arguments are the ones that cause a code crash. The default naming scheme is
not always helpful for figuring out the argument combination.

For example:

```console
$ make test s='src/ksp/ksp/tests/ex9.c' i='*1'
Using MAKEFLAGS: i=*1 s=src/ksp/ksp/tests/ex9.c
        TEST arch-osx-pkgs-opt-new/tests/counts/ksp_ksp_tests-ex9_1.counts
 ok ksp_ksp_tests-ex9_1+pc_fieldsplit_diag_use_amat-0_pc_fieldsplit_diag_use_amat-0_pc_fieldsplit_type-additive
 not ok diff-ksp_ksp_tests-ex9_1+pc_fieldsplit_diag_use_amat-0_pc_fieldsplit_diag_use_amat-0_pc_fieldsplit_type-additive
 ok ksp_ksp_tests-ex9_1+pc_fieldsplit_diag_use_amat-0_pc_fieldsplit_diag_use_amat-0_pc_fieldsplit_type-multiplicative
 ...
```

In this case, the trick is to use the verbose option, `V=1` (or for the shell script workflows, `-v`) to have it show the commands:

```console
$ make test s='src/ksp/ksp/tests/ex9.c' i='*1' V=1
Using MAKEFLAGS: V=1 i=*1 s=src/ksp/ksp/tests/ex9.c
arch-osx-pkgs-opt-new/tests/ksp/ksp/tests/runex9_1.sh  -v
 ok ksp_ksp_tests-ex9_1+pc_fieldsplit_diag_use_amat-0_pc_fieldsplit_diag_use_amat-0_pc_fieldsplit_type-additive # mpiexec  -n 1 ../ex9 -ksp_converged_reason -ksp_error_if_not_converged  -pc_fieldsplit_diag_use_amat 0 -pc_fieldsplit_diag_use_amat 0 -pc_fieldsplit_type additive > ex9_1.tmp 2> runex9_1.err
...
```

This can still be hard to read and pick out what you want. So use the fact that you want `not ok`
combined with the fact that `#` is the delimiter:

```console
$ make test s='src/ksp/ksp/tests/ex9.c' i='*1' v=1 | grep 'not ok' | cut -d# -f2
mpiexec  -n 1 ../ex9 -ksp_converged_reason -ksp_error_if_not_converged  -pc_fieldsplit_diag_use_amat 0 -pc_fieldsplit_diag_use_amat 0 -pc_fieldsplit_type multiplicative > ex9_1.tmp 2> runex9_1.err
```

## PETSC Test Harness

The goals of the PETSc test harness are threefold.

1. Provide standard output used by other testing tools
2. Be as lightweight as possible and easily fit within the PETSc build chain
3. Provide information on all tests, even those that are not built or run because they do not meet the configuration requirements

Before understanding the test harness, you should first understand the
desired requirements for reporting and logging.

### Testing the Parsing

After inserting the language into the file, you can test the parsing by
executing

A dictionary will be pretty-printed. From this dictionary printout, any
problems in the parsing are is usually obvious. This python file is used
by

in generating the test harness.

## Test Output Standards: TAP

The PETSc test system is designed to be compliant with the [Test Anything Protocol (TAP)](https://testanything.org/tap-specification.html).

This is a simple standard designed to allow testing tools to work
together easily. There are libraries to enable the output to be used
easily, including sharness, which is used by the Git team. However, the
simplicity of the PETSc tests and TAP specification means that we use
our own simple harness given by a single shell script that each file
sources: `$PETSC_DIR/config/petsc_harness.sh`.

As an example, consider this test input:

```yaml
test:
  suffix: 2
  output_file: output/ex138.out
  args: -f ${DATAFILESPATH}/matrices/small -mat_type {{aij baij sbaij}} -matload_block_size {{2 3}}
  requires: datafilespath
```

A sample output from this would be:

```
ok 1 In mat...tests: "./ex138 -f ${DATAFILESPATH}/matrices/small -mat_type aij -matload_block_size 2"
ok 2 In mat...tests: "Diff of ./ex138 -f ${DATAFILESPATH}/matrices/small -mat_type aij -matload_block_size 2"
ok 3 In mat...tests: "./ex138 -f ${DATAFILESPATH}/matrices/small -mat_type aij -matload_block_size 3"
ok 4 In mat...tests: "Diff of ./ex138 -f ${DATAFILESPATH}/matrices/small -mat_type aij -matload_block_size 3"
ok 5 In mat...tests: "./ex138 -f ${DATAFILESPATH}/matrices/small -mat_type baij -matload_block_size 2"
ok 6 In mat...tests: "Diff of ./ex138 -f ${DATAFILESPATH}/matrices/small -mat_type baij -matload_block_size 2"
...

ok 11 In mat...tests: "./ex138 -f ${DATAFILESPATH}/matrices/small -mat_type saij -matload_block_size 2"
ok 12 In mat...tests: "Diff of ./ex138 -f ${DATAFILESPATH}/matrices/small -mat_type aij -matload_block_size 2"
```

## Test Harness Implementation

Most of the requirements for being TAP-compliant lie in the shell
scripts, so we focus on that description.

A sample shell script is given the following.

```sh
#!/bin/sh
. petsc_harness.sh

petsc_testrun ./ex1 ex1.tmp ex1.err
petsc_testrun 'diff ex1.tmp output/ex1.out' diff-ex1.tmp diff-ex1.err

petsc_testend
```

`petsc_harness.sh` is a small shell script that provides the logging and reporting
functions `petsc_testrun` and `petsc_testend`.

A small sample of the output from the test harness is as follows.

```none
ok 1 ./ex1
ok 2 diff ex1.tmp output/ex1.out
not ok 4 ./ex2
#   ex2: Error: cannot read file
not ok 5 diff ex2.tmp output/ex2.out
ok 7 ./ex3 -f /matrices/small -mat_type aij -matload_block_size 2
ok 8 diff ex3.tmp output/ex3.out
ok 9 ./ex3 -f /matrices/small -mat_type aij -matload_block_size 3
ok 10 diff ex3.tmp output/ex3.out
ok 11 ./ex3 -f /matrices/small -mat_type baij -matload_block_size 2
ok 12 diff ex3.tmp output/ex3.out
ok 13 ./ex3 -f /matrices/small -mat_type baij -matload_block_size 3
ok 14 diff ex3.tmp output/ex3.out
ok 15 ./ex3 -f /matrices/small -mat_type sbaij -matload_block_size 2
ok 16 diff ex3.tmp output/ex3.out
ok 17 ./ex3 -f /matrices/small -mat_type sbaij -matload_block_size 3
ok 18 diff ex3.tmp output/ex3.out
# FAILED   4 5
# failed 2/16 tests; 87.500% ok
```

For developers, modifying the lines that get written to the file can be
done by modifying `$PETSC_DIR/config/example_template.py`.

To modify the test harness, you can modify `$PETSC_DIR/config/petsc_harness.sh`.

### Additional Tips

To rerun just the reporting use

```console
$ config/report_tests.py
```

To see the full options use

```console
$ config/report_tests.py -h
```

To see the full timing information for the five most expensive tests use

```console
$ config/report_tests.py -t 5
```
