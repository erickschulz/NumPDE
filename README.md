# NPDECODES [![pipeline status](https://gitlab.math.ethz.ch/ralfh/NPDERepo/badges/master/pipeline.svg)](https://gitlab.math.ethz.ch/ralfh/NPDERepo/-/pipelines)
This repository contains the codes for the homework problems of the recurring course **Numerical Methods for Partial Differential Equations** at [ETH Zurich](https://ethz.ch/en.html). The course treats finite element methods (FEMs) using the C++ library [LehrFEM++](https://github.com/craffael/lehrfempp) and relies on the following material:
* [lecture notes](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
* [homework problems](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/HOMEWORK/NPDEFL_Problems.pdf)

Moreover, enrolled students can access the [moodle](https://moodle-app2.let.ethz.ch/course/view.php?id=12060) page of the course.

## Requirements
Currently, only UNIX based operating systems are supported. Moreover, you need to have the following installed on your machine:
* C++20 compiler (e.g. gcc, clang)
* CMake (at least VERSION 3.10)
* python3
* A reader for .vtk files (e.g. paraview)
* git (not strictly needed, you could also download the repo as .zip file)

## Getting started
This section is suited only for your own computer. To build the codes on the student computers of ETH see below. Open a terminal and type
```
git clone https://gitlab.math.ethz.ch/ralfh/NPDERepo.git
cd NPDERepo/
mkdir build
cd build/
cmake ..
```
This will install LehrFEM++ and its dependencies into a folder `~/.hunter/`.

**Note**: The first compilation can take a significant amount of time. It is highly recommended to use parallel compilation with the `-j` flag. For example, to use 8 parallel jobs:
```
make -j8
```
This will significantly speed up the build process.

### Build Configurations

For optimal performance, always specify the build type during the `cmake` configuration step.

*   **Release Build (Optimized & Fast):**
    ```bash
    # Inside your build directory
    cmake -DCMAKE_BUILD_TYPE=Release ..
    ```
    This is recommended for running computationally intensive problems and matches the CI pipeline's configuration.

*   **Debug Build (For Development & Debugging):**
    ```bash
    # Inside your build directory
    cmake -DCMAKE_BUILD_TYPE=Debug ..
    ```
    This builds with debug symbols and no optimizations, making it suitable for use with a debugger.

To build a specific problem, say `TestQuadratureRules`, proceed as follows:
```
cd homeworks/TestQuadratureRules/
make
```
This will build from the source files in `NPDECRepo/homeworks/TestQuadratureRules/`, where the subfolder `mysolution/` contains templates to be changed by the students. Recompilation is done by invoking `make` again. The following executables are generated:
* `./TestQuadratureRules_mastersolution`: Runs the mastersolution.
* `./TestQuadratureRules_test_mastersolution`: Runs unit tests on all important functions of the mastersolution.
* `./TestQuadratureRules_mysolution`: Runs the students code, i.e. the one in `mysolution/`.
* `./TestQuadratureRules_test_mysolution`: Runs unit tests the students code, i.e. the one in `mysolution/`.

There is two folders called `homeworks/`. One contains the source files and one contains the executables:
```
.
├── build (was created by you)
│   ├── homeworks
│   :   ├── TestQuadratureRules
│       :   ├── TestQuadratureRules_mastersolution      (executable)
│           ├── TestQuadratureRules_mysolution          (executable)
│           ├── TestQuadratureRules_test_mastersolution (executable)
│           ├── TestQuadratureRules_test_mysolution     (executable)
│           :
│
├── homeworks
:   ├── TestQuadratureRules
    :   ├── mastersolution (folder containing source files)
        ├── mysolution     (folder containing source files, to be modified by you)
        ├── templates      (folder containing source files)
        :
```

## On the student computers
LehrFEM++ is already installed on the linux student computers in the ETH main building. You can even use them remotely by typing in your terminal
```
ssh -X <nethz_username>@slab1.ethz.ch
```
where `<nethz_username>` has the be replaced by your ETH username. To set up your local repository on the student computers, type:
```
cd /tmp
git clone https://gitlab.math.ethz.ch/ralfh/NPDERepo.git
mv NPDERepo ~
cd ~/NPDERepo
mkdir build
cd build
export HUNTER_ROOT=/opt/libs/NumPDE
cmake ..
```
The first four lines are due to limited resources on the student computers. Setting the environment variable `HUNTER_ROOT` tells CMake where to look for the preinstalled libraries. This environment variable is local to your terminal, i.e. has to be redefined if you start a new terminal. Apart from this, you can use the folder `~/NPDERepo` in the same way you would for the approach in the previous section. However, you have only very little memory available on the student computers. We therefore recommend to only build one problem at a time.

## FAQ

### clang: error: unknown argument: '-fcoalesce-templates'
Mac users, after updating to macOS Catalina 10.15.4, are receiving this error. The workaround is as follows: Navigate  your terminal into the folder `NPDECODES/` and type:
```
brew install gcc@8    #install gcc version 8, needs brew to be installed first...
gcc-8 --version       #check if gcc-8 was installed properly
g++-8 --version       #check if g++-8 was installed properly
rm -rf build          #delete the old build folder
mkdir build           #recreate it
```
If this has succeeded, you need to build the codes using the gcc compiler by defining the environment variables `CC` and `CXX`. This is done by navigating a terminal into `NPDECODES/build/` and running:
```
export CC=gcc-8
export CXX=g++-8
cmake ..
```
If the installation is successful, you can than build your codes using `make` as before. Note that the gcc version under OSX usulally just links to clang. However, the procedure above installs the actual gcc compiler.

# Notes for Developers
## Creating a new Homework
New homeworks should be created in the `developers` folder. When the homework is ready to be deployed,
it can be done by running `scripts/deploy_npde.py`.
Please adhere to the directory layout:
```
ProblemName/
└── mastersolution/
│   ├── scripts/            (scripts, bash or python, that are part of the mastersolution)
│   │   └── ...
│   ├── problemname_main.cc (contains main function, SHOULD NOT #include "problemname.cc", only problemname.h!)
│   ├── problemname.cc      (further source code if needed)
│   └── problemname.h       (header for classes and other declarations)
├── meshes/
│   └── ...                 (any mesh files if applicable)
├── scripts/                (scripts, bash or python, that are for use by students solving the homework)
│   └── ...
├── CMakeLists.txt
└── README.md
```
You can create a template using `scripts/python/new_problem.py <Problemname> <problemname>`. Feel free to delete any folders you don't need.

### Build System and Dependencies
The build system automatically discovers source files (`.cc` and `.h`) in the problem directory and provides default libraries (`Eigen3::Eigen` and `LF_ALL` for main sources, `GTest::gtest_main` for tests). In most cases, you don't need to specify dependencies manually.

However, if your problem requires:
- Additional libraries (e.g., `Boost::program_options` in `LeapfrogDissipativeWave`)
- Source files from subdirectories (e.g., `meshes/mesh.cc` in `ElementMatrixComputation`)
- Exclusion of certain files from compilation (e.g., `GaussLobattoParabolic`)

you can create a `dependencies.cmake` file in your `mastersolution/` directory with:
```cmake
set(SOURCES
  ${DIR}/problemname_main.cc
  ${DIR}/problemname.cc
  meshes/mesh.cc
)
set(LIBRARIES Eigen3::Eigen LF_ALL Boost::program_options)
```
Similarly, for tests, create `mastersolution/test/dependencies.cmake` if needed.

### Solution Tags

In the files of `./developers/mastersolution/` we put the following tags
```
#if SOLUTION
  <only in mastersolution>
#else
  <only in template>
#endif
```
to indicate what belongs to mastersolution and/or template. Based on these tags, the file `./scripts/deploy_npde.py` generates a directory `./homeworks/<ProblemName>/` containg the directories `mastersolution`, `mysolution`, `templates` with the corresponding content. The students work exclusively in `./homeworks/<ProblemName>/`.

### Working with and scripts meshes in source code
`CMake` creates dynamic links to `meshes` and `scripts` folders. Below are example usages.

#### Accessing a mesh file:
```cpp
std::string mesh_file = "meshes/hexagonal.msh";
```
#### Calling a script that is part of the mastersolution (i.e. in `ProblemName/mastersolution/scripts`):
```cpp
#include "systemcall.h"
...
systemcall::execute(
      "python3 ms_scripts/viswave.py solution.csv solution.eps");
```

#### Calling a script that is not part of the mastersolution (i.e. in `ProblemName/scripts`):
```cpp
#include "systemcall.h"
...
systemcall::execute(
      "python3 scripts/plot.py solution.csv");
```
The function `systemcall::execute` is a wrapper of `std::system` that checks whether the system call was successful.
The pipeline forbids `std::system` calls for this reason!

**NOTE: for these relative symbolic paths to work, you must call the executable from the problem folder in `build`!**

## Continuous Integration
For every commit, a pipeline will:
- check whether the code still adheres to the formatting standard set by clang-format (mostly follows google standard)
- compile the `developers` and `lecturecodes` folder on both mac and linux
- run a static analysis on `developers` and `lecturecodes` folder (clang-tidy)
- run all mastersolution executables and tests on both mac an linux

### clang-format
Before pushing, you can automatically apply the formatting standard by running `scripts/apply-clang-format.sh` *from the root of the repository*.

**NOTE**:
> Different versions of `clang-format` will result in different formattings given the same configuration (in file `.clang-format`)
> Because of this **the pipeline may fail if you use a different version than the gitlab runners**.
> The easiest way to avoid this is by using [pre-commit](https://pre-commit.com) with the supplied `.pre-commit-config.yaml`.

#### Installing pre-commit hook for automatic formatting
[pre-commit](https://pre-commit.com) is a lightweight python package. Install it with `pip`:
```bash
pip install pre-commit
```
Afterwards, you should be able to install the configured git hook by running
```bash
pre-commit install
```
from the **root of the repository**. This also installs the current version of `clang-format` used by this repo.

From now on *whenever you commit new code*, the changes will be checked for formatting.
If there are format issues, *will cancel* and the the files will be reformatted.
If the automatic format didn't break anything, you can just redo the commit.

### clang-tidy
The list of checks done by `clang-tidy` can be seen in the `.clang-tidy` config file. We treat warnings as errors so expect the pipeline to fail when pushing code for the first time.
The list of checks is a work in progress, so feel free to suggest adding/removing one (e.g. if a check causes a lot of false positives or seems overly strict).

#### In case of False Positives
Sometimes, `clang-tidy` is wrong. A prominent example is marking variables as unused when they're only called within an `LF_ASSERT_MSG` macro function in test executables.
When you are confident that `clang-tidy` is giving a false positive, you can tell it to ignore the line(s) as follows:

```cpp
  // NOLINTBEGIN(clang-analyzer-deadcode.DeadStores)
  const lf::fe::ScalarReferenceFiniteElement<double> *rsf_edge_p =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
  // NOLINTEND(clang-analyzer-deadcode.DeadStores)
```
Note that `clang-analyzer-deadcode.DeadStores` is the check warns about unused variables. Replace this with the appropriate check you want to ignore.
For single line `NOLINT` you can also do
```cpp
int n = 1; // NOLINT(<check>)

// NOLINTNEXTLINE(<check>)
int n = 1;
```
**But note** that sometimes `clang-format` will force a single line statement onto two lines, making the single line `NOLINT` not work!