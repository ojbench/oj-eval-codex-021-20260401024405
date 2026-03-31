# Problem 021 - Rectangle Counting

**ACMOJ Problem ID**: 2981

## Table of Contents

- [Problem 021 - Rectangle Counting](#problem-021-rectangle-counting)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
    - [Background](#background)
  - [Assignment Description](#assignment-description)
    - [Grade Composition](#grade-composition)
  - [Assignment Requirements](#assignment-requirements)
    - [Input Format](#input-format)
    - [Output Format](#output-format)
    - [Samples](#samples)
    - [Data Constraints](#data-constraints)
  - [Per-Testcase Resource Limits](#per-testcase-resource-limits)
  - [Test Data](#test-data)
  - [Submission Requirements](#submission-requirements)
    - [OJ Git Repository Compilation Process](#oj-git-repository-compilation-process)
    - [Git Configuration Requirements](#git-configuration-requirements)
    - [Submission Guidelines](#submission-guidelines)
    - [Evaluation Notes](#evaluation-notes)
    - [Academic Integrity](#academic-integrity)

## Introduction

### Background

Counting geometric objects under constraints is a fundamental class of problems in computational geometry. In this problem, you are asked to count the number of axis-aligned empty rectangles whose corners lie on a given set of marked points. The challenge lies in efficiently determining whether each candidate rectangle is "empty" (contains no other marked points), which requires combining geometric reasoning with efficient algorithmic techniques such as divide-and-conquer or sweep-line algorithms.

## Assignment Description

### Grade Composition

| Grading Component | Percentage |
| :--: | :--: |
| Pass **2981. Rectangle Counting** | 80% |
| Code Review | 20% |

Here are several points that need clarification:

- In the Code Review, we will **strictly examine your code style and repository organization structure, etc.**

- This assignment provides some sample data for testing, stored in the `/workspace/data/021/data_test/` directory. Note that these are not the test cases on the Online Judge. Passing all local test cases does not guarantee that you will pass the OJ tests.

- Besides the provided sample data, we also encourage you to design your own test data based on your program logic to assist debugging.

## Assignment Requirements

### Problem Description

Given $n$ marked points in a 2D plane, count the number of rectangles that satisfy all of the following conditions:

- The sides of the rectangle are parallel to the coordinate axes.
- The bottom-left corner and the top-right corner of the rectangle are marked points.
- There are no marked points on the interior or the boundary of the rectangle (except for the bottom-left and top-right corners).

### Input Format

The first line contains an integer $n$, the number of marked points.

The following $n$ lines each contain two integers $x, y$, representing the coordinates of a marked point.

### Output Format

Output a single integer: the number of rectangles that satisfy the above conditions.

### Samples

#### Sample 1

**Input:**
```
4
0 0
2 2
3 4
4 3
```

**Output:**
```
3
```

#### Sample 2

**Input:**
```
10
2 1
3 0
6 3
10 2
16 4
0 8
8 12
11 14
14 11
18 10
```

**Output:**
```
15
```

### Data Constraints

- Subtask 1 (20 pts): $n \leq 400$
- Subtask 2 (30 pts): $n \leq 5000$
- Subtask 3 (50 pts): $n \leq 2 \times 10^5$

All data guarantees: $2 \leq n \leq 2 \times 10^5$, $0 \leq x, y \leq 10^9$, all $x$-coordinates are distinct, all $y$-coordinates are distinct.

## Per-Testcase Resource Limits

- **Time Limit (per test case)**: 1000 ms
- **Memory Limit (per test case)**: 256 MiB
- **Disk Usage**: No disk usage is permitted.

## Test Data

The test data for this problem is located at `/workspace/data/021/data_test/`.

## Submission Requirements

### OJ Git Repository Compilation Process

For Git compilation, we will first clone the repository using a command similar to:
```bash
git clone <repo_url> . --depth 1 --recurse-submodules --shallow-submodules --no-local
```

Then we check if there is a `CMakeLists.txt` file. If it exists, we run (if not, a warning message will be displayed):
```bash
cmake .
```

Finally, we check if there is any of `GNUmakefile`/`makefile`/`Makefile` (if cmake was run previously, this will be the generated Makefile). If it exists, we run (if not, a warning message will be displayed):
```bash
make
```

After this process is complete, we will use the `code` file in the project root directory as the compilation result.

A `CMakeLists.txt` file is provided in the project. You can use or modify it as needed. The local environment has gcc-13 and g++-13 available.

### Git Configuration Requirements

**IMPORTANT**: You must create a `.gitignore` file in your project root directory to avoid OJ evaluation conflicts.

The `.gitignore` file should include at least the following entries:

```gitignore
CMakeFiles/
CMakeCache.txt
```

### Submission Guidelines

- The submitted code must be able to compile successfully through the above compilation process
- The compiled executable file name must be `code`
- The program needs to be able to read data from standard input and write results to standard output
- Please ensure the code runs correctly within the given time and space limits
- **You must use C++ or C language** to implement this assignment

### Evaluation Notes

- The evaluation system will test your program using the provided test data
- The program output must exactly match the expected output (including format)
- Exceeding time or memory limits will be judged as the corresponding error type
- Please pay attention to the overall time performance of your code and the time complexity of each part of your algorithm.

### Academic Integrity

If any violations are found during evaluation or code review (including but not limited to using unconventional methods to pass test cases), your final score may be significantly reduced or become **0 points**.
