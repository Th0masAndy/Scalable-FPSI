# Towards Scalable Fuzzy PSI via Efficient Fuzzy Matching

This repository provides a research implementation of fuzzy private set intersection (FPSI) protocols.

> Note: This project is experimental and primarily intended for research use. Parameters should be chosen according to the available memory, CPU features, and intended benchmark size.


## Location of Main Functionality
- `src/main.cpp` parses command-line flags and dispatches to the selected protocol variant.
- `src/fpsiLowPx.cpp` and `src/fpsiHighPx.cpp` implement the main FPSI variants.
- `src/bp25Px.cpp` and `src/bp25High.cpp` implement baseline variants used for comparison.
- `src/opprf.cpp`, `src/cmp.cpp`, `src/eq.cpp`, `src/mul.cpp`, `src/mux.cpp`, and `src/permute.cpp` contain reusable building blocks such as OPPRF, comparison, equality test, multiplication, mux, and permutation routines.
- `include/*.h` contains the corresponding public declarations, shared parameters, and utility functions.
- `test/opprf_test.cpp` contains standalone test helpers for OPPRF and comparison components; it is not currently wired into `CMakeLists.txt`.

## Requirements

- Linux on **AMD64/x86_64**. The CMake target enables AES/PCLMUL/SSE flags on x86 systems.
- A C++20 compiler, CMake, Make, Git, and Python 3.
- Third-party libraries installed under `thirdparty/out/install` by [build.sh](./build.sh):
  - [volePSI](thirdparty/volepsi), including cryptoTools, libOTe, coproto, macoro, and sparsehash dependencies fetched by volePSI's build. (We slightly modified the source code of RsOpprf and RsOprf to enable large inputs)
  - [BLAKE3](https://github.com/BLAKE3-team/BLAKE3).

- **System packages commonly needed for local builds:**

```bash
build-essential
cmake
git
python3
nasm
libssl-dev
libgmp-dev
libtool
pkg-config
```

## Local build

From the project root directory:

```bash
mkdir -p thirdparty
./build.sh

cmake -S . -B build
cmake --build build -j

# The executable will be located at ./build/fpsi
```

`build.sh` clones fixed revisions of volePSI and BLAKE3, installs them into `thirdparty/out/install`, and removes the temporary source checkouts afterwards.

## Docker (optional)

Use the provided Dockerfile for an isolated build environment:

```bash
docker build -t fpsi .
docker run -it --name fpsi --cap-add=NET_ADMIN --memory=512g fpsi
```

For large benchmarks, allocate enough memory and enable CPU features required by the third-party cryptographic libraries.

## Command-line Options

The first flag selects the protocol variant. Flags use a leading dash, for example `-high` or `-nn 10`.

| Flag | Meaning | Values / Notes |
|---|---|---|
| `-low` | Run the low-dimensional/main FPSI framework | Dispatches to `fpsiLowLpPx` in the current code |
| `-high` | Run the high-dimensional FPSI framework | Uses `fpsiHighPx` when `-p 0`, otherwise `fpsiHighLpPx` |
| `-bp25low` | Run the low-dimensional BP25-style baseline | Uses `bp25LowPx` when `-p 0`, otherwise `bp25LowLpPx` |
| `-bp25high` | Run the high-dimensional BP25-style baseline | Uses `bp25High` when `-p 0`, otherwise `bp25HighLp` |
| `-n` | Input set size | Overrides `-nn`; default is `1 << nn` |
| `-nn` | log2 of input set size | default `10` |
| `-d` | Dimension | default `2` |
| `-delta` | Distance threshold | default `2` |
| `-p` | Distance mode for LP variants | `0`: L-infinity-style path, `1`: L1, `2`: L2 where supported |
| `-try` | Number of repeated runs | default `1` |
| `-v` | Verbose timing / match output | `0`: off, `1`: on |

The program generates synthetic sender and receiver inputs internally. It does not currently read datasets from files.

## Usage Examples

Build first, then run from the build directory or use `./build/fpsi` from the repository root.

Run the high-dimensional default variant:

```bash
./build/fpsi -high -nn 10 -d 2 -delta 2 -try 1
```

Run the high-dimensional L2 path:

```bash
./build/fpsi -high -p 2 -nn 10 -d 2 -delta 2 -try 1
```

Run the low-dimensional framework with verbose output:

```bash
./build/fpsi -low -p 2 -nn 8 -d 2 -delta 2 -v 1
```

Run baseline variants:

```bash
./build/fpsi -bp25low -nn 8 -d 2 -delta 2
./build/fpsi -bp25high -p 2 -nn 8 -d 2 -delta 2
```

------------------------------------------------------------------------

## Baseline Implementations

The repository includes [BP25](https://eprint.iacr.org/2025/911) baseline in `src/bp25Px.cpp` and `src/bp25High.cpp`.

The external artifact used for comparison: 

### Piske et al. 

[Code](https://github.com/asu-crypto/daOT-fuzzyPSI) | [Paper](https://eprint.iacr.org/2025/996)

------------------------------------------------------------------------

## Acknowledgements

This project builds on the following open-source libraries and research codebases:

- [volePSI](https://github.com/ladnir/volepsi), libOTe, cryptoTools, coproto, and macoro for cryptographic protocol primitives.
- [BLAKE3](https://github.com/BLAKE3-team/BLAKE3) for hashing.
- Components inspired by [OpenCheetah](https://github.com/Alibaba-Gemini-Lab/OpenCheetah) millionaire comparison code.


## Citation

If you make use of our work, please consider citing us:

```bibtex
@inproceedings{hao2026scalable,
    author = {Hao, Meng and Yang, Xinpeng and Chen, Hanxiao and Zhang, Tianwei and Xue, Haiyang and Yang, Guomin and Li, Hongwei and Deng, Robert H.},
    title = {{Towards Scalable Fuzzy PSI via Efficient Fuzzy Matching}},
    booktitle = {Proceedings of the 2026 ACM SIGSAC Conference on Computer and Communications Security (CCS)},
    year = {2026},
}
```
