# |Ψ⟩CC — Quantum Circuits in Compiler

A fully `constexpr` quantum circuit emulator with zero hardcoded gate logic, executing quantum algorithms in **constant time** directly in your favourite modern C++ compiler.

## Vision

|Ψ⟩CC reimagines quantum computing by pushing quantum circuit simulation into the compile-time domain. By leveraging C++20's `constexpr` capabilities, I've built a quantum emulator that:

- **Compiles to constant-time execution** — Quantum algorithms run entirely during compilation with zero runtime overhead. For that reason, no measurement is possible, everything runs fully deterministic (hence no noise or multiple shots implemented), the output is an idealistic probability distribution. Of course the project is mere a toy and made for the love of maths behind quantum mechanics without any serious intention.
- **Zero hardcoded gate logic** — Gate implementations and applications are performed via a custom linear algebra engine, without classically implemented logical operations. 
- **Template-driven design** — Type-safe quantum operations with full compile-time verification, including the check for unitary matrices via concepts.
- **No external dependencies** — Pure C++ implementation with custom `constexpr` math utilities

## Getting Started

### Build

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

## Project Structure

```
include/
├── constexprmath/       # Compile-time math utilities
├── engine/              # Core circuit and gate infrastructure
└── gates/               # Quantum gate implementations

src/
├── bell_state.cpp       # Bell state example
├── ghz_state.cpp        # GHZ state example
└── shor_21_2.cpp        # Shor's algorithm example (21 = 3×7)
```
