// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Samples {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;

    @EntryPoint()
    operation SampleRandomNumber(nQubits : Int) : Result[] {

        // We prepare a register of qubits in a uniform
        // superposition state, such that when we measure,
        // all bitstrings occur with equal probability.
        use register = Qubit[nQubits];

        // Set qubits in superposition.
        ApplyToEachA(H, register);

        // Measure all qubits and return.
        return ForEach(MResetZ, register);

    }
}
