namespace Microsoft.Quantum.Samples.Chemistry.SimpleVQE {
    open Microsoft.Quantum.Core;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Characterization;

    operation RunVQE(nQubits: Int, hamiltonian: JWOptimizedHTerms, inputType: Int, inputState: JordanWignerInputState[], energyOffset: Double, theta1: Double, theta2: Double, theta3: Double, termType: Int, nHam: Int, nOp: Int) : Result {
        let (ZData, ZZData, PQandPQQRData, h0123Data) = hamiltonian!;
        let hamiltonianTermArray = [ZData, ZZData, PQandPQQRData, h0123Data];
        let hamiltonianTerms = hamiltonianTermArray[termType];
        let hamiltonianTerm = hamiltonianTerms[nHam];
        let (qubitIndices, coefficient) = hamiltonianTerm!;
        let measOps = VQEMeasurementOperators(nQubits, qubitIndices, termType);
        let op = measOps[nOp];
        use register = Qubit[nQubits];
        PrepareTrialStateUCCSD(inputState, register);
        let result = Measure(op, register);
        ResetAll(register);
        return result;
    }
}
