namespace Microsoft.Quantum.Samples.Chemistry.SimpleVQE {
    open Microsoft.Quantum.Core;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.JordanWigner;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.StatePrep;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.Utils;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.ChemUtils;
    open Microsoft.Quantum.Arrays;

    @EntryPoint()
    operation Test(theta1: Double, theta2: Double, theta3: Double, termType: Int, nHam: Int, nOp: Int) : Result {
    // operation Test() : Result {
    //     let theta1 = 0.001;
    //     let theta2 = -0.001;
    //     let theta3 = 0.001;
    //     let termType = 0;
    //     let nHam = 0;
    //     let nOp = 0;
        let nQubits = 4;
        let hamiltonian = JWOptimizedHTerms(
        [
            HTerm([0], [0.17120128499999998]),
            HTerm([1], [0.17120128499999998]),
            HTerm([2], [-0.222796536]),
            HTerm([3], [-0.222796536])],
        [
            HTerm([0, 1], [0.1686232915]),
            HTerm([0, 2], [0.12054614575]),
            HTerm([0, 3], [0.16586802525]),
            HTerm([1, 2], [0.16586802525]),
            HTerm([1, 3], [0.12054614575]),
            HTerm([2, 3], [0.1743495025])
        ],
        new HTerm[0],
        [
            HTerm([0, 1, 2, 3], [0.0, -0.0453218795, 0.0, 0.0453218795])
        ]
        );
        let inputState = (
            2,
            [
                JordanWignerInputState((theta1, 0.0), [2, 0]),
                JordanWignerInputState((theta2, 0.0), [3, 1]),
                JordanWignerInputState((theta3, 0.0), [2, 3, 1, 0])
                // JordanWignerInputState((1.0, 0.0), [0, 1])
            ]
        );
        let (ZData, ZZData, PQandPQQRData, h0123Data) = hamiltonian!;
        let hamiltonianTermArray = [ZData, ZZData, PQandPQQRData, h0123Data];
        let hamiltonianTerms = hamiltonianTermArray[termType];
        // Message($"Total ham terms: {Length(hamiltonianTerms)}");
        let hamiltonianTerm = hamiltonianTerms[nHam];
        let (qubitIndices, coefficient) = hamiltonianTerm!;
        let measOps = VQEMeasurementOperators(nQubits, qubitIndices, termType);
        // Message($"Total measurement operations: {Length(measOps)}");
        let op = measOps[nOp];
        use register = Qubit[nQubits];
        PrepareTrialStateUCCSD(inputState, register);
        let result = Measure(op, register);
        ResetAll(register);
        return result;
    }
}