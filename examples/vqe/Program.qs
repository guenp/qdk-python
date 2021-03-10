namespace Microsoft.Quantum.Samples.Chemistry.SimpleVQE {
    open Microsoft.Quantum.Core;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.JordanWigner;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.StatePrep;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.Utils;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.ChemUtils;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.EstimateEnergy;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Characterization;

    @EntryPoint()
    // operation Test(theta1: Double, theta2: Double, theta3: Double, nSamples: Int) : Double {
    // operation Test(theta1: Double, theta2: Double, theta3: Double, termType: Int, nHam: Int, nOp: Int) : Result {
    operation Test(theta1: Double, theta2: Double, theta3: Double, nSamples: Int) : Double {
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
                JordanWignerInputState((theta3, 0.0), [2, 3, 1, 0]),
                JordanWignerInputState((1.0, 0.0), [0, 1])
            ]
        );
        let (ZData, ZZData, PQandPQQRData, h0123Data) = hamiltonian!;
        let hamiltonianTermArray = [ZData, ZZData, PQandPQQRData, h0123Data];
        let energyOffset = -0.09883444600000002;
        mutable energy = 0.0;
        // let hamiltonianTerms = hamiltonianTermArray[termType];
        // // Message($"Total ham terms: {Length(hamiltonianTerms)}");
        // let hamiltonianTerm = hamiltonianTerms[nHam];
        // let (qubitIndices, coefficient) = hamiltonianTerm!;
        // let measOps = VQEMeasurementOperators(nQubits, qubitIndices, termType);
        // Message($"Total measurement operations: {Length(measOps)}");
        // let op = measOps[nOp];
        // use register = Qubit[nQubits];
        // PrepareTrialState(inputState, register);
        // let result = Measure(op, register);
        // ResetAll(register);
        // return result;

        mutable nResults = 0;
        mutable nCoeffs = 0;
        for termType in 0..Length(hamiltonianTermArray)-1 {
            let hamiltonianTerms = hamiltonianTermArray[termType];
            for hamiltonianTerm in hamiltonianTerms {
                set nResults += Length(hamiltonianTerms);
                let (qubitIndices, coefficient) = hamiltonianTerm!;
                let coefficients = ExpandedCoefficients(coefficient, termType);
                for coeff in coefficients {
                    if (coeff >= 1e-10 or coeff <= -1e-10) {
                        set nCoeffs += 1;
                    }
                }
            }
        }
        mutable results = new Result[nCoeffs];
        mutable coeffs = new Double[nCoeffs];
        let nMeasurements = nSamples;
        for termType in 0..Length(hamiltonianTermArray)-1 {
            let hamiltonianTerms = hamiltonianTermArray[termType];
            // Message($"Total ham terms: {Length(hamiltonianTerms)}");
            for hamiltonianTerm in hamiltonianTerms {
                let (qubitIndices, coefficient) = hamiltonianTerm!;
                let coefficients = ExpandedCoefficients(coefficient, termType);
                let measOps = VQEMeasurementOperators(nQubits, qubitIndices, termType);
                mutable i = 0;
                mutable j = 0;
                for op in measOps {
                    let coeff = coefficients[i];
                    set i += 1;
                    if (coeff >= 1e-10 or coeff <= -1e-10) {
                        // use register = Qubit[nQubits];
                        // PrepareTrialState(inputState, register);
                        // let result = Measure(op, register);
                        let energy_ = EstimateFrequencyA(_prepareTrialStateWrapper(inputState, _), Measure(op, _), nQubits, nMeasurements);
                        // ResetAll(register);
                        // set results w/= j <- result;
                        // set coeffs w/= j <- coeff;
                        // set j += 1;
                        set energy += (2. * energy_ - 1.) * coeff;
                    }
                }
            }
        }

        return energy + -0.09883444600000002;

        // Old code
        // for termType_ in 0..Length(hamiltonianTermArray)-1 {
        //     let hamiltonianTerms = hamiltonianTermArray[termType_];
        //     for hamiltonianTerm in hamiltonianTerms {
        //         let (qubitIndices, coefficient) = hamiltonianTerm!;
        //         let measOps = VQEMeasurementOperators(nQubits, qubitIndices, termType_);
        //         let coefficients = ExpandedCoefficients(coefficient, termType_);
        //         let jwTermEnergy = SumTermExpectation(inputState, measOps, coefficients, nQubits, nSamples);
        //         set energy += jwTermEnergy;
        //     }
        // }
        // return energy;
    }
}
