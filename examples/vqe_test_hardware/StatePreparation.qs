// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Samples.Chemistry.SimpleVQE.StatePrep {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Samples.Chemistry.SimpleVQE.JordanWigner;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Simulation;

    /// Mock operation to prepare a variational trial state.
    operation PrepareTrialState (stateData : (Int, JordanWignerInputState[]), qubits : Qubit[]) : Unit {
        let (stateType, terms) = stateData;
        for (term in terms) {
            let (coefficient, excitation) = term!;
            let (theta, phi) = coefficient;
            // Prepare an arbitrary state
            for (i in excitation)
            {
                X(qubits[i]);
            }
            // Apply Ry and Rz rotations to mimick a real variational state
            if (theta < 1.0) {
                for (qubit in qubits) {
                    // Do a Y rotation
                    within {
                        S(qubit);
                    }
                    apply {
                        Rx(theta, qubit);
                    }
                    // Then an X rotation
                    Rz(theta, qubit);
                }
            }
        }
    }

    function _JordanWignerClusterOperatorGeneratorIndex(data: JordanWignerInputState): ((Int[], Double[]), Int[]) {
        let ((real, imaginary), idxFermions) = data!;
        if (Length(idxFermions) == 2) {
            // PQ term
            return (([0],[real]),idxFermions);
        }
        elif(Length(idxFermions) == 4){
            // PQRS term
            return (([2],[real]),idxFermions);
        }
        else{
            // Any other term in invalid
            return (([-1],[0.0]),[0]);
        }
    }

    operation TrotterStepImpl (terms: JordanWignerInputState[], idx: Int, trotterStepSize: Double, qubits: Qubit[]) : Unit {
        let ((real, imaginary), idxFermions) = terms[idx]!;
        if (Length(idxFermions) == 2) {
            // PQ term
            let generatorIndex = (([0],[real]),idxFermions);
            _ApplyJordanWignerClusterOperatorPQTerm (GeneratorIndex(generatorIndex), trotterStepSize, qubits);
        }
        elif(Length(idxFermions) == 4){
            // PQRS term
            let generatorIndex = (([2],[real]),idxFermions);
            _ApplyJordanWignerClusterOperatorPQRSTerm (GeneratorIndex(generatorIndex), trotterStepSize, qubits);
        }
        // let generatorIndex = _JordanWignerClusterOperatorGeneratorIndex(terms[idx]);
        // let ((idxTermType, idxDoubles), idxFermions) = generatorIndex;
        // let termType = idxTermType[0];
        // if (termType == 0) {
        //     _ApplyJordanWignerClusterOperatorPQTerm (GeneratorIndex(generatorIndex), trotterStepSize, qubits);
        // } elif (termType == 2) {
        //     _ApplyJordanWignerClusterOperatorPQRSTerm (GeneratorIndex(generatorIndex), trotterStepSize, qubits);
        // }
    }

    operation PrepareTrialStateUCCSD (stateData : (Int, JordanWignerInputState[]), qubits : Qubit[]) : Unit {
        // Apply X to qubits 0 and 1
        X(qubits[0]);
        X(qubits[1]);
        let (stateType, initTerms) = stateData;
        let nTerms = Length(initTerms);
        let terms = initTerms[0..nTerms-1];
        let trotterStepSize = 1.0;
        let trotterOrder = 1;
        let maxTime = 1.0;
        for idx in 0 .. nTerms - 1 {
            TrotterStepImpl(terms, idx, trotterStepSize, qubits);
        }
    }

    /// # Summary
    /// Private wrapper around PrepareTrialState to make it compatible with EstimateFrequencyA by defining an adjoint.
    /// EstimateFrequencyA has built-in emulation feature when targeting the QuantumSimulator, which speeds up its execution.
    ///
    /// # Input
    /// ## inputState
    /// The Jordan-Wigner input required for PrepareTrialState to run.
    /// ## qubits
    /// A qubit register.
    // operation _prepareTrialStateWrapper(inputState : (Int, JordanWignerInputState[]), qubits : Qubit[]) : Unit is Adj {

    //     body (...) {
    //         Microsoft.Quantum.Diagnostics.AssertAllZero(qubits);
    //         PrepareTrialStateUCCSD(inputState, qubits);
    //     }

    //     // Define a non-matching adjoint body for compliance with EstimateFrequencyA
    //     adjoint (...) {
    //         ResetAll(qubits);
    //     }
    // }
}
