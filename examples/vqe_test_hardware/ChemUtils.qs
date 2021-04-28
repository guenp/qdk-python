// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Samples.Chemistry.SimpleVQE.ChemUtils {

    
    /// # Summary
    /// Format of data passed from C# to Q# to represent a term of the Hamiltonian.
    // The list of integers indicates the hamiltonian term indices and the list
    // of doubles contains the coefficient.
    newtype HTerm = (Int[], Double[]);
    
}


